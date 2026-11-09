#!/usr/bin/env python3
"""Tasmanian bushfire impacts from Sentinel-2.

Pick a recent Tasmanian bushfire (a built-in site, a fire from the Tasmanian Fire History layer, or
any longitude/latitude) and this script:

1. fetches the official fire boundary and dates from theLIST Fire History layer (State of Tasmania);
2. builds cloud-free median Sentinel-2 composites from Digital Earth Australia (DEA) Analysis Ready
   Data for: before the fire, just after it, and each summer since;
3. maps burn severity with the differenced Normalised Burn Ratio (dNBR), corrected with nearby
   unburnt land, and checks it against the official boundary;
4. tracks recovery each summer by severity class, against unburnt land;
5. writes figures and a summary table to outputs/tasmania/<site>/.

Examples
    python tas_fires.py --list
    python tas_fires.py --site riveaux-road
    python tas_fires.py --search "pine tier"
    python tas_fires.py --fire-id 10137
    python tas_fires.py --lon 146.7 --lat -43.1
    python tas_fires.py --all            # every built-in site, then the overview
    python tas_fires.py --overview       # overview + interactive map from finished sites

Downloads (fire boundaries, composites) are cached in downloads/tasmania/ (not committed).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import math
import os
import re
import sys
import time
import warnings
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# GDAL settings for reading cloud-optimised GeoTIFFs over HTTPS. Set as environment variables so
# they apply in every worker thread.
for _k, _v in {
    "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
    "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif",
    "GDAL_HTTP_MULTIRANGE": "YES",
    "GDAL_HTTP_MERGE_CONSECUTIVE_RANGES": "YES",
    "GDAL_HTTP_MAX_RETRY": "4",
    "GDAL_HTTP_RETRY_DELAY": "2",
    "VSI_CACHE": "TRUE",
}.items():
    os.environ.setdefault(_k, _v)

ROOT = Path(__file__).resolve().parent
DOWNLOADS = ROOT / "downloads" / "tasmania"
OUTPUTS = ROOT / "outputs" / "tasmania"

FIRE_LAYER = ("https://services.thelist.tas.gov.au/arcgis/rest/services/Public/"
              "EmergencyManagementPublic/MapServer/1")
DEA_STAC_SEARCH = "https://explorer.dea.ga.gov.au/stac/search"
DEA_COLLECTIONS = ["ga_s2am_ard_3", "ga_s2bm_ard_3", "ga_s2cm_ard_3"]
COASTLINE_URL = "https://naciscdn.org/naturalearth/10m/physical/ne_10m_land.zip"

CRS = "EPSG:3577"          # GDA94 / Australian Albers, equal-area, metres (as in the workshop data)
RES = 20.0                 # metres; the NIR-narrow and SWIR bands are 20 m native
COARSE = 2                 # recovery summers are read at RES * COARSE (40 m) to save memory
SMALL_FIRE_HA = 1000       # fires smaller than this (official area) keep RES for every summer:
                           # at 40 m a small fire has too few pixels per severity class (>= 50)
BUFFER_M = 3000.0          # study area = fire extent plus this margin
NODATA = -999              # DEA ARD nodata
SCALE = 10000.0            # DEA ARD reflectance scale
FMASK_CLEAR = 1            # fmask: 0 nodata, 1 clear, 2 cloud, 3 shadow, 4 snow, 5 water
BANDS = {                  # our name -> DEA asset (Sentinel-2 band)
    "red": "nbart_red",        # B04, 10 m
    "nir": "nbart_nir_2",      # B8A, 20 m
    "swir1": "nbart_swir_2",   # B11, 20 m
    "swir2": "nbart_swir_3",   # B12, 20 m
}
FMASK = "oa_fmask"
MIN_SEASON = 20152016      # Sentinel-2 era

# Key & Benson (2006) dNBR thresholds (generic, not calibrated for Tasmanian vegetation)
SEVERITY = [
    ("Unburnt", -math.inf, 0.10),
    ("Low", 0.10, 0.27),
    ("Moderate-low", 0.27, 0.44),
    ("Moderate-high", 0.44, 0.66),
    ("High", 0.66, math.inf),
]
SEVERITY_COLOURS = ["#d9d9d9", "#fed976", "#fd8d3c", "#e31a1c", "#800026"]  # grey + YlOrRd


def _load_sites(path: Path) -> dict:
    """Built-in sites from sites.csv (shared with tas_fires.R): key -> fire_id, groups, region."""
    with open(path, newline="") as f:
        return {r["key"]: dict(fire_id=int(r["fire_id"]), groups=r["groups"].split(";"),
                               region=r["region"]) for r in csv.DictReader(f)}


# Tabs on the landing page: group key -> (tab label, heading, lede, overview figure)
GROUPS = {
    "tasmania": ("Major fires", "Recent Tasmanian bushfires",
                 "large and significant fires across the state", "overview.png"),
    "fern-tree": ("Near Fern Tree", "Bushfires near Fern Tree",
                  "smaller fires within about 20 km of Fern Tree, on kunanyi/Mount Wellington",
                  "overview_fern_tree.png"),
    "uncertain": ("Uncertain results", "Fires with weak or noisy results",
                  "small fires near Fern Tree whose burn signal is hard to tell apart from "
                  "change in unburnt land nearby, or is only a speckle of low severity",
                  "overview_uncertain.png"),
}
FERN_TREE = (147.262, -42.923)   # longitude, latitude
SITES = _load_sites(ROOT / "sites.csv")

ATTRIBUTION = ("Sentinel-2: Digital Earth Australia ARD (Geoscience Australia, CC BY 4.0; contains "
               "Copernicus data). Fire boundaries: Fire History, theLIST (c) State of Tasmania, "
               "CC BY 3.0 AU.")
DISCLAIMER = ("Not an official resource: exploratory results made with an AI coding agent at a "
              "workshop, not checked by a bushfire or fire-ecology expert.")
DISCLAIMER_OFFICIAL = "For official fire information, see the Tasmania Fire Service."
TFS_URL = "https://www.fire.tas.gov.au/"


# ----------------------------------------------------------------------------------------------
# Small pure functions (tested in tests/test_tas_fires.py)
# ----------------------------------------------------------------------------------------------

def to_reflectance(dn: np.ndarray) -> np.ndarray:
    """DEA stored values (reflectance x 10,000, nodata -999) to reflectance 0-1, NaN for nodata."""
    a = np.asarray(dn, dtype="float32")
    out = a / SCALE
    out[(a == NODATA) | (a < 0)] = np.nan
    return out


def norm_diff(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """(a - b) / (a + b), NaN where undefined."""
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (a - b) / (a + b)
    out[~np.isfinite(out)] = np.nan
    return out


def nbr(c: dict) -> np.ndarray:
    """Normalised Burn Ratio = (NIR - SWIR2) / (NIR + SWIR2)."""
    return norm_diff(c["nir"], c["swir2"])


def ndvi(c: dict) -> np.ndarray:
    """Normalised Difference Vegetation Index = (NIR - red) / (NIR + red)."""
    return norm_diff(c["nir"], c["red"])


def classify_dnbr(dnbr: np.ndarray) -> np.ndarray:
    """dNBR to severity class index 0..4 (see SEVERITY); -1 where dNBR is NaN."""
    edges = [s[2] for s in SEVERITY[:-1]]
    cls = np.digitize(np.nan_to_num(dnbr, nan=0.0), edges).astype("int8")
    cls[~np.isfinite(dnbr)] = -1
    return cls


def pre_window(ignition: dt.date, days: int = 120) -> tuple[dt.date, dt.date]:
    return ignition - dt.timedelta(days=days), ignition - dt.timedelta(days=1)


def immediate_window(end: dt.date, days: int = 120) -> tuple[dt.date, dt.date]:
    return end + dt.timedelta(days=1), end + dt.timedelta(days=days)


def summer_window(year: int) -> tuple[dt.date, dt.date]:
    """Summer labelled by the year it ends in: 1 Dec (year-1) to 31 Mar (year)."""
    return dt.date(year - 1, 12, 1), dt.date(year, 3, 31)


def first_summer_after(end: dt.date) -> int:
    """First summer that starts after the fire was declared out."""
    year = end.year + 1
    while summer_window(year)[0] <= end:
        year += 1
    return year


def latest_summer(today: dt.date) -> int:
    """Most recent summer that has finished by `today`."""
    return today.year if today >= dt.date(today.year, 4, 1) else today.year - 1


def pre_summer_year(ignition: dt.date) -> int:
    """Last full summer (1 Dec - 31 Mar) that ended before the fire started."""
    year = ignition.year
    while summer_window(year)[1] >= ignition:
        year -= 1
    return year


# ----------------------------------------------------------------------------------------------
# GDAL (osgeo.gdal / osr / ogr): all raster reading, warping, writing and coordinate transforms
# ----------------------------------------------------------------------------------------------

def _gdal():
    from osgeo import gdal
    gdal.UseExceptions()
    return gdal


def _srs(code: str):
    """osr.SpatialReference in longitude/latitude (x/y) axis order, whatever the CRS says."""
    from osgeo import osr
    osr.UseExceptions()
    s = osr.SpatialReference()
    s.SetFromUserInput(code)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return s


def _ct(src: str, dst: str):
    from osgeo import osr
    return osr.CoordinateTransformation(_srs(src), _srs(dst))


def true_north_angle(x: float, y: float) -> float:
    """Angle of true north, in degrees clockwise from grid up, at EPSG:3577 point (x, y)."""
    lon, lat, _ = _ct(CRS, "EPSG:4326").TransformPoint(x, y)
    x2, y2, _ = _ct("EPSG:4326", CRS).TransformPoint(lon, lat + 0.05)
    return math.degrees(math.atan2(x2 - x, y2 - y))


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


@dataclass(frozen=True)
class Grid:
    """A north-up raster grid in EPSG:3577."""
    x0: float   # left edge
    y0: float   # top edge
    res: float
    width: int
    height: int

    @property
    def geotransform(self) -> tuple[float, ...]:
        """GDAL geotransform: (left, pixel width, 0, top, 0, -pixel height)."""
        return (self.x0, self.res, 0.0, self.y0, 0.0, -self.res)

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        return (self.x0, self.y0 - self.height * self.res, self.x0 + self.width * self.res, self.y0)

    @property
    def extent(self) -> list[float]:
        """For matplotlib imshow: left, right, bottom, top."""
        l, b, r, t = self.bounds
        return [l, r, b, t]

    @property
    def pixel_ha(self) -> float:
        return self.res * self.res / 10000.0

    def coarsen(self, f: int) -> "Grid":
        return Grid(self.x0, self.y0, self.res * f, self.width // f, self.height // f)

    def bounds_lonlat(self) -> list[float]:
        return list(_ct(CRS, "EPSG:4326").TransformBounds(*self.bounds, 21))


def grid_for_bounds(bounds: tuple[float, float, float, float], buffer_m: float = BUFFER_M,
                    res: float = RES) -> Grid:
    """Grid covering `bounds` (EPSG:3577) plus a buffer, snapped to 2*res so it can be coarsened."""
    step = res * COARSE
    minx, miny, maxx, maxy = bounds
    x0 = math.floor((minx - buffer_m) / step) * step
    x1 = math.ceil((maxx + buffer_m) / step) * step
    y0 = math.floor((miny - buffer_m) / step) * step
    y1 = math.ceil((maxy + buffer_m) / step) * step
    return Grid(x0, y1, res, int(round((x1 - x0) / res)), int(round((y1 - y0) / res)))


# ----------------------------------------------------------------------------------------------
# Fire History (theLIST)
# ----------------------------------------------------------------------------------------------

def _get_json(url: str, params: dict | None = None, body: dict | None = None,
              tries: int = 4) -> dict:
    import requests
    for i in range(tries):
        try:
            if body is not None:
                r = requests.post(url, json=body, timeout=180)
            else:
                r = requests.get(url, params=params, timeout=180)
            r.raise_for_status()
            d = r.json()
            if isinstance(d, dict) and "error" in d:
                raise RuntimeError(f"{url}: {d['error']}")
            return d
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 ** (i + 1))
    raise RuntimeError("unreachable")


def _ms_to_date(ms) -> dt.date | None:
    if ms is None:
        return None
    return dt.datetime.fromtimestamp(ms / 1000, dt.timezone.utc).date()


FIRE_FIELDS = "FIRE_ID,FIRE_NAME,FIRE_AR_HA,IGN_DATE,INCD_END,IGN_SEASON,IGN_CAUSE1,FIRE_TYPE"


def _fire_record(attrs: dict) -> dict:
    return dict(
        fire_id=attrs["FIRE_ID"], name=(attrs["FIRE_NAME"] or "").strip(), area_ha=attrs["FIRE_AR_HA"],
        ignition=_ms_to_date(attrs["IGN_DATE"]), end=_ms_to_date(attrs.get("INCD_END")),
        season=int(attrs["IGN_SEASON"]) if attrs.get("IGN_SEASON") else None,
        cause=attrs.get("IGN_CAUSE1"), fire_type=attrs.get("FIRE_TYPE"))


def search_fires(text: str = "", lon: float | None = None, lat: float | None = None,
                 limit: int = 25) -> list[dict]:
    """Bushfires since 2015-16 by name and/or containing a point, largest first."""
    where = f"FIRE_TYPE='Bushfire' AND IGN_SEASON>={MIN_SEASON}"
    if text:
        safe = text.upper().replace("'", "''")
        where += f" AND UPPER(FIRE_NAME) LIKE '%{safe}%'"
    params = dict(where=where, outFields=FIRE_FIELDS, returnGeometry="false",
                  orderByFields="FIRE_AR_HA DESC", resultRecordCount=limit, f="json")
    if lon is not None and lat is not None:
        params.update(geometry=f"{lon},{lat}", geometryType="esriGeometryPoint", inSR=4326,
                      spatialRel="esriSpatialRelIntersects", orderByFields="IGN_DATE DESC")
    d = _get_json(FIRE_LAYER + "/query", params)
    return [_fire_record(f["attributes"]) for f in d.get("features", [])]


def get_fire(fire_id: int) -> tuple[dict, object]:
    """Fire record and boundary (shapely, EPSG:4326). Cached in downloads/tasmania/fires/."""
    from shapely.geometry import shape
    from shapely.ops import unary_union
    from shapely.validation import make_valid

    path = DOWNLOADS / "fires" / f"{fire_id}.geojson"
    if not path.exists():
        params = dict(where=f"FIRE_ID={int(fire_id)}", outFields=FIRE_FIELDS, returnGeometry="true",
                      outSR=4326, maxAllowableOffset=0.0001, geometryPrecision=6, f="geojson")
        d = _get_json(FIRE_LAYER + "/query", params)
        if not d.get("features"):
            raise SystemExit(f"Fire {fire_id} not found in Fire History")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(d))
    d = json.loads(path.read_text())
    rec = _fire_record(d["features"][0]["properties"])
    geom = unary_union([make_valid(shape(f["geometry"])).buffer(0) for f in d["features"]])
    return rec, geom


def fires_in_area(bounds_ll: list[float], start: dt.date, end: dt.date) -> list:
    """Boundaries (EPSG:4326) of every fire (any type) in an area that started in [start, end]."""
    from shapely.geometry import shape
    params = dict(
        where=f"IGN_DATE >= DATE '{start.isoformat()}' AND IGN_DATE <= DATE '{end.isoformat()}'",
        geometry=",".join(f"{v:.5f}" for v in bounds_ll), geometryType="esriGeometryEnvelope",
        inSR=4326, spatialRel="esriSpatialRelIntersects", outFields="FIRE_ID,FIRE_NAME",
        returnGeometry="true", outSR=4326, maxAllowableOffset=0.0003, f="geojson")
    d = _get_json(FIRE_LAYER + "/query", params)
    out = []
    for f in d.get("features", []):
        if f.get("geometry"):
            out.append(shape(f["geometry"]))
    return out


# ----------------------------------------------------------------------------------------------
# Sentinel-2 from Digital Earth Australia
# ----------------------------------------------------------------------------------------------

def stac_search(bbox_ll: list[float], start: dt.date, end: dt.date) -> list[dict]:
    body = {"collections": DEA_COLLECTIONS, "bbox": bbox_ll, "limit": 250,
            "datetime": f"{start.isoformat()}T00:00:00Z/{end.isoformat()}T23:59:59Z"}
    import requests
    items, url, method = [], DEA_STAC_SEARCH, "POST"
    for _ in range(100):
        if method == "POST":
            d = _get_json(url, body=body)
        else:
            d = requests.get(url, timeout=180).json()
        items += d.get("features", [])
        nxt = [l for l in d.get("links", []) if l.get("rel") == "next"]
        if not nxt:
            break
        link = nxt[0]
        url, method = link["href"], link.get("method", "GET").upper()
        if method == "POST":
            body = {**body, **link.get("body", {})} if link.get("merge") else link.get("body", body)
    return items


def _href(item: dict, asset: str) -> str:
    h = item["assets"][asset]["href"]
    # Straight from the public S3 bucket (Sydney): the data.dea.ga.gov.au CloudFront front end
    # rate-limits bursts of range requests with HTTP 403.
    return h.replace("s3://dea-public-data/", "https://dea-public-data.s3.ap-southeast-2.amazonaws.com/")


def group_by_date(items: list[dict]) -> list[tuple[str, list[dict], float]]:
    """[(date, items, mean tile cloud %)], least cloudy first."""
    groups: dict[str, list[dict]] = {}
    for it in items:
        groups.setdefault(it["properties"]["datetime"][:10], []).append(it)
    out = []
    for d, its in groups.items():
        cc = [it["properties"].get("eo:cloud_cover") for it in its]
        cc = [c for c in cc if c is not None]
        out.append((d, its, float(np.mean(cc)) if cc else 100.0))
    return sorted(out, key=lambda t: (t[2], t[0]))


def _warp(url: str, grid: Grid) -> np.ndarray:
    """Band 1 of `url` (a path or an https URL) warped onto `grid`, nearest neighbour.

    Only the bytes covering the grid are read (GDAL /vsicurl/ range requests). Areas the file
    doesn't cover come back as the file's nodata value. Overviews are not used, so pixels are
    always sampled from full resolution."""
    gdal = _gdal()
    src = "/vsicurl/" + url if url.startswith(("http://", "https://")) else url
    for i in range(4):
        try:
            ds = gdal.Warp("", src, format="MEM", dstSRS=CRS, outputBounds=grid.bounds,
                           width=grid.width, height=grid.height, resampleAlg="near",
                           overviewLevel="NONE", multithread=False)
            return ds.GetRasterBand(1).ReadAsArray()
        except Exception as exc:
            missing = "403" in str(exc) or "404" in str(exc)   # S3 answers 403 for a missing file
            if i == 3 or (missing and i == 2):
                raise
            time.sleep(3 * (i + 1))
    raise RuntimeError("unreachable")


def read_date(items: list[dict], grid: Grid) -> tuple[dict, float, np.ndarray, np.ndarray]:
    """Mosaic one date's tiles. Returns clear (fmask) reflectance per band, the clear fraction of
    what fmask saw as not water, and two masks: seen-not-water (clear, cloud, shadow, snow) and
    water. Across dates, a pixel ever labelled water is treated as sea, so a coastal study area is
    not penalised for its sea."""
    out = {b: np.full((grid.height, grid.width), np.nan, "float32") for b in BANDS}
    filled = np.zeros((grid.height, grid.width), bool)
    land = np.zeros((grid.height, grid.width), bool)
    water = np.zeros((grid.height, grid.width), bool)
    for it in items:
        try:
            fm = _warp(_href(it, FMASK), grid)
            land |= (fm >= 1) & (fm <= 4)
            water |= fm == 5
            ok = (fm == FMASK_CLEAR) & ~filled
            if not ok.any():
                continue
            refl = {b: to_reflectance(_warp(_href(it, a), grid)) for b, a in BANDS.items()}
        except Exception as exc:   # e.g. a file missing from the archive (S3 answers 403)
            print(f"    skipped tile {it['id']} ({it['properties']['datetime'][:10]}): "
                  f"{str(exc)[:80]}")
            continue
        for r in refl.values():
            ok &= np.isfinite(r)
        for b in BANDS:
            out[b][ok] = refl[b][ok]
        filled |= ok
    return out, float(filled.sum() / max(land.sum(), 1)), land, water


def _write_composite(path: Path, comp: dict, grid: Grid, meta: dict) -> None:
    gdal = _gdal()
    names = list(BANDS) + ["n_clear"]
    path.parent.mkdir(parents=True, exist_ok=True)
    ds = gdal.GetDriverByName("GTiff").Create(
        str(path), grid.width, grid.height, len(names), gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"])
    ds.SetGeoTransform(grid.geotransform)
    ds.SetSpatialRef(_srs(CRS))
    for i, n in enumerate(names, 1):
        band = ds.GetRasterBand(i)
        band.SetNoDataValue(float("nan"))
        band.SetDescription(n)
        band.WriteArray(comp[n].astype("float32"))
    ds = None   # closes and flushes the file
    path.with_suffix(".json").write_text(json.dumps(meta, indent=1))


def _read_composite(path: Path) -> tuple[dict, dict]:
    gdal = _gdal()
    ds = gdal.Open(str(path))
    comp = {}
    for i in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(i)
        comp[band.GetDescription()] = band.ReadAsArray()
    ds = None
    return comp, json.loads(path.with_suffix(".json").read_text())


def composite(site_dir: Path, label: str, bbox_ll: list[float], start: dt.date, end: dt.date,
              grid: Grid, max_dates: int = 4, min_clear: float = 0.25, refresh: bool = False,
              min_cover: float = 0.9, log=print) -> tuple[dict | None, dict]:
    """Median of up to `max_dates` clear-ish dates in [start, end]. Cached as a GeoTIFF."""
    path = site_dir / f"{label}_{int(grid.res)}m.tif"
    if path.exists() and not refresh:
        comp, meta = _read_composite(path)
        # accept if well covered, or if the last build already tried every candidate date
        # (composites built before coverage was measured over land lack "land_fraction": a
        # patchy one of those is rebuilt once, as the sea may be what made it look patchy)
        if meta.get("valid_fraction", 0) >= min_cover - 0.05 or (
                meta.get("searched_all") and meta.get("coverage_basis") == "land-v2"):
            log(f"  {label}: cached, dates {', '.join(meta['dates'])}")
            return comp, meta
        log(f"  {label}: cached composite covers only {meta.get('valid_fraction')}: rebuilding")
    groups = group_by_date(stac_search(bbox_ll, start, end))
    candidates = [g for g in groups if g[2] < 80]
    log(f"  {label}: {start} to {end}: {len(groups)} dates found, {len(candidates)} under 80% tile "
        f"cloud; reading up to {max_dates}")
    stacks = {b: [] for b in BANDS}
    used, fracs = [], []
    pos = 0
    with ThreadPoolExecutor(max_workers=3) as ex:
        # Least-cloudy dates first, but keep going until the study area is covered: dates from
        # one satellite pass can all cover the same part of a site that spans two passes.
        cover = np.zeros((grid.height, grid.width), bool)
        seen = np.zeros((grid.height, grid.width), bool)
        water_any = np.zeros((grid.height, grid.width), bool)
        land_eff = lambda: (seen & ~water_any) | cover
        covered = lambda: float(cover.sum() / max(land_eff().sum(), 1))
        limit = min(len(candidates), 3 * max_dates)
        while pos < limit and (len(used) < max_dates or covered() < min_cover):
            batch = candidates[pos:pos + 3]
            pos += len(batch)
            for (d, its, cc), (arrs, frac, land, water) in zip(
                    batch, ex.map(lambda g: read_date(g[1], grid), batch)):
                seen |= land
                water_any |= water
                log(f"    {d}: tile cloud {cc:4.0f}%, clear over the land {frac:5.1%}")
                if frac >= min_clear:
                    used.append(d)
                    fracs.append(round(frac, 3))
                    cover |= np.isfinite(arrs["nir"])
                    for b in BANDS:
                        stacks[b].append(arrs[b])
        log(f"    {len(used)} dates used; land covered {covered():.1%} "
            f"(land is {land_eff().mean():.0%} of the study area)")
    meta = dict(label=label, window=[start.isoformat(), end.isoformat()], dates=sorted(used),
                clear_fraction=dict(zip(used, fracs)), res_m=grid.res,
                collections=DEA_COLLECTIONS, n_dates_found=len(groups), searched_all=True)
    if not used:
        log(f"  {label}: no usable dates")
        return None, meta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        comp = {b: np.nanmedian(np.stack(stacks[b]), axis=0).astype("float32") for b in BANDS}
    comp["n_clear"] = np.isfinite(np.stack(stacks["nir"])).sum(axis=0).astype("float32")
    land_mask = land_eff()
    meta["valid_fraction"] = round(float(np.isfinite(comp["nir"])[land_mask].mean()), 3)
    meta["land_fraction"] = round(float(land_mask.mean()), 3)
    meta["coverage_basis"] = "land-v2"
    _write_composite(path, comp, grid, meta)
    return comp, meta


# ----------------------------------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------------------------------

def to_albers(geom):
    """A shapely geometry from longitude/latitude (EPSG:4326) to EPSG:3577 metres."""
    import shapely
    ct = _ct("EPSG:4326", CRS)

    def f(xy: np.ndarray) -> np.ndarray:
        if len(xy) == 0:
            return xy
        return np.asarray(ct.TransformPoints(np.ascontiguousarray(xy, "float64")))[:, :2]
    return shapely.transform(geom, f)


def rasterize(geoms: list, grid: Grid) -> np.ndarray:
    """True inside any of the (EPSG:3577) geometries: pixels whose centre is inside."""
    from osgeo import ogr
    gdal = _gdal()
    geoms = [g for g in geoms if g is not None and not g.is_empty]
    if not geoms:
        return np.zeros((grid.height, grid.width), bool)
    drv = ogr.GetDriverByName("MEM") or ogr.GetDriverByName("Memory")   # "MEM" from GDAL 3.11
    vds = drv.CreateDataSource("")
    lyr = vds.CreateLayer("geoms", srs=_srs(CRS))
    for g in geoms:
        feat = ogr.Feature(lyr.GetLayerDefn())
        feat.SetGeometry(ogr.CreateGeometryFromWkb(g.wkb))
        lyr.CreateFeature(feat)
    rds = gdal.GetDriverByName("MEM").Create("", grid.width, grid.height, 1, gdal.GDT_Byte)
    rds.SetGeoTransform(grid.geotransform)
    rds.SetSpatialRef(_srs(CRS))
    gdal.RasterizeLayer(rds, [1], lyr, burn_values=[1])
    return rds.GetRasterBand(1).ReadAsArray().astype(bool)


def _rings(geom):
    polys = getattr(geom, "geoms", [geom])
    for p in polys:
        if p.geom_type == "Polygon":
            yield p.exterior
            yield from p.interiors
        elif hasattr(p, "geoms"):
            yield from _rings(p)


# ----------------------------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------------------------

def add_map_furniture(ax, grid_or_extent, location="lower left"):
    """Scale bar (metres, EPSG:3577) and a true-north arrow."""
    from matplotlib_scalebar.scalebar import ScaleBar
    ext = grid_or_extent.extent if isinstance(grid_or_extent, Grid) else grid_or_extent
    ax.add_artist(ScaleBar(1, "m", location=location, box_alpha=0.75, font_properties={"size": 8}))
    cx, cy = (ext[0] + ext[1]) / 2, (ext[2] + ext[3]) / 2
    ang = math.radians(true_north_angle(cx, cy))
    bx, by, L = 0.93, 0.80, 0.11
    tip = (bx + L * math.sin(ang), by + L * math.cos(ang))
    ax.annotate("", xy=tip, xytext=(bx, by), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5), zorder=10)
    ax.text(tip[0], tip[1] + 0.015, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=9, fontweight="bold", zorder=10,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))


def plot_outline(ax, geom3577, **kw):
    style = dict(color="#00e5ff", lw=0.8)
    style.update(kw)
    for ring in _rings(geom3577):
        x, y = ring.xy
        ax.plot(x, y, **style)


def false_colour(c: dict, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """RGB = SWIR2, NIR, red, with a fixed per-band stretch; NaN -> white."""
    rgb = np.dstack([c["swir2"], c["nir"], c["red"]])
    out = np.clip((rgb - lo) / (hi - lo), 0, 1)
    out[~np.isfinite(rgb).all(axis=2)] = 1.0
    return out


def _fmt_dates(meta: dict) -> str:
    d = meta.get("dates", [])
    if not d:
        return "no clear images"
    a, b = dt.date.fromisoformat(d[0]), dt.date.fromisoformat(d[-1])
    span = a.strftime("%d %b %Y") if a == b else f"{a.strftime('%d %b %Y')} to {b.strftime('%d %b %Y')}"
    return f"{span} (median of {len(d)} date{'s' if len(d) > 1 else ''})"


# ----------------------------------------------------------------------------------------------
# Analysis of one fire
# ----------------------------------------------------------------------------------------------

def analyse_fire(key: str, fire_id: int, region: str = "", max_dates: int = 4, refresh: bool = False,
                 today: dt.date | None = None) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
    from matplotlib.patches import Patch

    today = today or dt.date.today()
    rec, geom_ll = get_fire(fire_id)
    if rec["ignition"] is None:
        raise SystemExit(f"Fire {fire_id} has no ignition date")
    end = rec["end"] or (rec["ignition"] + dt.timedelta(days=30))
    print(f"\n== {key}: {rec['name']} (FIRE_ID {fire_id}), ignited {rec['ignition']}, "
          f"out {end}, official area {rec['area_ha']:,.0f} ha ==")
    site_dir = DOWNLOADS / key
    out_dir = OUTPUTS / key
    out_dir.mkdir(parents=True, exist_ok=True)

    geom = to_albers(geom_ll)
    grid = grid_for_bounds(geom.bounds)
    if max(grid.width, grid.height) > 4500:   # very large fire: use coarser pixels
        grid = grid_for_bounds(geom.bounds, res=RES * 2)
    small = rec["area_ha"] < SMALL_FIRE_HA
    cgrid = grid if small else grid.coarsen(COARSE)
    bbox_ll = grid.bounds_lonlat()
    print(f"  grid {grid.width} x {grid.height} px at {grid.res:.0f} m "
          f"({grid.width * grid.res / 1000:.1f} x {grid.height * grid.res / 1000:.1f} km), EPSG:3577; "
          f"later summers at {cgrid.res:.0f} m" + (" (small fire)" if small else ""))

    # --- composites -----------------------------------------------------------------------------
    s, e = pre_window(rec["ignition"])
    pre, pre_m = composite(site_dir, "pre", bbox_ll, s, e, grid, max_dates, refresh=refresh)
    if pre is None or pre_m.get("valid_fraction", 0) < 0.6:
        s, e = pre_window(rec["ignition"], days=365)
        pre, pre_m = composite(site_dir, "pre365", bbox_ll, s, e, grid, max_dates, refresh=refresh)
    # Severity baseline: the last full summer before the fire, so before and after are compared in
    # the same season (spring-to-summer drying of grass and crops otherwise looks like a burn).
    py = pre_summer_year(rec["ignition"])
    base, base_m = None, {}
    if summer_window(py)[0] >= dt.date(2015, 7, 1):        # Sentinel-2 data in Tasmania
        s, e = summer_window(py)
        base, base_m = composite(site_dir, f"presummer{py}", bbox_ll, s, e, grid, max_dates,
                                 refresh=refresh)
    if base is None or base_m.get("valid_fraction", 0) < 0.8:
        why = ("no earlier Sentinel-2 summer" if base is None and not base_m else
               f"summer {py - 1}-{str(py)[2:]} too cloudy")
        base, base_m = pre, pre_m
        base_name = f"just before the fire ({pre_m['dates'][0][:4]}; {why})"
    else:
        base_name = f"summer {py - 1}-{str(py)[2:]}"
    if pre is None:
        raise SystemExit("No clear images in the year before the fire")
    s, e = immediate_window(end)
    imm, imm_m = composite(site_dir, "immediate", bbox_ll, s, e, grid, max_dates, refresh=refresh)
    first = first_summer_after(end)
    last = latest_summer(today)
    initial = first > last     # no full summer since the fire was out: initial assessment
    if initial:
        if imm is None:
            raise SystemExit("No clear imagery yet from after the fire")
        s1, s1_m = imm, imm_m
        nxt = end + dt.timedelta(days=1)
        first = nxt.year + (1 if nxt.month == 12 else 0)
        sev_name = "just after the fire (initial assessment)"
        print(f"  no full summer since the fire was out: severity from the images just after it")
    else:
        s, e = summer_window(first)
        s1, s1_m = composite(site_dir, f"summer{first}", bbox_ll, s, e, grid, max_dates,
                             refresh=refresh)
        sev_name = f"summer {first - 1}-{str(first)[2:]}"
    if pre is None or s1 is None:
        raise SystemExit("Not enough clear imagery before the fire or in the first summer after it")
    summers = {first: (s1, s1_m, grid)}
    for y in range(first + 1, last + 1):
        s, e = summer_window(y)
        c, m = composite(site_dir, f"summer{y}", bbox_ll, s, e, cgrid, max(2, max_dates - 1),
                         refresh=refresh)
        if c is not None:
            summers[y] = (c, m, cgrid)

    # --- burn severity --------------------------------------------------------------------------
    inside = rasterize([geom], grid)
    others_ll = fires_in_area(bbox_ll, pre_window(rec["ignition"], 365)[0], summer_window(last)[1])
    others = [to_albers(g).buffer(200) for g in others_ll]
    burnt_any = rasterize(others + [geom.buffer(500)], grid)
    control = ~burnt_any

    nbr_pre, nbr_s1 = nbr(base), nbr(s1)
    dnbr_raw = nbr_pre - nbr_s1
    ok_ctrl = control & np.isfinite(dnbr_raw)
    offset = float(np.nanmedian(dnbr_raw[ok_ctrl])) if ok_ctrl.sum() > 1000 else 0.0
    dnbr = dnbr_raw - offset
    cls = classify_dnbr(dnbr)
    px = grid.pixel_ha
    valid_in = inside & np.isfinite(dnbr)
    burnt = dnbr >= SEVERITY[1][1]
    stats = dict(
        site=key, fire_id=fire_id, fire_name=rec["name"], region=region,
        ignition=rec["ignition"].isoformat(), declared_out=end.isoformat(), cause=rec["cause"],
        official_area_ha=round(rec["area_ha"], 1),
        boundary_area_ha=round(float(inside.sum() * px), 1),
        valid_inside_pct=round(100 * valid_in.sum() / max(inside.sum(), 1), 1),
        dnbr_offset=round(offset, 4),
        burnt_inside_ha=round(float((burnt & valid_in).sum() * px), 1),
        burnt_inside_pct=round(100 * (burnt & valid_in).sum() / max(valid_in.sum(), 1), 1),
        burnt_outside_ha=round(float((burnt & control & np.isfinite(dnbr)).sum() * px), 1),
        control_ha=round(float(ok_ctrl.sum() * px), 1),
        pre_dates=";".join(pre_m["dates"]), baseline=base_name,
        baseline_dates=";".join(base_m["dates"]), immediate_dates=";".join(imm_m.get("dates", [])),
        severity_basis="initial (just after the fire)" if initial else f"summer {first}",
        severity_summer=first, severity_dates=";".join(s1_m["dates"]),
        res_m=grid.res, recovery_res_m=cgrid.res,
    )
    for i, (name, _, _) in enumerate(SEVERITY):
        n = int(((cls == i) & valid_in).sum())
        stats[f"ha_{slugify(name)}"] = round(n * px, 1)
        stats[f"pct_{slugify(name)}"] = round(100 * n / max(valid_in.sum(), 1), 1)

    print(f"  reflectance ranges (pre): " + ", ".join(
        f"{b} {np.nanpercentile(pre[b], 1):.3f}-{np.nanpercentile(pre[b], 99):.3f}" for b in BANDS))
    print(f"  severity baseline: {base_name}")
    print(f"  NBR pre {np.nanpercentile(nbr_pre, 1):.2f} to {np.nanpercentile(nbr_pre, 99):.2f}; "
          f"dNBR (offset {offset:+.3f} removed) {np.nanpercentile(dnbr, 1):.2f} to "
          f"{np.nanpercentile(dnbr, 99):.2f}")
    print(f"  boundary {stats['boundary_area_ha']:,.0f} ha (official {stats['official_area_ha']:,.0f}); "
          f"valid {stats['valid_inside_pct']}%; burnt (dNBR >= 0.10) {stats['burnt_inside_ha']:,.0f} ha "
          f"= {stats['burnt_inside_pct']}% of the boundary; dNBR >= 0.10 outside any fire "
          f"{stats['burnt_outside_ha']:,.0f} ha of {stats['control_ha']:,.0f} ha unburnt land")

    # --- recovery -------------------------------------------------------------------------------
    rows = []
    for y, (c, m, g) in sorted(summers.items()):
        f = 1 if g.res == grid.res else COARSE
        cls_g = cls[::f, ::f][:g.height, :g.width]
        in_g = inside[::f, ::f][:g.height, :g.width]
        ctrl_g = control[::f, ::f][:g.height, :g.width]
        pre_g = nbr_pre[::f, ::f][:g.height, :g.width]
        d = pre_g - nbr(c)                       # NBR drop from before the fire
        ctrl_drop = float(np.nanmedian(d[ctrl_g])) if np.isfinite(d[ctrl_g]).sum() > 100 else 0.0
        row = dict(site=key, summer=y, years_since=y - rec["ignition"].year,
                   dates=";".join(m["dates"]), control_drop=round(ctrl_drop, 4),
                   control_ndvi=round(float(np.nanmedian(ndvi(c)[ctrl_g])), 4))
        for i, (name, _, _) in enumerate(SEVERITY):
            sel = in_g & (cls_g == i) & np.isfinite(d)
            if sel.sum() < 50:
                row[f"deficit_{slugify(name)}"] = ""
                row[f"ndvi_{slugify(name)}"] = ""
                continue
            row[f"deficit_{slugify(name)}"] = round(float(np.nanmedian(d[sel])) - ctrl_drop, 4)
            row[f"ndvi_{slugify(name)}"] = round(float(np.nanmedian(ndvi(c)[sel])), 4)
        rows.append(row)
    # pre-fire NDVI per class (reference line)
    ndvi_pre = ndvi(base)
    pre_ndvi = {slugify(n): float(np.nanmedian(ndvi_pre[inside & (cls == i)]))
                for i, (n, _, _) in enumerate(SEVERITY) if (inside & (cls == i)).sum() >= 50}
    pre_ndvi_ctrl = float(np.nanmedian(ndvi_pre[control]))
    hi = slugify("High")
    mh = slugify("Moderate-high")
    for k in (hi, mh):
        first_def = rows[0].get(f"deficit_{k}")
        last_def = rows[-1].get(f"deficit_{k}")
        if first_def not in ("", None) and last_def not in ("", None) and first_def > 0:
            stats[f"recovered_{k}_pct"] = round(100 * (1 - last_def / first_def), 1)
        else:
            stats[f"recovered_{k}_pct"] = ""
    stats["latest_summer"] = max(summers)

    # --- figures --------------------------------------------------------------------------------
    caption = (DISCLAIMER + f" {rec['name']} fire (FIRE_ID {fire_id}), ignited {rec['ignition']:%d %b %Y}, "
               f"declared out {end:%d %b %Y}. " + ATTRIBUTION)
    title_site = f"{rec['name']} fire, {rec['ignition'].year}" + (f": {region}" if region else "")

    # 1. false colour before / after
    panels = [("Before the fire", pre, pre_m, grid)]
    if imm is not None:
        panels.append(("Just after the fire", imm, imm_m, grid))
    if not initial:
        panels.append((f"First summer after ({first - 1}-{str(first)[2:]})", s1, s1_m, grid))
    if last != first and last in summers:
        panels.append((f"Latest summer ({last - 1}-{str(last)[2:]})", *summers[last]))
    lo = np.array([np.nanpercentile(pre[b], 2) for b in ("swir2", "nir", "red")])
    hi_ = np.array([np.nanpercentile(pre[b], 98) for b in ("swir2", "nir", "red")])
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n, 4.6 * grid.height / grid.width + 1.4),
                             constrained_layout=True)
    for ax, (t, c, m, g) in zip(np.atleast_1d(axes), panels):
        ax.imshow(false_colour(c, lo, hi_), extent=g.extent, interpolation="nearest")
        plot_outline(ax, geom)
        ax.set_title(f"{t}\n{_fmt_dates(m)}", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        add_map_furniture(ax, g)
    fig.suptitle(f"{title_site}\nSentinel-2 false colour (R = SWIR2 B12, G = NIR B8A, B = red B04), "
                 f"same stretch for all panels (reflectance {lo.min():.2f}-{hi_.max():.2f}); "
                 "burnt ground shows magenta-pink, healthy forest green, buttongrass and sedgeland purple; cyan = official boundary",
                 fontsize=9)
    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=7, wrap=True)
    fig.savefig(out_dir / "before_after.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 2. dNBR and severity
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0 * grid.height / grid.width + 1.6),
                             constrained_layout=True)
    im = axes[0].imshow(dnbr, extent=grid.extent, cmap="BrBG_r",
                        norm=TwoSlopeNorm(vcenter=0, vmin=-0.8, vmax=0.8), interpolation="nearest")
    plot_outline(axes[0], geom, color="#0072B2", lw=0.6)
    cb = fig.colorbar(im, ax=axes[0], shrink=0.8, extend="both")
    cb.set_label("dNBR (unitless; > 0 = burn, < 0 = greener than before)")
    axes[0].set_title(f"dNBR: {base_name} minus {sev_name}"
                      f"\noffset {offset:+.3f} from unburnt land removed", fontsize=9)
    cmap = ListedColormap(SEVERITY_COLOURS)
    cmap.set_bad("white")
    axes[1].imshow(np.ma.masked_less(cls, 0), extent=grid.extent, cmap=cmap,
                   norm=BoundaryNorm(np.arange(-0.5, len(SEVERITY)), len(SEVERITY)),
                   interpolation="nearest")
    plot_outline(axes[1], geom, color="#0072B2", lw=0.6)
    axes[1].set_title(f"Burn severity class (Key & Benson dNBR thresholds), {sev_name}"
                      f"\n{stats['burnt_inside_pct']}% of the official "
                      f"boundary burnt (dNBR >= 0.10); white = no clear data", fontsize=9)
    handles = [Patch(color=c, label=f"{n} ({stats['pct_' + slugify(n)]}%)")
               for (n, _, _), c in zip(SEVERITY, SEVERITY_COLOURS)]
    handles.append(Patch(facecolor="none", edgecolor="#0072B2", label="Official boundary"))
    axes[1].legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3,
                   fontsize=7, title="% of the official boundary area", title_fontsize=7,
                   frameon=False)
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])
        add_map_furniture(ax, grid)
    fig.suptitle(title_site, fontsize=10)
    fig.text(0.5, -0.02, caption + " Thresholds are generic, not calibrated for Tasmanian "
             "vegetation.", ha="center", va="top", fontsize=7, wrap=True)
    fig.savefig(out_dir / "dnbr_severity.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 3. recovery
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    years = [r["summer"] for r in rows]
    for i, ((name, _, _), col) in enumerate(zip(SEVERITY, SEVERITY_COLOURS)):
        k = slugify(name)
        vals = [r[f"deficit_{k}"] for r in rows]
        if all(v == "" for v in vals):
            continue
        y = [np.nan if v == "" else v for v in vals]
        nd = [np.nan if r[f"ndvi_{k}"] == "" else r[f"ndvi_{k}"] for r in rows]
        c = "#737373" if i == 0 else col
        axes[0].plot(years, y, "o-", color=c, label=name)
        axes[1].plot(years, nd, "o-", color=c, label=name)
        if k in pre_ndvi:
            axes[1].axhline(pre_ndvi[k], color=c, ls=":", lw=1)
    axes[1].plot(years, [r["control_ndvi"] for r in rows], "s--", color="#0072B2",
                 label="Unburnt land nearby")
    axes[1].axhline(pre_ndvi_ctrl, color="#0072B2", ls=":", lw=1)
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].set_ylabel("NBR deficit vs before the fire\n(unitless; 0 = fully recovered)")
    axes[0].set_title("Burn-ratio deficit by severity class, relative to unburnt land", fontsize=9)
    axes[1].set_ylabel("Median NDVI (unitless)")
    axes[1].set_title("Greenness (NDVI) by severity class; dotted = before the fire", fontsize=9)
    for ax in axes:
        ax.set_xlabel("Summer (Dec-Mar, labelled by the year it ends)")
        ax.set_xticks(years)
        ax.set_xlim(years[0] - 0.5, years[-1] + 0.5)   # both panels, even with one point plotted
        ax.tick_params(axis="x", labelrotation=45)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.suptitle(f"{title_site}: recovery, summers {years[0] - 1}-{str(years[0])[2:]} to "
                 f"{years[-1] - 1}-{str(years[-1])[2:]}", fontsize=10)
    fig.text(0.5, -0.02, caption, ha="center", va="top", fontsize=7, wrap=True)
    fig.savefig(out_dir / "recovery.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 4. severity overlay for the interactive map (EPSG:4326 PNG)
    _severity_overlay(cls, grid, out_dir / "severity_overlay.png", out_dir / "severity_overlay.json")
    (out_dir / "boundary.geojson").write_text(json.dumps(
        dict(type="FeatureCollection", features=[dict(type="Feature", properties=dict(
            site=key, name=rec["name"]), geometry=json.loads(json.dumps(
                geom_ll.simplify(0.0005).__geo_interface__)))])))

    _write_csv(out_dir / "summary.csv", [stats])
    _write_csv(out_dir / "recovery.csv", rows)
    print(f"  wrote {out_dir.relative_to(ROOT)}/: before_after.png, dnbr_severity.png, "
          f"recovery.png, summary.csv, recovery.csv")
    return stats


def _severity_overlay(cls: np.ndarray, grid: Grid, png: Path, meta: Path) -> None:
    from matplotlib.colors import to_rgba
    from PIL import Image
    gdal = _gdal()
    src = gdal.GetDriverByName("MEM").Create("", grid.width, grid.height, 1, gdal.GDT_Int16)
    src.SetGeoTransform(grid.geotransform)
    src.SetSpatialRef(_srs(CRS))
    src.GetRasterBand(1).SetNoDataValue(-1)
    src.GetRasterBand(1).WriteArray(cls.astype("int16"))
    # GDAL's default longitude/latitude grid for this raster (a VRT: nothing is computed yet),
    # then the same extent origin with pixels `scale` times larger, so the PNG is at most ~1500 px
    vrt = gdal.Warp("", src, format="VRT", dstSRS="EPSG:4326")
    w, h = vrt.RasterXSize, vrt.RasterYSize
    west, dx, _, north, _, dy = vrt.GetGeoTransform()
    scale = max(1, int(math.ceil(max(w, h) / 1500)))
    w, h = w // scale, h // scale
    east, south = west + dx * scale * w, north + dy * scale * h
    out = gdal.Warp("", src, format="MEM", dstSRS="EPSG:4326", outputBounds=(west, south, east, north),
                    width=w, height=h, resampleAlg="near", srcNodata=-1, dstNodata=-1)
    dst = out.GetRasterBand(1).ReadAsArray()
    rgba = np.zeros((h, w, 4), "uint8")
    for i, c in enumerate(SEVERITY_COLOURS):
        if i == 0:
            continue   # unburnt transparent
        rgba[dst == i] = [int(255 * v) for v in to_rgba(c, 0.8)]
    Image.fromarray(rgba).save(png)
    meta.write_text(json.dumps(dict(bounds=[[south, west], [north, east]])))


def _write_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _read_csv(path: Path) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------------------------------
# Overview across sites
# ----------------------------------------------------------------------------------------------

def _coastline():
    """Natural Earth 1:10m land polygons clipped to Tasmania (EPSG:4326, shapely)."""
    import shapefile
    import requests
    from shapely.geometry import box, shape
    path = DOWNLOADS / "ne_10m_land.zip"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(requests.get(COASTLINE_URL, timeout=300).content)
    with zipfile.ZipFile(path) as z:
        r = shapefile.Reader(shp=io.BytesIO(z.read("ne_10m_land.shp")),
                             shx=io.BytesIO(z.read("ne_10m_land.shx")),
                             dbf=io.BytesIO(z.read("ne_10m_land.dbf")))
        clip = box(143.5, -44.0, 149.0, -39.2)
        geoms = [shape(s.__geo_interface__) for s in r.shapes()]
    return [g.intersection(clip) for g in geoms if g.intersects(clip)]


def overview() -> None:
    """summary_all.csv, one overview figure per site group, the interactive map and the landing page."""
    summaries = []
    for d in sorted(OUTPUTS.glob("*/summary.csv")):
        summaries += _read_csv(d)
    if not summaries:
        raise SystemExit("No finished sites yet: run --site first")
    summaries.sort(key=lambda r: r["ignition"])
    for r in summaries:   # sites run with --fire-id or --lon/--lat go with the major fires
        r["groups"] = ";".join(SITES.get(r["site"], {}).get("groups", ["tasmania"]))
    _write_csv(OUTPUTS / "summary_all.csv", summaries)
    colours = _year_colours(summaries)[0]   # colour = fire year, as in the overview figures
    groups = {}
    for g in GROUPS:
        rows = [r for r in summaries if g in r["groups"].split(";")]
        if rows:
            groups[g] = rows
            _overview_figure(rows, g)
    _folium_map(summaries, colours)
    _landing_page(groups)
    print(f"wrote outputs/tasmania/: {', '.join(GROUPS[g][3] for g in groups)}, sites_map.html, "
          f"index.html and summary_all.csv ({len(summaries)} sites)")


def _year_colours(summaries: list[dict]):
    """Colour = fire year (viridis, trimmed before the pale yellow end): meaningful, colour-blind
    safe, and it scales to any number of fires. Returns site -> hex, the colormap, the norm and
    the years."""
    from matplotlib import colormaps
    from matplotlib.colors import Normalize, to_hex
    years = [int(r["ignition"][:4]) for r in summaries]
    year_norm = Normalize(min(years), max(max(years), min(years) + 1))
    year_cmap = colormaps["viridis"].resampled(256)
    colours = {r["site"]: to_hex(year_cmap(0.88 * year_norm(int(r["ignition"][:4]))))
               for r in summaries}
    return colours, year_cmap, year_norm, years


def _overview_figure(summaries: list[dict], group: str) -> None:
    """Map of the sites, the severity mix of each fire, and recovery curves, for one group."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from shapely.geometry import shape

    local = group != "tasmania"     # zoomed map with numbered dots, not the whole state
    colours, year_cmap, year_norm, years = _year_colours(summaries)
    number = {r["site"]: i + 1 for i, r in enumerate(summaries)}
    tag = (lambda r: f"{number[r['site']]}. ") if local else (lambda r: "")
    fig = plt.figure(figsize=(15, max(6.2, 2.5 + 0.28 * len(summaries))), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.15])
    ax = fig.add_subplot(gs[0])
    small_key = ""
    if local:
        geoms = {r["site"]: to_albers(shape(json.loads((OUTPUTS / r["site"] / "boundary.geojson")
                 .read_text())["features"][0]["geometry"])) for r in summaries}
        fx, fy, _ = _ct("EPSG:4326", CRS).TransformPoint(*FERN_TREE)
        xs = [g.representative_point().x for g in geoms.values()] + [fx]
        ys = [g.representative_point().y for g in geoms.values()] + [fy]
        pad = 3000
        x0, x1, y0, y1 = min(xs) - pad, max(xs) + pad, min(ys) - pad, max(ys) + 3 * pad   # room for N
        for g in _coastline():
            for ring in _rings(to_albers(g)):
                x, y = ring.xy
                ax.fill(x, y, color="#f0f0f0", ec="#636363", lw=0.5)
        for r in summaries:
            g = geoms[r["site"]]
            for ring in _rings(g):
                x, y = ring.xy
                ax.fill(x, y, color=colours[r["site"]], lw=0)
            c = g.representative_point()
            ax.plot(c.x, c.y, "o", ms=7, mfc=colours[r["site"]], mec="black", mew=0.6, zorder=3)
            ax.annotate(str(number[r["site"]]), (c.x, c.y), xytext=(5, 4), textcoords="offset points",
                        fontsize=8, fontweight="bold", zorder=4,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
        ax.plot(fx, fy, "*", ms=14, mfc="black", mec="white", zorder=5)
        ax.annotate("Fern Tree", (fx, fy), xytext=(-8, 6), textcoords="offset points", fontsize=8,
                    fontstyle="italic", ha="right", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        title = ("Fires near Fern Tree (numbered dots; Fire History, theLIST)\n"
                 "EPSG:3577; coastline is 1:10 million, so only approximate here")
    else:
        small: list[dict] = []
        label_offset = {"riveaux-road": (-70, -22), "gell-river": (-62, -6),
                        "lake-mackenzie": (-100, -34), "mangana-road": (-40, 28)}
        for g in _coastline():
            for ring in _rings(to_albers(g)):
                x, y = ring.xy
                ax.fill(x, y, color="#f0f0f0", ec="#636363", lw=0.5)
        for r in summaries:
            b = json.loads((OUTPUTS / r["site"] / "boundary.geojson").read_text())
            g = to_albers(shape(b["features"][0]["geometry"]))
            for ring in _rings(g):
                x, y = ring.xy
                ax.fill(x, y, color=colours[r["site"]], alpha=0.85, lw=0)
            c = g.representative_point()
            if float(r["official_area_ha"]) < 10000:     # number the smaller fires, name the large ones
                ax.plot(c.x, c.y, "o", ms=6, mfc=colours[r["site"]], mec="black", mew=0.6)
                small.append(r)
                continue
            ax.annotate(f"{r['fire_name']}\n{r['ignition'][:4]}", (c.x, c.y), fontsize=7,
                        xytext=label_offset.get(r["site"], (6, 4)), textcoords="offset points",
                        arrowprops=(dict(arrowstyle="-", color="#636363", lw=0.6)
                                    if r["site"] in label_offset else None),
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=colours[r["site"]], alpha=0.85))
        if small:
            # number the small fires (any region) and list them in a key under the map
            key = []
            for i, r in enumerate(small, 1):
                g = to_albers(shape(json.loads((OUTPUTS / r["site"] / "boundary.geojson").read_text())
                                    ["features"][0]["geometry"]))
                c = g.representative_point()
                ax.annotate(str(i), (c.x, c.y), xytext=(4, 2), textcoords="offset points",
                            fontsize=6.5, fontweight="bold")
                key.append(f"{i}  {r['fire_name']} {r['ignition'][:4]}")
            small_key = "Fires under 10,000 ha (numbered dots):\n" + "\n".join(key)
        title = "Fire boundaries studied (Fire History, theLIST)\nEPSG:3577"
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=9)
    add_map_furniture(ax, list(ax.get_xlim()) + list(ax.get_ylim()))

    ax = fig.add_subplot(gs[1])
    labels = [f"{tag(r)}{r['fire_name']} ({r['ignition'][:4]})"
              + (" *" if str(r.get("severity_basis", "")).startswith("initial") else "")
              for r in summaries]
    left = np.zeros(len(summaries))
    for (name, _, _), col in zip(SEVERITY, SEVERITY_COLOURS):
        v = np.array([float(r[f"pct_{slugify(name)}"]) for r in summaries])
        ax.barh(labels, v, left=left, color=col, label=name, edgecolor="white", lw=0.5)
        left += v
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of the official boundary with clear data")
    ax.set_title("Burn severity mix, first summer after each fire\n(dNBR, Key & Benson classes; * = initial assessment)",
                 fontsize=9)
    ax.legend(fontsize=7, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.32))

    ax = fig.add_subplot(gs[2])
    max_years = 1
    n_lines = len(ax.lines)
    if local:
        # Small fires near Fern Tree are mostly low and moderate-low severity: plot those classes
        pair = ("Moderate-low", "Low")
        for r in summaries:
            rec = _read_csv(OUTPUTS / r["site"] / "recovery.csv")
            for k, ls in ((slugify(pair[0]), "-"), (slugify(pair[1]), "--")):
                pts = [(int(x["years_since"]), float(x[f"deficit_{k}"])) for x in rec
                       if x.get(f"deficit_{k}") not in ("", None)]
                if pts:
                    xs, ys = zip(*pts)
                    max_years = max(max_years, max(xs))
                    ax.plot(xs, ys, ls, marker="o", ms=3, color=colours[r["site"]])
                    # label each line with its map number
                    ax.annotate(str(number[r["site"]]), (xs[-1], ys[-1]), xytext=(4, 0),
                                textcoords="offset points", fontsize=7, va="center",
                                color=colours[r["site"]], fontweight="bold")
        rec_title = (f"{pair[0]} (solid) and {pair[1].lower()} (dashed) severity recovery\n"
                     "vs unburnt land; 0 = pre-fire; numbers as on the map")
        right_pad = 0.5
    else:
        ends = []
        k = slugify("High")
        for r in summaries:
            rec = _read_csv(OUTPUTS / r["site"] / "recovery.csv")
            pts = [(int(x["years_since"]), float(x[f"deficit_{k}"])) for x in rec
                   if x.get(f"deficit_{k}") not in ("", None)]
            if pts:
                xs, ys = zip(*pts)
                max_years = max(max_years, max(xs))
                ax.plot(xs, ys, "-", marker="o", ms=3, color=colours[r["site"]])
                if len(xs) >= 3:     # label only lines with a few summers; newer fires bunch at the start
                    ends.append([xs[-1], ys[-1], r["fire_name"], colours[r["site"]]])
        # direct labels at line ends, nudged apart vertically where they would overlap
        ends.sort(key=lambda e: (e[0], e[1]))
        for i in range(1, len(ends)):
            if ends[i][0] == ends[i - 1][0] and ends[i][1] - ends[i - 1][1] < 0.03:
                ends[i][1] = ends[i - 1][1] + 0.03
        for x, y, name, c in ends:
            ax.annotate(name, (x, y), xytext=(4, 0), textcoords="offset points", fontsize=6.5,
                        color="#333333", va="center")
        rec_title = "High-severity recovery vs unburnt land\n0 = pre-fire burn ratio"
        right_pad = 2.2      # room for the names
    if len(ax.lines) == n_lines:   # nothing to plot: say why rather than leave a blank panel
        ax.text(0.5, 0.5, "No severity class covers 2 ha (50 pixels) or more\nin any of these "
                "fires, so recovery cannot be tracked", transform=ax.transAxes, ha="center",
                va="center", fontsize=9, color="#333333", wrap=True)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("Summers since the fire (years)")
    ax.set_xticks(range(1, max_years + 1))
    ax.set_xlim(0.5, max_years + 0.5 + right_pad)
    ax.set_ylabel("NBR deficit vs before the fire (unitless)")
    ax.set_title(rec_title + "; colour = fire year", fontsize=9)
    ax.grid(alpha=0.3)
    sm = plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=year_cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.01)
    yrs = sorted(set(years))
    cb.set_ticks([0.88 * year_norm(y) for y in yrs])
    cb.set_ticklabels([str(y) for y in yrs])
    cb.ax.tick_params(labelsize=7)
    cb.set_label("Fire year (ignition)", fontsize=8)
    cb.ax.set_ylim(0, 0.88)
    fig.suptitle(f"{GROUPS[group][1]} from Sentinel-2 ({len(summaries)} fires, "
                 f"{summaries[0]['ignition'][:4]}-{summaries[-1]['ignition'][:4]}): where they burnt, "
                 "how severely, and how fast the vegetation is recovering", fontsize=11)
    fig.text(0.5, -0.02, DISCLAIMER + " Each fire: median Sentinel-2 composites before the fire and each summer "
             "after; dNBR corrected with unburnt land nearby. " + ATTRIBUTION +
             " Coastline: Natural Earth (public domain).", ha="center", va="top", fontsize=7,
             wrap=True)
    if small_key:   # outside the constrained layout, in the empty space under the map
        fig.text(0.012, 0.235, small_key, ha="left", va="top", fontsize=6.5, in_layout=False)
    fig.savefig(OUTPUTS / GROUPS[group][3], dpi=130, bbox_inches="tight")
    plt.close(fig)


LANDING_CSS = """
:root{--bg:#fbfaf7;--card:#ffffff;--ink:#1d1d1b;--muted:#5d5d58;--line:#e2e0da;--accent:#0072B2;
--chip:#f1efe9;--sev-low:#fed976;--sev-ml:#fd8d3c;--sev-mh:#e31a1c;--sev-hi:#800026;--sev-un:#d9d9d9}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#151514;--card:#1f1f1d;
--ink:#ecebe6;--muted:#a8a69f;--line:#34332f;--accent:#56B4E9;--chip:#2a2926}}
:root[data-theme="dark"]{--bg:#151514;--card:#1f1f1d;--ink:#ecebe6;--muted:#a8a69f;--line:#34332f;
--accent:#56B4E9;--chip:#2a2926}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.55 system-ui,-apple-system,
"Segoe UI",Roboto,sans-serif}
main{max-width:1100px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:clamp(1.6rem,4vw,2.3rem);line-height:1.15;margin:0 0 8px}
h2{font-size:1.25rem;margin:40px 0 12px}
p{margin:0 0 12px;max-width:72ch}
a{color:var(--accent)}
.lede{color:var(--muted);font-size:1.05rem}
.links{display:flex;flex-wrap:wrap;gap:8px;margin:20px 0 8px}
.links a{display:inline-block;padding:8px 14px;border:1px solid var(--line);border-radius:999px;
background:var(--card);text-decoration:none;font-weight:600}
figure{margin:16px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;
padding:8px}
figure img{width:100%;height:auto;display:block;border-radius:6px}
figcaption{font-size:.85rem;color:var(--muted);padding:8px 4px 2px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:16px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;overflow:hidden;
display:flex;flex-direction:column}
.card img{width:100%;aspect-ratio:16/10;object-fit:cover;object-position:right center;
display:block;border-bottom:1px solid var(--line);background:var(--chip)}
.card .pic{position:relative}
.inset{position:absolute;top:6px;left:6px;width:96px;height:auto;padding:3px;
background:var(--card);border:1px solid var(--line);
border-radius:6px;pointer-events:none}
.card .body{padding:14px 16px 16px;display:flex;flex-direction:column;gap:8px;flex:1}
.card h3{margin:0;font-size:1.1rem}
.when{color:var(--muted);font-size:.9rem;margin:0}
.region{font-size:.92rem;margin:0}
.stats{display:grid;grid-template-columns:1fr 1fr;gap:6px 12px;margin:4px 0;font-size:.9rem}
.stats div{display:flex;flex-direction:column}
.stats b{font-size:1.15rem;font-variant-numeric:tabular-nums}
.stats span{color:var(--muted);font-size:.8rem}
.bar{display:flex;height:10px;border-radius:5px;overflow:hidden;background:var(--chip)}
.bar i{display:block;height:100%}
.legend{display:flex;flex-wrap:wrap;gap:4px 12px;font-size:.8rem;color:var(--muted)}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;
vertical-align:-1px}
.warn{border:1px solid var(--sev-ml);border-left:6px solid var(--sev-ml);border-radius:8px;
background:var(--card);padding:10px 14px;margin:14px 0;max-width:72ch;font-size:.95rem}
.warn p{margin:0}
.note{font-size:.82rem;color:var(--muted);border-left:3px solid var(--sev-ml);padding-left:8px}
.figlinks{margin-top:auto;font-size:.9rem;display:flex;gap:12px;flex-wrap:wrap}
ul{padding-left:20px;max-width:72ch}
li{margin-bottom:6px}
.tabs>input{position:absolute;opacity:0;pointer-events:none}
.tablist{display:flex;gap:4px;margin-top:28px;border-bottom:1px solid var(--line)}
.tablist label{padding:10px 18px;border:1px solid transparent;border-radius:10px 10px 0 0;
margin-bottom:-1px;cursor:pointer;font-weight:600;color:var(--muted)}
.tabs>input:focus-visible~.tablist label{outline:2px solid var(--accent)}
.panel{display:none}
h3.sub{font-size:1.1rem;margin:28px 0 10px}
footer{margin-top:40px;font-size:.85rem;color:var(--muted);border-top:1px solid var(--line);
padding-top:16px}
"""


def _esc(text) -> str:
    import html
    return html.escape(str(text), quote=True)


def _num(v, fmt="{:,.0f}", suffix="") -> str:
    try:
        return fmt.format(float(v)) + suffix
    except (TypeError, ValueError):
        return "n/a"


INSET_W_PX = 96       # displayed width of the card locator inset (CSS px)
INSET_MIN_PX = 2.0    # draw the fire's own shape only if it spans at least this many px


def _inset_svgs(sites: list[str]) -> tuple[str, dict[str, str], float]:
    """Locator insets for the landing-page cards: one Tasmania outline (an SVG <symbol>, drawn
    once and reused, so it is identical on every card), plus per fire a horizontal and a vertical
    line crossing at the fire and, when it would span >= INSET_MIN_PX, the fire's official
    boundary filled in the same blue (#0072B2) as the boundary line on the card's maps. All in
    EPSG:3577 at one scale (km units), with a 100 km scale bar and a true-north arrow."""
    import shapely
    from shapely.geometry import box, shape
    clip = box(143.8, -43.7, 148.55, -39.55)   # main island, King Island and the Furneaux Group
    land = shapely.union_all([g.intersection(clip) for g in _coastline() if g.intersects(clip)])
    land = to_albers(land).simplify(800)
    x0, y0, x1, y1 = land.bounds
    pad = 6000
    x0, y0, x1, y1 = x0 - pad, y0 - pad, x1 + pad, y1 + 45000   # room at the top for N
    W, H = (x1 - x0) / 1000, (y1 - y0) / 1000                   # viewBox in km
    km_per_px = W / INSET_W_PX

    def path(geom) -> str:
        parts = []
        for ring in _rings(geom):
            xy = np.asarray(ring.coords)
            pts = " ".join(f"{(x - x0) / 1000:.1f},{(y1 - y) / 1000:.1f}" for x, y in xy)
            parts.append(f"M{pts}Z")
        return "".join(parts)

    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    ang = true_north_angle(cx, cy)
    nx, ny, L = W - 30, 66, 26                                  # N arrow base, length (km)
    sb_x, sb_y = 14, H - 12                                     # 100 km scale bar
    symbol = (f'<svg width="0" height="0" style="position:absolute" aria-hidden="true">'
              f'<symbol id="tas-inset" viewBox="0 0 {W:.1f} {H:.1f}">'
              f'<path d="{path(land)}" fill="var(--chip)" stroke="var(--muted)" stroke-width="2.2" '
              f'stroke-linejoin="round"/>'
              f'<g transform="rotate({ang:.1f} {nx} {ny})" fill="var(--ink)">'
              f'<line x1="{nx}" y1="{ny}" x2="{nx}" y2="{ny - L + 8}" stroke="var(--ink)" '
              f'stroke-width="3.5"/><path d="M{nx - 7},{ny - L + 10}L{nx},{ny - L}L{nx + 7},{ny - L + 10}Z"/>'
              f'<text x="{nx}" y="{ny - L - 4}" font-size="30" font-weight="700" '
              f'text-anchor="middle">N</text></g>'
              f'<path d="M{sb_x},{sb_y - 7}V{sb_y}H{sb_x + 100}V{sb_y - 7}" fill="none" '
              f'stroke="var(--ink)" stroke-width="3.5"/>'
              f'<text x="{sb_x + 50}" y="{sb_y - 12}" font-size="26" text-anchor="middle" '
              f'fill="var(--ink)">100 km</text>'
              f'</symbol></svg>')
    out = {}
    for site in sites:
        gj = json.loads((OUTPUTS / site / "boundary.geojson").read_text())
        g = to_albers(shape(gj["features"][0]["geometry"]))
        p = g.representative_point()
        fx, fy = (p.x - x0) / 1000, (y1 - p.y) / 1000
        bx0, by0, bx1, by1 = g.bounds
        span_px = max(bx1 - bx0, by1 - by0) / 1000 / km_per_px
        shape_svg = ""
        if span_px >= INSET_MIN_PX:
            shape_svg = (f'<path d="{path(g.simplify(1500))}" fill="#0072B2" stroke="#ffffff" '
                         f'stroke-width="1.5" fill-rule="evenodd"/>')
        what = ("its official boundary in blue, as on the maps" if shape_svg
                else "too small to show its boundary at this scale")
        out[site] = (f'<svg class="inset" viewBox="0 0 {W:.1f} {H:.1f}" role="img" '
                     f'aria-label="Locator map of Tasmania: the lines cross at the fire; {what}">'
                     f'<use href="#tas-inset"/>'
                     f'<g stroke="#D55E00" stroke-width="2.4" opacity="0.9">'
                     f'<line x1="0" y1="{fy:.1f}" x2="{W:.1f}" y2="{fy:.1f}"/>'
                     f'<line x1="{fx:.1f}" y1="0" x2="{fx:.1f}" y2="{H:.1f}"/></g>'
                     f'{shape_svg}</svg>')
    return symbol, out, INSET_MIN_PX * km_per_px


def _landing_page(groups: dict[str, list[dict]]) -> None:
    """outputs/tasmania/index.html: a static landing page for GitHub Pages (no external assets),
    with one tab per site group (CSS only: radio buttons, no scripts)."""
    sev_vars = ["--sev-un", "--sev-low", "--sev-ml", "--sev-mh", "--sev-hi"]
    summaries = list({r["site"]: r for rows in groups.values() for r in rows}.values())
    inset_symbol, insets, inset_min_km = _inset_svgs([r["site"] for r in summaries])
    card_html = {}
    for r in summaries:
        site = r["site"]
        ign = dt.date.fromisoformat(r["ignition"])
        out = dt.date.fromisoformat(r["declared_out"])
        bar = "".join(
            f'<i style="width:{float(r["pct_" + slugify(n)]):.1f}%;background:var({v})" '
            f'title="{n}: {float(r["pct_" + slugify(n)]):.1f}%"></i>'
            for (n, _, _), v in zip(SEVERITY, sev_vars))
        hi = float(r["pct_high"]) + float(r["pct_moderate-high"])
        rec = r.get("recovered_high_pct", "")
        initial = str(r.get("severity_basis", "")).startswith("initial")
        if initial:
            rec_txt, rec_lab = "n/a", "too soon to measure recovery"
        elif str(r["latest_summer"]) == str(r["severity_summer"]):
            rec_txt, rec_lab = "n/a", "only one summer since the fire so far"
        else:
            rec_txt = _num(rec, "{:.0f}", "%")
            rec_lab = f"high-severity recovery by summer {int(r['latest_summer']) - 1}-" \
                      f"{str(r['latest_summer'])[2:]}"
        notes = []
        if initial:
            notes.append("Initial assessment from images just after the fire: no full summer has "
                         "passed yet, so seasonal drying adds some false low severity.")
        base_txt = str(r.get("baseline", ""))
        if base_txt.startswith("just before the fire"):
            why = base_txt.split("; ", 1)[-1].rstrip(")") if "; " in base_txt else ""
            notes.append(f"Baseline is the images just before the fire ({why}), not the summer "
                         "before, so seasonal change may add some false low severity.")
        area = float(r["official_area_ha"])
        if area < 10:
            notes.append(f"Very small fire: {area:.1f} ha is only about {area / 0.04:.0f} pixels of "
                         "20 m, so the severity and recovery numbers are indicative only.")
        elif area < 1000:
            notes.append("Small fire: near the limit of 20 m pixels.")
        false_pct = 100 * float(r["burnt_outside_ha"]) / max(float(r["control_ha"]), 1)
        if false_pct > 8:
            notes.append(f"{false_pct:.0f}% of unburnt land nearby also shows a burn-like drop "
                         "(clearing, development, harvesting or drying): read low severity with "
                         "caution.")
        if "uncertain" in SITES.get(site, {}).get("groups", []):
            notes.insert(0, f"Weak or noisy result: {float(r['burnt_inside_pct']):.1f}% of the "
                         f"official boundary shows a burn (dNBR >= 0.10), against {false_pct:.1f}% "
                         f"of unburnt land nearby ({float(r['burnt_inside_ha']):.1f} ha burnt "
                         "inside). The fire may not be visible at 20 m.")
        note_html = "".join(f'<p class="note">{_esc(n)}</p>' for n in notes)
        card_html[site] = (f"""
<article class="card">
  <div class="pic"><a href="{_esc(site)}/dnbr_severity.png"><img loading="lazy" src="{_esc(site)}/dnbr_severity.png"
     alt="Burn severity map of the {_esc(r['fire_name'])} fire"></a>{insets[site]}</div>
  <div class="body">
    <h3>{_esc(r['fire_name'])} fire, {ign.year}</h3>
    <p class="when">Ignited {ign:%d %b %Y}, declared out {out:%d %b %Y} ({_esc(str(r['cause']).lower())})</p>
    <p class="region">{_esc(r['region'])}</p>
    <div class="stats">
      <div><b>{_num(r['official_area_ha'])} ha</b><span>official area</span></div>
      <div><b>{_num(r['burnt_inside_pct'], '{:.0f}', '%')}</b><span>of it burnt (dNBR &gt;= 0.10)</span></div>
      <div><b>{hi:.0f}%</b><span>moderate-high or high severity</span></div>
      <div><b>{rec_txt}</b><span>{_esc(rec_lab)}</span></div>
    </div>
    <div class="bar" role="img" aria-label="Severity mix">{bar}</div>
    {note_html}
    <div class="figlinks">
      <a href="{_esc(site)}/before_after.png">Before and after</a>
      <a href="{_esc(site)}/dnbr_severity.png">Severity</a>
      <a href="{_esc(site)}/recovery.png">Recovery</a>
      <a href="{_esc(site)}/summary.csv">Numbers (CSV)</a>
    </div>
  </div>
</article>""")
    legend = "".join(f'<span><i style="background:var({v})"></i>{n}</span>'
                     for (n, _, _), v in zip(SEVERITY, sev_vars))
    inputs, labels, panels = [], [], []
    for i, (g, rows) in enumerate(groups.items()):
        tab, heading, blurb, fig = GROUPS[g]
        inputs.append(f'<input type="radio" name="tab" id="tab-{g}"{" checked" if i == 0 else ""}>')
        labels.append(f'<label for="tab-{g}" role="tab">{_esc(tab)} ({len(rows)})</label>')
        cards = "".join(card_html[r["site"]] for r in sorted(rows, key=lambda r: r["ignition"],
                                                               reverse=True))
        panels.append(f"""<section class="panel" id="panel-{g}">
<h2>{_esc(heading)}: {_esc(blurb)}</h2>
<figure>
  <a href="{fig}"><img src="{fig}" alt="Overview of the {_esc(tab.lower())}: a map of the fires, the burn severity mix of each fire, and recovery of high-severity areas by years since the fire"></a>
  <figcaption>Left: the fires. Middle: severity mix inside each official boundary.
  Right: how fast {"high-severity areas regain" if g == "tasmania" else "moderate-low (solid) and low (dashed) severity areas regain"} their pre-fire
  burn ratio, relative to unburnt land nearby.</figcaption>
</figure>
<h3 class="sub">The fires, most recent first</h3>
<p class="note">{_esc(DISCLAIMER)}</p>
<div class="legend">{legend}</div>
<p class="when">Inset on each map: all of Tasmania at the same scale; the orange lines cross at
the fire, and the fire's official boundary is filled in the maps' boundary blue when it is big
enough to see (about {inset_min_km:.0f} km across or more).</p>
<div class="grid">{cards}
</div>
</section>""")
    tab_inputs, tab_labels, tab_panels = "\n".join(inputs), "".join(labels), "\n".join(panels)
    tab_css = "".join(f"#tab-{g}:checked~.tablist label[for=tab-{g}]{{background:var(--card);"
                      f"border-color:var(--line);border-bottom-color:var(--card);color:var(--ink)}}"
                      f"#tab-{g}:checked~#panel-{g}{{display:block}}" for g in groups)
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tasmanian Bushfires from Sentinel-2</title>
<meta name="description" content="Where recent Tasmanian bushfires burnt, how severely, and how fast the vegetation is recovering, from Sentinel-2 satellite images.">
<style>{LANDING_CSS}{tab_css}</style>
</head>
<body>
{inset_symbol}
<main>
<h1>Tasmanian bushfires from Sentinel-2</h1>
<p class="lede">Where {len(summaries)} Tasmanian bushfires burnt, how severely, and how fast the
vegetation is recovering, measured from free Sentinel-2 satellite images and checked against the
State's official fire boundaries.</p>
<div class="warn" role="note"><p><b>{_esc(DISCLAIMER.split(":")[0])}.</b>
{_esc(DISCLAIMER.split(": ", 1)[1][0].upper() + DISCLAIMER.split(": ", 1)[1][1:])} Severity classes
use generic thresholds and the numbers may contain errors. Don't use them for safety, planning
or management decisions. <a href="{TFS_URL}">{_esc(DISCLAIMER_OFFICIAL[:-1])}</a>.</p></div>
<div class="links">
  <a href="sites_map.html">Interactive map (all fires)</a>
  <a href="summary_all.csv">All numbers (CSV)</a>
</div>
<div class="tabs">
{tab_inputs}
<div class="tablist" role="tablist">{tab_labels}</div>
{tab_panels}
</div>

<h2>How to read this</h2>
<ul>
<li><b>Burn ratio (NBR)</b> compares near-infrared light (high from healthy leaves) with
short-wave infrared (high from bare, burnt ground). <b>dNBR</b> is the drop in NBR from the summer
before the fire to the first summer after it: the bigger the drop, the more severe the burn.</li>
<li>Severity classes use the generic Key &amp; Benson (2006) thresholds. They are not calibrated
against field plots in Tasmanian vegetation, so treat them as relative.</li>
<li><b>Recovery</b> is how much of the first-summer drop in NBR has been regained since, after
allowing for year-to-year change in unburnt land nearby. A recovered burn ratio is not a recovered
ecosystem: fire-sensitive communities such as pencil pine can be replaced by plants that look just
as green from space.</li>
<li>Thin cloud, smoke, forestry clearfelling and planned burns can also look like burns. Check the
image dates and the white (no clear data) areas in each figure.</li>
</ul>
<p>Method, caveats and how to run it for any Tasmanian fire: <a href="../../TASMANIA.md">TASMANIA.md</a>
and <a href="../../tas_fires.py">tas_fires.py</a>.</p>

<footer>
<p><b>{_esc(DISCLAIMER)}</b> <a href="{TFS_URL}">{_esc(DISCLAIMER_OFFICIAL[:-1])}</a>.</p>
<p>Sentinel-2: Digital Earth Australia Analysis Ready Data, Geoscience Australia, CC BY 4.0
(contains modified Copernicus Sentinel data). Fire boundaries: Fire History, theLIST, (c) State of
Tasmania, CC BY 3.0 AU. Coastline: Natural Earth (public domain). Map tiles: (c) OpenStreetMap
contributors. Built with Claude Code at the AEO26 workshop, Nipaluna/Hobart, November 2026.</p>
</footer>
</main>
</body>
</html>
"""
    (OUTPUTS / "index.html").write_text(page)


def _folium_map(summaries: list[dict], colours: dict) -> None:
    import folium
    m = folium.Map(location=[-42.0, 146.6], zoom_start=7, tiles="OpenStreetMap",
                   control_scale=True)
    for r in summaries:
        site = r["site"]
        fg = folium.FeatureGroup(name=f"{r['fire_name']} ({r['ignition'][:4]})")
        ov = json.loads((OUTPUTS / site / "severity_overlay.json").read_text())
        folium.raster_layers.ImageOverlay(str(OUTPUTS / site / "severity_overlay.png"),
                                          bounds=ov["bounds"],
                                          opacity=0.85).add_to(fg)
        html = (f"<b>{r['fire_name']}</b> (FIRE_ID {r['fire_id']})<br>{r['region']}<br>"
                f"Ignited {r['ignition']}, out {r['declared_out']}, cause {r['cause']}<br>"
                f"Official area {float(r['official_area_ha']):,.0f} ha; burnt by dNBR "
                f"{r['burnt_inside_pct']}%<br>High {r['pct_high']}%, moderate-high "
                f"{r['pct_moderate-high']}%<br>"
                f"High-severity recovery to {r['latest_summer']}: {r['recovered_high_pct']}%<br>"
                f"<a href='{site}/before_after.png' target='_blank'>before/after</a> | "
                f"<a href='{site}/dnbr_severity.png' target='_blank'>severity</a> | "
                f"<a href='{site}/recovery.png' target='_blank'>recovery</a>"
                f"<br><i style='font-size:11px'>{DISCLAIMER}</i>")
        folium.GeoJson(json.loads((OUTPUTS / site / "boundary.geojson").read_text()),
                       style_function=lambda _, c=colours[site]: dict(color=c, weight=2,
                                                                      fillOpacity=0.0),
                       tooltip=r["fire_name"], popup=folium.Popup(html, max_width=320)).add_to(fg)
        fg.add_to(m)
    legend = "".join(f"<div><span style='background:{c};width:12px;height:12px;display:inline-block;"
                     f"margin-right:4px'></span>{n}</div>"
                     for (n, _, _), c in zip(SEVERITY[1:], SEVERITY_COLOURS[1:]))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;bottom:30px;right:10px;z-index:9999;background:white;"
        "padding:6px;font-size:12px;border:1px solid #999'><b>Burn severity (dNBR)</b>"
        f"{legend}<div style='font-size:10px;max-width:220px'>{ATTRIBUTION}</div></div>"))
    m.get_root().html.add_child(folium.Element(
        "<div style='position:fixed;top:10px;left:60px;right:60px;z-index:9999;max-width:640px;"
        "margin:0 auto;background:#fff8e6;border:1px solid #e69f00;border-left:6px solid #e69f00;"
        f"padding:6px 10px;font:13px system-ui,sans-serif'><b>{DISCLAIMER.split(':')[0]}.</b> "
        f"{DISCLAIMER.split(': ', 1)[1].capitalize()} <a href='{TFS_URL}' target='_blank'>"
        f"{DISCLAIMER_OFFICIAL[:-1]}</a>.</div>"))
    folium.Marker([FERN_TREE[1], FERN_TREE[0]], tooltip="Fern Tree",
                  icon=folium.Icon(color="black", icon="star")).add_to(m)
    folium.LayerControl(collapsed=True).add_to(m)
    m.save(str(OUTPUTS / "sites_map.html"))


# ----------------------------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------------------------

FERN_TREE = (147.263, -42.923)   # lon, lat: distance reference for --candidates


def region_of(lon: float, lat: float) -> str:
    """Rough Tasmanian region from a fire's centre (for triage only)."""
    if lat > -40.6 or lon < 144.3:
        return "Bass Strait islands"
    if lat > -41.6 and lon < 146.3:
        return "North-west"
    if lat > -41.75 and lon >= 147.3:
        return "North-east"
    if lon >= 147.6:
        return "East"
    if -42.35 < lat < -41.55 and 146.0 < lon < 147.0:
        return "Central Plateau"
    if lat < -42.3 and lon < 146.7:
        return "South-west"
    if lat < -42.55 and lon >= 146.7:
        return "South / Hobart"
    if lon < 145.9:
        return "West coast"
    return "Midlands / Derwent"


def candidates(min_ha: float = 100, today: dt.date | None = None) -> list[dict]:
    """Every Fire History bushfire since 2015-16 >= min_ha not yet a site, tagged for triage."""
    today = today or dt.date.today()
    last = latest_summer(today)
    done = {v["fire_id"] for v in SITES.values()}
    rows, offset = [], 0
    while True:
        d = _get_json(FIRE_LAYER + "/query", dict(
            where=f"FIRE_TYPE='Bushfire' AND IGN_SEASON>={MIN_SEASON} AND FIRE_AR_HA>={min_ha}",
            outFields=FIRE_FIELDS, returnGeometry="true", outSR=4326, maxAllowableOffset=0.01,
            orderByFields="FIRE_AR_HA DESC", resultOffset=offset, resultRecordCount=1000, f="json"))
        for f in d.get("features", []):
            r = _fire_record(f["attributes"])
            if r["fire_id"] in done or r["ignition"] is None:
                continue
            pts = [p for ring in f["geometry"]["rings"] for p in ring]
            lon = (min(p[0] for p in pts) + max(p[0] for p in pts)) / 2
            lat = (min(p[1] for p in pts) + max(p[1] for p in pts)) / 2
            end = r["end"] or r["ignition"]
            km = math.hypot((lon - FERN_TREE[0]) * 81.6, (lat - FERN_TREE[1]) * 111.0)
            summer_base = summer_window(pre_summer_year(r["ignition"]))[0] >= dt.date(2015, 7, 1)
            first = first_summer_after(end)
            rows.append(dict(
                fire_id=r["fire_id"], name=r["name"], area_ha=round(r["area_ha"]),
                ignition=r["ignition"].isoformat(), declared_out=end.isoformat(), cause=r["cause"],
                region=region_of(lon, lat), km_from_fern_tree=round(km), lon=round(lon, 4),
                lat=round(lat, 4), baseline="summer before" if summer_base else "pre-fire (fallback)",
                assessment="extended" if first <= last else "initial",
                recovery_summers=max(0, last - first + 1) if first <= last else 0,
                pixels_20m=round(r["area_ha"] / 0.04)))
        if not d.get("exceededTransferLimit"):
            break
        offset += 1000
    return sorted(rows, key=lambda x: -x["area_ha"])


def _print_fires(fires: list[dict]) -> None:
    for f in fires:
        print(f"  --fire-id {f['fire_id']:<6} {f['name'][:34]:34} {f['area_ha']:>9,.0f} ha  "
              f"ignited {f['ignition']}  out {f['end']}  ({f['cause']})")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__.split("Examples")[1])
    p.add_argument("--list", action="store_true", help="list the built-in sites")
    p.add_argument("--search", metavar="TEXT", help="search Fire History bushfires by name")
    p.add_argument("--site", action="append", choices=sorted(SITES), help="built-in site (repeatable)")
    p.add_argument("--fire-id", type=int, action="append", help="any FIRE_ID from Fire History")
    p.add_argument("--lon", type=float, help="longitude: analyse the latest fire covering this point")
    p.add_argument("--lat", type=float, help="latitude (with --lon)")
    p.add_argument("--group", choices=sorted(GROUPS), action="append",
                   help="every built-in site in a group (repeatable), then the overview")
    p.add_argument("--all", action="store_true", help="every built-in site, then the overview")
    p.add_argument("--overview", action="store_true", help="overview figure and interactive map")
    p.add_argument("--candidates", action="store_true",
                   help="triage every other bushfire since 2015-16 (>= 100 ha) to "
                        "outputs/tasmania/candidates.csv")
    p.add_argument("--max-dates", type=int, default=4, help="images per composite (default 4)")
    p.add_argument("--refresh", action="store_true", help="ignore cached composites")
    a = p.parse_args(argv)

    if a.list:
        print("Built-in sites (python tas_fires.py --site KEY):")
        for k, v in SITES.items():
            rec, _ = get_fire(v["fire_id"])
            print(f"  {k:24} {rec['name']:22} ignited {rec['ignition']}  {rec['area_ha']:>8,.0f} ha"
                  f"  {v['region']}")
        return 0
    if a.candidates:
        rows = candidates()
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        _write_csv(OUTPUTS / "candidates.csv", rows)
        print(f"{len(rows)} candidate fires written to outputs/tasmania/candidates.csv")
        return 0
    if a.search is not None:
        _print_fires(search_fires(a.search))
        return 0
    jobs = []
    keys = sorted(SITES) if a.all else list(a.site or [])
    for grp in a.group or []:
        keys += [k for k, v in SITES.items() if grp in v["groups"] and k not in keys]
    for k in keys:
        jobs.append((k, SITES[k]["fire_id"], SITES[k]["region"]))
    for fid in a.fire_id or []:
        rec, _ = get_fire(fid)
        jobs.append((f"{slugify(rec['name'])}-{rec['ignition'].year}", fid, ""))
    if a.lon is not None or a.lat is not None:
        if a.lon is None or a.lat is None:
            p.error("--lon and --lat go together")
        fires = search_fires(lon=a.lon, lat=a.lat)
        if not fires:
            raise SystemExit(f"No bushfire since 2015-16 covers {a.lon}, {a.lat} in Fire History")
        print(f"Bushfires covering {a.lon}, {a.lat} (latest first):")
        _print_fires(fires)
        f = fires[0]
        jobs.append((f"{slugify(f['name'])}-{f['ignition'].year}", f["fire_id"], ""))
    if not jobs and not a.overview:
        p.print_help()
        return 1
    for key, fid, region in jobs:
        analyse_fire(key, fid, region, max_dates=a.max_dates, refresh=a.refresh)
    if a.overview or a.all or a.group:
        overview()
    return 0


if __name__ == "__main__":
    sys.exit(main())
