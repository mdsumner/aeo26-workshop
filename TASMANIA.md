# Tasmanian bushfire impacts from Sentinel-2

`tas_fires.py` is a small, self-contained demonstration. You choose a recent Tasmanian bushfire, and
it shows **where the fire burnt, how severely, and how fast the vegetation has recovered**, using
free Sentinel-2 satellite images and the State's official fire boundaries.

## Run it

All raster reading, warping and writing goes through GDAL's own Python bindings (`osgeo.gdal`).
The image search is plain JSON from DEA's STAC API, fetched with `requests`: no STAC library is
needed. The bindings must match an installed GDAL (3.8 or later), so install them with GDAL, not
from PyPI on their own:

| Where | GDAL bindings |
|---|---|
| Docker | `ghcr.io/hypertidy/gdal-r-python` already has them (GDAL 3.13), and R's `gdalraster` and `terra` too |
| Ubuntu or Debian | `sudo apt install python3-gdal` (use the system `python3`) |
| conda, mamba or pixi | the `gdal` package from conda-forge |

Then the rest, with pip:

```bash
pip install numpy shapely matplotlib matplotlib-scalebar requests pyshp pillow folium pytest
```

```bash
python tas_fires.py --list
```

```bash
python tas_fires.py --site riveaux-road
```

Other ways to choose a fire:

| Command | What it does |
|---|---|
| `python tas_fires.py --site KEY` | One of the built-in sites below (repeatable) |
| `python tas_fires.py --search "pine tier"` | Search Fire History bushfires (2015-16 onwards) by name |
| `python tas_fires.py --fire-id 10137` | Any fire by its Fire History `FIRE_ID` |
| `python tas_fires.py --lon 146.7 --lat -43.1` | The most recent bushfire that covers this point |
| `python tas_fires.py --group fern-tree` | Every built-in site in a group (`tasmania` or `fern-tree`), then the overview |
| `python tas_fires.py --all` | Every built-in site, then the overview |
| `python tas_fires.py --overview` | Rebuild the overview and interactive map from finished sites |
| `python tas_fires.py --candidates` | Triage every other bushfire since 2015-16 (100 ha or more) into `outputs/tasmania/candidates.csv`: region, distance from Fern Tree, which baseline and assessment it would get, and pixel count |

The first run of a site downloads about 0.5-2 GB of image pieces and takes about 10-30 minutes.
Composites are cached in `downloads/tasmania/<site>/` (not committed), so later runs are fast.

## Built-in sites

The built-in sites are listed in `sites.csv` (key, `FIRE_ID`, groups, region), which `tas_fires.py`
and `tas_fires.R` both read. To add a site, add a row. Each group is a tab on the landing page.

### Major fires (group `tasmania`)

| Key | Fire | Ignited | Official area | Landscape |
|---|---|---|---|---|
| `lake-mackenzie` | Lake Mackenzie Road | Jan 2016 | 24,711 ha | Central Plateau, World Heritage Area: alpine and subalpine, pencil pine |
| `gell-river` | Gell River | Dec 2018 | 35,037 ha | Southwest wilderness: buttongrass, scrub and forest |
| `riveaux-road` | Riveaux Road | Jan 2019 | 63,895 ha | Southern forests west of Geeveston: tall wet eucalypt forest |
| `great-pine-tier` | Great Pine Tier | Jan 2019 | 51,228 ha | Central Plateau near Miena: subalpine woodland, grassland, heath |
| `mangana-road` | Mangana Road | Dec 2019 | 14,760 ha | North-east, Fingal Valley: dry eucalypt forest and farmland |
| `mount-donaldson` | Mount Donaldson | Feb 2025 | 57,949 ha | West coast near Corinna (takayna/Tarkine) |
| `proctors-road` | Proctors Road | Oct 2017 | 226 ha | Hobart's southern hills, about 3.5 km from Fern Tree |
| `margate-plain` | Margate Plain | Jan 2025 | 484 ha | Hills west of Margate, south of Hobart |
| `glenlusk` | Molesworth Road | Dec 2025 | 36 ha | Glenlusk, north-west of Hobart near Molesworth |
| `nelsons-tier` | Nelsons Tier | Oct 2025 | 1,286 ha | About 40 km north-east of Hobart |
| `labillardiere-peninsula` | Labillardiere Peninsula | Oct 2025 | 825 ha | South Bruny National Park, Bruny Island |
| `pelham-road` | Pelham Road | Dec 2019 | 2,115 ha | About 40 km north of Hobart near Pelham (lightning) |
| `thunder-and-lightning-bay` | Thunder And Lightning Bay | Dec 2016 | 31,770 ha | Flinders Island, Bass Strait |
| `dolphin-sands-2025` | Dolphin Sands Road | Dec 2025 | 691 ha | Dolphin Sands near Swansea, east coast (initial assessment) |
| `dolphin-sands-2023` | Dolphin Sands Road | Nov 2023 | 250 ha | The same spit, 5 km further east (no overlap with 2025) |
| `mcneills-road-2026` | Mcneills Road | Apr 2026 | 1,064 ha | Next to Dolphin Sands, autumn 2026 (initial; winter after-images) |
| `tasman-highway-2025` | Tasman Highway | Dec 2025 | 1,170 ha | North-east: the largest fire of December 2025 (initial assessment) |

The last four were picked with `--candidates` (see below). The first six are among the largest bushfires since Sentinel-2 began (2015-16), and they cover
the main Tasmanian landscapes. The last three are small fires on Hobart's urban edge. They are at
the limit of what 20 m pixels can show (36 ha is about 900 pixels), but they are the fires closest
to home. `glenlusk` was declared out in December 2025, so it has no full summer since the fire yet.
Its severity is an *initial assessment* from the images just after the fire, and the tool switches
to the first summer automatically once summer 2026-27 has passed.

### Near Fern Tree (group `fern-tree`)

Bushfires since 2015-16 within about 25 km of Fern Tree, on kunanyi/Mount Wellington: every fire of
5 ha or more within 25 km, plus fires of 2 ha or more within 10 km (distance from the centre of the
fire's boundary). Proctors Road, Margate Plain and
Glenlusk are in this group too. Most of the 99 bushfires recorded within 15 km are campfires under
1 ha, which 20 m pixels cannot measure (one pixel is 0.04 ha). Even the fires listed here are small:
2 ha is about 50 pixels, so their severity and recovery numbers are indicative only, and the landing
page says so on each card under 10 ha.

| Key | Fire | Ignited | Official area | From Fern Tree |
|---|---|---|---|---|
| `turn-creek-road-2016` | Turn Creek Road | Jan 2016 | 7.2 ha | 14 km west |
| `tianna-road-2017` | Tianna Road | Feb 2017 | 22.3 ha | 11 km north-east |
| `magnolia-road-2017` | Magnolia Road | Oct 2017 | 6.9 ha | 16 km north-east |
| `holmfield-avenue-2017` | Holmfield Avenue | Nov 2017 | 13.0 ha | 15 km east |
| `richmond-road-2019` | Richmond Road | Jan 2019 | 15.1 ha | 20 km north-east |
| `brogo-way-2019` | Brogo Way | Jan 2019 | 5.8 ha | 16 km east |
| `slatterys-road-2019` | Slatterys Road, Electrona | Oct 2019 | 6.2 ha | 14 km south |
| `mount-hull-road-2020` | Mount Hull Road | Nov 2020 | 9.9 ha | 9 km north-west |
| `downhams-road-2023-may` | Downhams Road | May 2023 | 9.3 ha | 15 km north-east |
| `downhams-road-2023-oct` | Downhams Road | Oct 2023 | 6.8 ha | 15 km north-east |
| `mount-hull-road-2024` | Mount Hull Road | Oct 2024 | 2.0 ha | 10 km north |
| `snowy-view-heights-2024` | Snowy View Heights | Oct 2024 | 53.6 ha | 18 km west |
| `glen-dhu-road-2026` | Glen Dhu Road | Feb 2026 | 5.6 ha | 18 km north-west |
| `hydehurst-road-2016` | Hydehurst Road | Dec 2016 | 8.9 ha | 19 km north-west |
| `glenfern-road-2019` | Glenfern Road, Glenfern | Oct 2019 | 210.8 ha | 25 km north-west |
| `cloverside-road-2020` | Cloverside Road | Nov 2020 | 9.1 ha | 18 km west |
| `back-river-road-2023` | Back River Road | Feb 2023 | 10.2 ha | 25 km north-west |
| `boyer-road-2026` | Boyer Road | Jan 2026 | 8.3 ha | 21 km north |

The last five were added in November 2026, when the radius went from 20 to 25 km, with four more
that went to Uncertain results (below). Hydehurst Road, Cloverside Road, Rockingham Drive and
Jeffreys Track already met the 20 km rule but had been missed.

### Uncertain results (group `uncertain`)

Fires near Fern Tree that meet the rules above, but whose results are weak or noisy, so they
have their own tab: the share of the boundary that shows a burn is close to, or below, the share of
unburnt land nearby that shows the same drop, or the burn is only a speckle of low severity. Each
card gives the two percentages.

| Key | Fire | Ignited | Official area | Why it is uncertain |
|---|---|---|---|---|
| `peter-murrell-2017` | Peter Murrell NR, Blackmans Bay | Nov 2017 | 4.9 ha | Less of the boundary shows a burn than of unburnt land nearby |
| `wyre-forest-road-2021` | Wyre Forest Road, Molesworth | Mar 2021 | 16.6 ha | Only about 10% of the boundary shows a burn |
| `channel-highway-2024` | Channel Highway | Feb 2024 | 2.4 ha | About 0.5 ha of low-severity speckle |
| `rockingham-drive-2015` | Rockingham Drive | Nov 2015 | 26.7 ha | Less of the boundary shows a burn (4.9%) than of unburnt land nearby (13.3%) |
| `ethan-court-2016` | Ethan Court | Apr 2016 | 7.0 ha | A thin strip inside a wide speckle of low severity (36.7% inside, 18.3% nearby) |
| `jeffreys-track-2019` | Jeffreys Track | Dec 2019 | 20.7 ha | Little inside the boundary (21.6%) beyond what unburnt land nearby shows (15.5%) |
| `shark-point-road-2019` | Shark Point Road | Dec 2019 | 14.0 ha | On the estuary edge; less shows a burn inside (11.3%) than nearby (14.4%) |

## Outputs

For each site, in `outputs/tasmania/<site>/`:

| File | What it shows |
|---|---|
| `before_after.png` | False colour (SWIR2, NIR, red) before the fire, just after it, the first summer after, and the latest summer, with one colour stretch for all panels and the official boundary |
| `dnbr_severity.png` | dNBR (the change in burn ratio) and burn-severity classes, with the official boundary |
| `recovery.png` | The burn-ratio deficit and greenness (NDVI) each summer, by severity class, against unburnt land nearby |
| `summary.csv` | Areas by severity class (ha), agreement with the official boundary, and the image dates used |
| `recovery.csv` | The numbers behind `recovery.png` |

Across all sites, in `outputs/tasmania/`:

| File | What it is |
|---|---|
| `index.html` | Landing page with a tab per group of sites, and one card per fire (most recent first) with links to every figure. An **All** tab lists every fire once and sorts by date, distance to a place (or your own location, worked out in the browser), size, or certainty (% of the boundary burnt minus % of unburnt land nearby with the same drop). Built for GitHub Pages, with no external files; one small inline script does the sorting |
| `sites_map.html` | Interactive map with the fire boundaries and severity overlays (embedded), and popups linking to the figures |
| `overview.png` | Major fires: map of the sites, the severity mix of each fire, and recovery curves |
| `overview_fern_tree.png` | The same for the fires near Fern Tree, on a zoomed map with numbered fires |
| `overview_uncertain.png` | The same for the fires with uncertain results |
| `summary_all.csv` | One row per fire |

`python tas_fires.py --overview` rebuilds all of these from the finished sites. The links are relative,
so the pages work on GitHub Pages (`/outputs/tasmania/`) or from a local copy of the repo.

## Method

1. **Fire boundary and dates:** Fire History (Tasmania Fire Service and Parks and Wildlife Service)
   via theLIST. The study area is the fire's extent plus 3 km, on a 20 m grid in **EPSG:3577**
   (Australian Albers, equal-area).
2. **Images:** Digital Earth Australia Sentinel-2 Analysis Ready Data (NBART surface reflectance,
   Sentinel-2A/B/C), read straight from the public `dea-public-data` S3 bucket. Only the part
   covering the study area is read. Values are reflectance x 10,000, with nodata -999. Clouds,
   shadows, snow and water are removed with DEA's `fmask` layer.
3. **Composites:** the per-pixel median of the clearest dates in each window, adding dates until at
   least 90% of the study area is covered:
   - *before:* the 120 days before ignition (365 days if that is too cloudy)
   - *just after:* the 120 days after the fire was declared out
   - *each summer* (1 Dec to 31 Mar) after the fire. Later summers are read at 40 m to save downloads.
4. **Burn severity:** NBR = (NIR - SWIR2) / (NIR + SWIR2), with B8A and B12. dNBR = NBR in the
   **last full summer before the fire** minus NBR in the **first summer after it** (an "extended
   assessment", about a year on). Both are summers, so the grass and crops that dry out from spring
   to summer don't show up as burns. That problem added 16% false burn around Great Pine Tier when
   the spring images just before the fire were used. If Sentinel-2 has no summer before the fire
   (Lake Mackenzie, January 2016), the images just before the fire are used instead, and the figure
   says so. The median dNBR of unburnt land nearby (outside every fire
   recorded in Fire History over the period) is subtracted, to remove differences in season,
   sun angle and processing. The classes use the Key & Benson (2006) thresholds: unburnt < 0.10,
   low 0.10-0.27, moderate-low 0.27-0.44, moderate-high 0.44-0.66, high >= 0.66.
5. **Recovery:** for each summer, the median NBR deficit (NBR in the summer before the fire minus NBR that summer) per
   severity class, minus the same for unburnt land. Zero means back to the pre-fire burn ratio.

## Caveats

- The Key & Benson thresholds are generic (developed in the USA). They are not calibrated against
  field plots in Tasmanian vegetation, so treat the classes as relative, not absolute.
- Severity is measured about a year after the fire. Fast-recovering vegetation, such as buttongrass
  moorland and grassland, may already have regreened, so its severity is underestimated.
- A recovered burn ratio is not a recovered ecosystem. Fire-sensitive communities, such as pencil
  pine and King Billy pine on the Central Plateau, can be replaced by other plants that look just as
  green from space.
- Official boundaries include unburnt patches, and they are drawn at different scales for different
  fires. The "% of boundary burnt" figure reflects both of these.
- `fmask` misses some thin cloud, haze and smoke. These lower NBR, so they can show up as false
  "low" severity outside the boundary (for example, the pale patches south of the Riveaux Road
  fire in summer 2019-20). Forestry clearfelling and planned burns also look like burns from space.
- Coastal sites: coverage is measured over land only. A pixel that fmask labels water on any date
  is treated as sea.
- Cloud is common in western Tasmania. Check the dates and the % covered in `summary.csv`, and the
  white (no data) areas in the figures.

## Data and attribution

- Sentinel-2: Digital Earth Australia, Geoscience Australia, CC BY 4.0 (contains modified
  Copernicus Sentinel data).
- Fire boundaries: Fire History, theLIST, (c) State of Tasmania, CC BY 3.0 AU.
- Coastline (overview map): Natural Earth, public domain.
- Interactive map tiles: (c) OpenStreetMap contributors.
