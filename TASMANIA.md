# Tasmanian bushfire impacts from Sentinel-2

`tas_fires.py` is a small, self-contained demonstration. You choose a recent Tasmanian bushfire, and
it shows **where the fire burnt, how severely, and how fast the vegetation has recovered**, using
free Sentinel-2 satellite images and the State's official fire boundaries.

## Run it

```bash
pip install -r requirements.txt pyproj shapely matplotlib-scalebar requests pyshp
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
| `python tas_fires.py --all` | Every built-in site, then the overview |
| `python tas_fires.py --overview` | Rebuild the overview and interactive map from finished sites |

The first run of a site downloads about 0.5-2 GB of image pieces and takes about 10-30 minutes.
Composites are cached in `downloads/tasmania/<site>/` (not committed), so later runs are fast.

## Built-in sites

| Key | Fire | Ignited | Official area | Landscape |
|---|---|---|---|---|
| `lake-mackenzie` | Lake Mackenzie Road | Jan 2016 | 24,711 ha | Central Plateau, World Heritage Area: alpine and subalpine, pencil pine |
| `gell-river` | Gell River | Dec 2018 | 35,037 ha | Southwest wilderness: buttongrass, scrub and forest |
| `riveaux-road` | Riveaux Road | Jan 2019 | 63,895 ha | Southern forests west of Geeveston: tall wet eucalypt forest |
| `great-pine-tier` | Great Pine Tier | Jan 2019 | 51,228 ha | Central Plateau near Miena: subalpine woodland, grassland, heath |
| `mangana-road` | Mangana Road | Dec 2019 | 14,760 ha | North-east, Fingal Valley: dry eucalypt forest and farmland |
| `mount-donaldson` | Mount Donaldson | Feb 2025 | 57,949 ha | West coast near Corinna (takayna/Tarkine) |

These are among the largest bushfires since Sentinel-2 began (2015-16), and they cover the main
Tasmanian landscapes and three fire seasons.

## Outputs

For each site, in `outputs/tasmania/<site>/`:

| File | What it shows |
|---|---|
| `before_after.png` | False colour (SWIR2, NIR, red) before the fire, just after it, the first summer after, and the latest summer, with one colour stretch for all panels and the official boundary |
| `dnbr_severity.png` | dNBR (the change in burn ratio) and burn-severity classes, with the official boundary |
| `recovery.png` | The burn-ratio deficit and greenness (NDVI) each summer, by severity class, against unburnt land nearby |
| `summary.csv` | Areas by severity class (ha), agreement with the official boundary, and the image dates used |
| `recovery.csv` | The numbers behind `recovery.png` |

Across all sites, in `outputs/tasmania/`: `overview.png` (map, severity mix, recovery curves),
`summary_all.csv`, and `sites_map.html`, an interactive map with the severity overlays. Open the
HTML file from a local copy of the repo: the links inside it are relative.

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
4. **Burn severity:** NBR = (NIR - SWIR2) / (NIR + SWIR2), with B8A and B12. dNBR = NBR before minus
   NBR in the **first summer after the fire**. That is the same season as the before image, about a
   year on (an "extended assessment"). The median dNBR of unburnt land nearby (outside every fire
   recorded in Fire History over the period) is subtracted, to remove differences in season,
   sun angle and processing. The classes use the Key & Benson (2006) thresholds: unburnt < 0.10,
   low 0.10-0.27, moderate-low 0.27-0.44, moderate-high 0.44-0.66, high >= 0.66.
5. **Recovery:** for each summer, the median NBR deficit (NBR before minus NBR that summer) per
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
- Cloud is common in western Tasmania. Check the dates and the % covered in `summary.csv`, and the
  white (no data) areas in the figures.

## Data and attribution

- Sentinel-2: Digital Earth Australia, Geoscience Australia, CC BY 4.0 (contains modified
  Copernicus Sentinel data).
- Fire boundaries: Fire History, theLIST, (c) State of Tasmania, CC BY 3.0 AU.
- Coastline (overview map): Natural Earth, public domain.
- Interactive map tiles: (c) OpenStreetMap contributors.
