#!/usr/bin/env Rscript
# Tasmanian bushfire impacts from Sentinel-2: the analysis in R, with gdalraster.
#
# This is a line-by-line twin of the analysis in tas_fires.py (fire boundary, composites, dNBR,
# severity areas, recovery), so the two can be compared. It writes the same summary.csv and
# recovery.csv, to outputs/tasmania-r/<site>/, from its own image cache in downloads/tasmania-r/.
# It does not draw the figures or the web pages: tas_fires.py does that.
#
# All raster and geometry work goes through GDAL (gdalraster). The fire boundaries come from
# theLIST's ArcGIS REST service and the image list from DEA's STAC API: both are plain JSON over
# HTTPS, read with curl and jsonlite.
#
# Usage
#   Rscript tas_fires.R --site glenlusk
#   Rscript tas_fires.R --site glenlusk --compare   # then compare with the Python outputs
#   Rscript tas_fires.R --compare-only --site glenlusk
#
# Packages: gdalraster, jsonlite, curl; matrixStats if installed (much faster medians).

suppressPackageStartupMessages({
  library(gdalraster)
  library(jsonlite)
})

# GDAL settings for reading cloud-optimised GeoTIFFs over HTTPS (as in tas_fires.py)
for (kv in list(c("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR"),
                c("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif"),
                c("GDAL_HTTP_MULTIRANGE", "YES"),
                c("GDAL_HTTP_MERGE_CONSECUTIVE_RANGES", "YES"),
                c("GDAL_HTTP_MAX_RETRY", "4"),
                c("GDAL_HTTP_RETRY_DELAY", "2"),
                c("VSI_CACHE", "TRUE"))) {
  if (!nzchar(Sys.getenv(kv[1]))) do.call(Sys.setenv, setNames(list(kv[2]), kv[1]))
}

script_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
ROOT <- script_dir()
DOWNLOADS <- file.path(ROOT, "downloads", "tasmania-r")
OUTPUTS <- file.path(ROOT, "outputs", "tasmania-r")

FIRE_LAYER <- paste0("https://services.thelist.tas.gov.au/arcgis/rest/services/Public/",
                     "EmergencyManagementPublic/MapServer/1")
DEA_STAC_SEARCH <- "https://explorer.dea.ga.gov.au/stac/search"
DEA_COLLECTIONS <- c("ga_s2am_ard_3", "ga_s2bm_ard_3", "ga_s2cm_ard_3")

CRS <- "EPSG:3577"   # GDA94 / Australian Albers, equal-area, metres
RES <- 20            # metres
COARSE <- 2L         # recovery summers are read at RES * COARSE
BUFFER_M <- 3000     # study area = fire extent plus this margin
NODATA <- -999       # DEA ARD nodata
SCALE <- 10000       # DEA ARD reflectance scale
FMASK_CLEAR <- 1     # fmask: 0 nodata, 1 clear, 2 cloud, 3 shadow, 4 snow, 5 water
BANDS <- c(red = "nbart_red", nir = "nbart_nir_2", swir1 = "nbart_swir_2", swir2 = "nbart_swir_3")
FMASK <- "oa_fmask"
FIRE_FIELDS <- "FIRE_ID,FIRE_NAME,FIRE_AR_HA,IGN_DATE,INCD_END,IGN_SEASON,IGN_CAUSE1,FIRE_TYPE"

# Key & Benson (2006) dNBR thresholds; class index 0..4
SEVERITY <- data.frame(name = c("Unburnt", "Low", "Moderate-low", "Moderate-high", "High"),
                       lo = c(-Inf, 0.10, 0.27, 0.44, 0.66),
                       hi = c(0.10, 0.27, 0.44, 0.66, Inf))

SITES <- list(
  "lake-mackenzie" = list(fire_id = 8075, region = paste(
    "Central Plateau, Tasmanian Wilderness World Heritage Area: alpine and subalpine",
    "vegetation, including fire-sensitive pencil pine")),
  "riveaux-road" = list(fire_id = 10137, region = paste(
    "Southern forests west of Geeveston and the Huon Valley: tall wet eucalypt forest",
    "next to Southwest National Park")),
  "gell-river" = list(fire_id = 10105, region =
    "Southwest wilderness (World Heritage Area): buttongrass moorland, scrub and forest"),
  "great-pine-tier" = list(fire_id = 10148, region =
    "Central Plateau near Miena: subalpine woodland, grassland and heath"),
  "mangana-road" = list(fire_id = 10768, region =
    "North-east, Fingal Valley: dry eucalypt forest and farmland"),
  "mount-donaldson" = list(fire_id = 14017, region = paste(
    "West coast near Corinna and the Pieman River (takayna/Tarkine): buttongrass,",
    "scrub and rainforest")),
  "glenlusk" = list(fire_id = 14725, region =
    "Glenlusk, north-west of Hobart near Molesworth: bushland at the urban edge"),
  "proctors-road" = list(fire_id = 9082, region =
    "Hobart's southern hills, about 3.5 km from Fern Tree: bushland at the urban edge"),
  "margate-plain" = list(fire_id = 14005, region =
    "Hills west of Margate, south of Hobart: bushland and forest")
)

# ---------------------------------------------------------------------------------------------
# Small pure functions (the same as in tas_fires.py)
# ---------------------------------------------------------------------------------------------

to_reflectance <- function(dn) {
  out <- dn / SCALE
  out[dn == NODATA | dn < 0] <- NaN
  out
}

norm_diff <- function(a, b) {
  out <- (a - b) / (a + b)
  out[!is.finite(out)] <- NaN
  out
}

nbr <- function(c) norm_diff(c$nir, c$swir2)
ndvi <- function(c) norm_diff(c$nir, c$red)

classify_dnbr <- function(dnbr) {
  edges <- SEVERITY$hi[-nrow(SEVERITY)]
  x <- dnbr
  x[!is.finite(x)] <- 0
  cls <- findInterval(x, edges)          # same as numpy.digitize: edges are lower-inclusive
  cls[!is.finite(dnbr)] <- -1L
  dim(cls) <- dim(dnbr)
  cls
}

pre_window <- function(ignition, days = 120) c(ignition - days, ignition - 1)
immediate_window <- function(end, days = 120) c(end + 1, end + days)
summer_window <- function(year) as.Date(c(sprintf("%d-12-01", year - 1), sprintf("%d-03-31", year)))

first_summer_after <- function(end) {
  year <- as.integer(format(end, "%Y")) + 1L
  while (summer_window(year)[1] <= end) year <- year + 1L
  year
}

latest_summer <- function(today) {
  y <- as.integer(format(today, "%Y"))
  if (today >= as.Date(sprintf("%d-04-01", y))) y else y - 1L
}

pre_summer_year <- function(ignition) {
  year <- as.integer(format(ignition, "%Y"))
  while (summer_window(year)[2] >= ignition) year <- year - 1L
  year
}

slugify <- function(text) gsub("(^-|-$)", "", gsub("[^a-z0-9]+", "-", tolower(text)))

# A north-up grid in EPSG:3577
make_grid <- function(x0, y0, res, width, height) {
  list(x0 = x0, y0 = y0, res = res, width = as.integer(width), height = as.integer(height))
}
grid_bounds <- function(g) c(g$x0, g$y0 - g$height * g$res, g$x0 + g$width * g$res, g$y0)
grid_gt <- function(g) c(g$x0, g$res, 0, g$y0, 0, -g$res)
grid_coarsen <- function(g, f) make_grid(g$x0, g$y0, g$res * f, g$width %/% f, g$height %/% f)
grid_pixel_ha <- function(g) g$res * g$res / 10000
grid_bounds_lonlat <- function(g) transform_bounds(grid_bounds(g), CRS, "EPSG:4326", densify_pts = 21L)

grid_for_bounds <- function(bounds, buffer_m = BUFFER_M, res = RES) {
  step <- res * COARSE
  x0 <- floor((bounds[1] - buffer_m) / step) * step
  x1 <- ceiling((bounds[3] + buffer_m) / step) * step
  y0 <- floor((bounds[2] - buffer_m) / step) * step
  y1 <- ceiling((bounds[4] + buffer_m) / step) * step
  make_grid(x0, y1, res, round((x1 - x0) / res), round((y1 - y0) / res))
}

# ---------------------------------------------------------------------------------------------
# Web services: plain JSON
# ---------------------------------------------------------------------------------------------

get_json <- function(url, query = NULL, body = NULL, tries = 4) {
  if (!is.null(query)) {
    url <- paste0(url, "?", paste(names(query), vapply(query, function(v)
      curl::curl_escape(as.character(v)), ""), sep = "=", collapse = "&"))
  }
  for (i in seq_len(tries)) {
    res <- tryCatch({
      h <- curl::new_handle(timeout = 180)
      if (!is.null(body)) {
        curl::handle_setheaders(h, "Content-Type" = "application/json")
        curl::handle_setopt(h, postfields = toJSON(body, auto_unbox = TRUE, digits = NA))
      }
      r <- curl::curl_fetch_memory(url, handle = h)
      if (r$status_code >= 400) stop("HTTP ", r$status_code)
      txt <- rawToChar(r$content)
      d <- fromJSON(txt, simplifyVector = FALSE)
      if (is.list(d) && !is.null(d$error)) stop(url, ": ", toJSON(d$error, auto_unbox = TRUE))
      list(d = d, txt = txt)
    }, error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (i == tries) stop(res)
    Sys.sleep(2^i)
  }
}

ms_to_date <- function(ms) {
  if (is.null(ms)) return(NULL)
  as.Date(as.POSIXct(ms / 1000, origin = "1970-01-01", tz = "UTC"))
}

fire_record <- function(a) {
  list(fire_id = a$FIRE_ID, name = trimws(if (is.null(a$FIRE_NAME)) "" else a$FIRE_NAME),
       area_ha = a$FIRE_AR_HA, ignition = ms_to_date(a$IGN_DATE), end = ms_to_date(a$INCD_END),
       cause = if (is.null(a$IGN_CAUSE1)) "" else a$IGN_CAUSE1)
}

# Read every geometry (as WKB) from a GeoJSON file, with GDAL
read_geoms <- function(path) {
  v <- new(GDALVector, path)
  on.exit(v$close())
  f <- v$fetch(-1)
  Filter(function(g) length(g) > 0, f[[attr(f, "gis")$geom_col_name %||% ncol(f)]])
}
`%||%` <- function(a, b) if (is.null(a)) b else a

get_fire <- function(fire_id) {
  path <- file.path(DOWNLOADS, "fires", sprintf("%d.geojson", as.integer(fire_id)))
  if (!file.exists(path)) {
    r <- get_json(paste0(FIRE_LAYER, "/query"), query = list(
      where = sprintf("FIRE_ID=%d", as.integer(fire_id)), outFields = FIRE_FIELDS,
      returnGeometry = "true", outSR = 4326, maxAllowableOffset = 0.0001, geometryPrecision = 6,
      f = "geojson"))
    if (!length(r$d$features)) stop("Fire ", fire_id, " not found in Fire History")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(r$txt, path)
  }
  d <- fromJSON(path, simplifyVector = FALSE)
  rec <- fire_record(d$features[[1]]$properties)
  parts <- lapply(read_geoms(path), function(g) g_buffer(g_make_valid(g), 0, quad_segs = 16L))
  list(rec = rec, geom = Reduce(g_union, parts))
}

fires_in_area <- function(bounds_ll, start, end) {
  r <- get_json(paste0(FIRE_LAYER, "/query"), query = list(
    where = sprintf("IGN_DATE >= DATE '%s' AND IGN_DATE <= DATE '%s'", start, end),
    geometry = paste(sprintf("%.5f", bounds_ll), collapse = ","),
    geometryType = "esriGeometryEnvelope", inSR = 4326, spatialRel = "esriSpatialRelIntersects",
    outFields = "FIRE_ID,FIRE_NAME", returnGeometry = "true", outSR = 4326,
    maxAllowableOffset = 0.0003, f = "geojson"))
  if (!length(r$d$features)) return(list())
  tmp <- tempfile(fileext = ".geojson")
  on.exit(unlink(tmp))
  writeLines(r$txt, tmp)
  read_geoms(tmp)
}

stac_search <- function(bbox_ll, start, end) {
  body <- list(collections = DEA_COLLECTIONS, bbox = bbox_ll, limit = 250,
               datetime = sprintf("%sT00:00:00Z/%sT23:59:59Z", start, end))
  items <- list()
  url <- DEA_STAC_SEARCH
  method <- "POST"
  for (i in 1:100) {
    d <- if (method == "POST") get_json(url, body = body)$d else get_json(url)$d
    items <- c(items, d$features)
    nxt <- Filter(function(l) identical(l$rel, "next"), d$links)
    if (!length(nxt)) break
    link <- nxt[[1]]
    url <- link$href
    method <- toupper(link$method %||% "GET")
    if (method == "POST") {
      body <- if (isTRUE(link$merge)) modifyList(body, link$body %||% list()) else (link$body %||% body)
    }
  }
  items
}

href <- function(item, asset) {
  sub("s3://dea-public-data/", "https://dea-public-data.s3.ap-southeast-2.amazonaws.com/",
      item$assets[[asset]]$href, fixed = TRUE)
}

# [(date, items, mean tile cloud %)], least cloudy first
group_by_date <- function(items) {
  d <- vapply(items, function(it) substr(it$properties$datetime, 1, 10), "")
  out <- lapply(unique(d), function(day) {
    its <- items[d == day]
    cc <- unlist(lapply(its, function(it) it$properties[["eo:cloud_cover"]]))
    list(date = day, items = its, cloud = if (length(cc)) mean(cc) else 100)
  })
  if (!length(out)) return(out)
  out[order(vapply(out, `[[`, 0, "cloud"), vapply(out, `[[`, "", "date"))]
}

# ---------------------------------------------------------------------------------------------
# Raster I/O with GDAL (gdalraster)
# ---------------------------------------------------------------------------------------------

# Band 1 of a file or URL warped onto the grid (nearest neighbour, full resolution, no overviews):
# a height x width matrix
warp_to_grid <- function(url, g) {
  src <- if (grepl("^https?://", url)) paste0("/vsicurl/", url) else url
  dst <- sprintf("/vsimem/warp_%d_%s.tif", Sys.getpid(), basename(tempfile()))
  on.exit(vsi_unlink(dst), add = TRUE)
  for (i in 1:4) {
    res <- tryCatch({
      warp(src, dst, t_srs = CRS, quiet = TRUE,
           cl_arg = c("-te", grid_bounds(g), "-ts", g$width, g$height, "-r", "near",
                      "-ovr", "NONE", "-of", "GTiff"))
      ds <- new(GDALRaster, dst)
      v <- as.numeric(read_ds(ds, bands = 1))
      nd <- ds$getNoDataValue(1)
      ds$close()
      v[is.na(v)] <- if (is.na(nd)) 0 else nd   # read_ds gives NA for nodata: put the raw value back
      matrix(v, nrow = g$height, ncol = g$width, byrow = TRUE)
    }, error = function(e) e)
    if (!inherits(res, "error")) return(res)
    missing <- grepl("403|404", conditionMessage(res))
    if (i == 4 || (missing && i == 3)) stop(res)
    Sys.sleep(3 * i)
  }
}

# Mosaic one date's tiles: clear (fmask) reflectance per band, and the clear fraction
read_date <- function(items, g) {
  out <- lapply(BANDS, function(b) matrix(NaN, g$height, g$width))
  filled <- matrix(FALSE, g$height, g$width)
  for (it in items) {
    r <- tryCatch({
      fm <- warp_to_grid(href(it, FMASK), g)
      ok <- (fm == FMASK_CLEAR) & !filled
      if (!any(ok)) NULL else
        list(ok = ok, refl = lapply(BANDS, function(a) to_reflectance(warp_to_grid(href(it, a), g))))
    }, error = function(e) {
      message(sprintf("    skipped tile %s (%s): %s", it$id, substr(it$properties$datetime, 1, 10),
                      substr(conditionMessage(e), 1, 80)))
      NULL
    })
    if (is.null(r)) next
    ok <- r$ok
    for (x in r$refl) ok <- ok & is.finite(x)
    for (b in names(BANDS)) out[[b]][ok] <- r$refl[[b]][ok]
    filled <- filled | ok
  }
  list(arrs = out, frac = mean(filled))
}

write_composite <- function(path, comp, g, meta) {
  nms <- c(names(BANDS), "n_clear")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ds <- create("GTiff", path, g$width, g$height, length(nms), "Float32",
               options = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"), return_obj = TRUE)
  ds$setGeoTransform(grid_gt(g))
  ds$setProjection(srs_to_wkt(CRS))
  for (i in seq_along(nms)) {
    ds$setNoDataValue(i, NaN)
    ds$setDescription(i, nms[i])
    ds$write(i, 0, 0, g$width, g$height, as.vector(t(comp[[nms[i]]])))
  }
  ds$close()
  writeLines(toJSON(meta, auto_unbox = TRUE, pretty = TRUE, digits = NA),
             sub("\\.tif$", ".json", path))
}

read_composite <- function(path) {
  ds <- new(GDALRaster, path)
  on.exit(ds$close())
  w <- ds$getRasterXSize()
  h <- ds$getRasterYSize()
  comp <- list()
  for (i in seq_len(ds$getRasterCount())) {
    comp[[ds$getDescription(i)]] <- matrix(as.numeric(read_ds(ds, bands = i)), h, w, byrow = TRUE)
  }
  list(comp = comp, meta = fromJSON(sub("\\.tif$", ".json", path), simplifyVector = TRUE))
}

row_nanmedian <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) return(matrixStats::rowMedians(m, na.rm = TRUE))
  apply(m, 1, function(v) stats::median(v, na.rm = TRUE))
}

# Median of up to `max_dates` clear-ish dates in [start, end]; cached as a GeoTIFF
composite <- function(site_dir, label, bbox_ll, start, end, g, max_dates = 4, min_clear = 0.25,
                      refresh = FALSE, min_cover = 0.9) {
  path <- file.path(site_dir, sprintf("%s_%dm.tif", label, as.integer(g$res)))
  if (file.exists(path) && !refresh) {
    rc <- read_composite(path)
    if ((rc$meta$valid_fraction %||% 0) >= min_cover - 0.05 || isTRUE(rc$meta$searched_all)) {
      cat(sprintf("  %s: cached, dates %s\n", label, paste(rc$meta$dates, collapse = ", ")))
      return(rc)
    }
    cat(sprintf("  %s: cached composite covers only %s: rebuilding\n", label, rc$meta$valid_fraction))
  }
  groups <- group_by_date(stac_search(bbox_ll, start, end))
  candidates <- Filter(function(x) x$cloud < 80, groups)
  cat(sprintf("  %s: %s to %s: %d dates found, %d under 80%% tile cloud; reading up to %d\n",
              label, start, end, length(groups), length(candidates), max_dates))
  stacks <- lapply(BANDS, function(b) list())
  used <- character()
  fracs <- numeric()
  pos <- 0
  cover <- matrix(FALSE, g$height, g$width)
  limit <- min(length(candidates), 3 * max_dates)
  while (pos < limit && (length(used) < max_dates || mean(cover) < min_cover)) {
    batch <- candidates[(pos + 1):min(pos + 3, length(candidates))]
    pos <- pos + length(batch)
    reads <- parallel::mclapply(batch, function(x) read_date(x$items, g), mc.cores = 3)
    for (k in seq_along(batch)) {
      r <- reads[[k]]
      if (inherits(r, "try-error")) stop(r)
      cat(sprintf("    %s: tile cloud %4.0f%%, clear over study area %5.1f%%\n",
                  batch[[k]]$date, batch[[k]]$cloud, 100 * r$frac))
      if (r$frac >= min_clear) {
        used <- c(used, batch[[k]]$date)
        fracs <- c(fracs, round(r$frac, 3))
        cover <- cover | is.finite(r$arrs$nir)
        for (b in names(BANDS)) stacks[[b]][[length(stacks[[b]]) + 1]] <- r$arrs[[b]]
      }
    }
  }
  cat(sprintf("    %d dates used; study area covered %.1f%%\n", length(used), 100 * mean(cover)))
  o <- order(used)
  meta <- list(label = label, window = c(format(start), format(end)), dates = used[o],
               clear_fraction = as.list(setNames(fracs, used)), res_m = g$res,
               collections = DEA_COLLECTIONS, n_dates_found = length(groups), searched_all = TRUE)
  if (!length(used)) {
    cat(sprintf("  %s: no usable dates\n", label))
    return(list(comp = NULL, meta = meta))
  }
  comp <- lapply(stacks, function(s) {
    m <- row_nanmedian(do.call(cbind, lapply(s, as.vector)))
    m[is.na(m)] <- NaN
    matrix(m, g$height, g$width)
  })
  comp$n_clear <- Reduce(`+`, lapply(stacks$nir, function(x) is.finite(x) * 1))
  meta$valid_fraction <- round(mean(is.finite(comp$nir)), 3)
  write_composite(path, comp, g, meta)
  read_composite(path)   # analyse the stored (Float32) values, as Python does
}

# ---------------------------------------------------------------------------------------------
# Geometry with GDAL (gdalraster's g_* functions, WKB in and out)
# ---------------------------------------------------------------------------------------------

to_albers <- function(wkb) g_transform(wkb, "EPSG:4326", CRS)

geom_bounds <- function(wkb) {
  e <- g_envelope(wkb)                     # xmin, xmax, ymin, ymax
  unname(c(e[1], e[3], e[2], e[4]))        # xmin, ymin, xmax, ymax
}

# TRUE where a pixel centre is inside any of the (EPSG:3577) geometries
rasterize_mask <- function(geoms, g) {
  geoms <- Filter(function(x) length(x) > 0 && !g_is_empty(x), geoms)
  if (!length(geoms)) return(matrix(FALSE, g$height, g$width))
  csv <- tempfile(fileext = ".csv")
  dst <- sprintf("/vsimem/mask_%d_%s.tif", Sys.getpid(), basename(tempfile()))
  on.exit({ unlink(csv); vsi_unlink(dst) }, add = TRUE)
  writeLines(c("id,WKT", sprintf('%d,"%s"', seq_along(geoms),
                                 vapply(geoms, function(x) g_wk2wk(x), ""))), csv)
  rasterize(csv, dst, burn_value = 1, te = grid_bounds(g), ts = c(g$width, g$height),
            dtName = "Byte", init = 0, fmt = "GTiff", add_options = c("-a_srs", CRS), quiet = TRUE)
  ds <- new(GDALRaster, dst)
  v <- read_ds(ds, bands = 1)
  ds$close()
  matrix(!is.na(v) & v == 1, g$height, g$width, byrow = TRUE)   # background 0 may read as NA
}

# ---------------------------------------------------------------------------------------------
# Analysis of one fire
# ---------------------------------------------------------------------------------------------

nanmedian <- function(x) if (any(!is.na(x))) stats::median(x, na.rm = TRUE) else NaN
pctl <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE, type = 7))   # = numpy default

analyse_fire <- function(key, fire_id, region = "", max_dates = 4, refresh = FALSE, today = Sys.Date()) {
  f <- get_fire(fire_id)
  rec <- f$rec
  geom_ll <- f$geom
  if (is.null(rec$ignition)) stop("Fire ", fire_id, " has no ignition date")
  end <- rec$end %||% (rec$ignition + 30)
  cat(sprintf("\n== %s: %s (FIRE_ID %d), ignited %s, out %s, official area %s ha ==\n", key,
              rec$name, as.integer(fire_id), rec$ignition, end,
              format(round(rec$area_ha), big.mark = ",")))
  site_dir <- file.path(DOWNLOADS, key)
  out_dir <- file.path(OUTPUTS, key)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  geom <- to_albers(geom_ll)
  grid <- grid_for_bounds(geom_bounds(geom))
  if (max(grid$width, grid$height) > 4500) grid <- grid_for_bounds(geom_bounds(geom), res = RES * 2)
  cgrid <- grid_coarsen(grid, COARSE)
  bbox_ll <- grid_bounds_lonlat(grid)
  cat(sprintf("  grid %d x %d px at %.0f m (%.1f x %.1f km), EPSG:3577\n", grid$width, grid$height,
              grid$res, grid$width * grid$res / 1000, grid$height * grid$res / 1000))

  # --- composites ------------------------------------------------------------------------------
  w <- pre_window(rec$ignition)
  pre <- composite(site_dir, "pre", bbox_ll, w[1], w[2], grid, max_dates, refresh = refresh)
  if (is.null(pre$comp) || (pre$meta$valid_fraction %||% 0) < 0.6) {
    w <- pre_window(rec$ignition, days = 365)
    pre <- composite(site_dir, "pre365", bbox_ll, w[1], w[2], grid, max_dates, refresh = refresh)
  }
  py <- pre_summer_year(rec$ignition)
  base <- list(comp = NULL, meta = list())
  if (summer_window(py)[1] >= as.Date("2015-07-01")) {
    w <- summer_window(py)
    base <- composite(site_dir, sprintf("presummer%d", py), bbox_ll, w[1], w[2], grid, max_dates,
                      refresh = refresh)
  }
  if (is.null(base$comp) || (base$meta$valid_fraction %||% 0) < 0.8) {
    base <- pre
    base_name <- sprintf("just before the fire (%s; no earlier summer available)",
                         substr(pre$meta$dates[1], 1, 4))
  } else {
    base_name <- sprintf("summer %d-%s", py - 1, substr(py, 3, 4))
  }
  w <- immediate_window(end)
  imm <- composite(site_dir, "immediate", bbox_ll, w[1], w[2], grid, max_dates, refresh = refresh)
  first <- first_summer_after(end)
  last <- latest_summer(today)
  initial <- first > last
  if (initial) {
    if (is.null(imm$comp)) stop("No clear imagery yet from after the fire")
    s1 <- imm
    nxt <- end + 1
    first <- as.integer(format(nxt, "%Y")) + (format(nxt, "%m") == "12")
    cat("  no full summer since the fire was out: severity from the images just after it\n")
  } else {
    w <- summer_window(first)
    s1 <- composite(site_dir, sprintf("summer%d", first), bbox_ll, w[1], w[2], grid, max_dates,
                    refresh = refresh)
  }
  if (is.null(pre$comp) || is.null(s1$comp)) {
    stop("Not enough clear imagery before the fire or in the first summer after it")
  }
  summers <- list()
  summers[[as.character(first)]] <- list(comp = s1$comp, meta = s1$meta, grid = grid)
  if (last >= first + 1) for (y in (first + 1):last) {
    w <- summer_window(y)
    cm <- composite(site_dir, sprintf("summer%d", y), bbox_ll, w[1], w[2], cgrid,
                    max(2, max_dates - 1), refresh = refresh)
    if (!is.null(cm$comp)) summers[[as.character(y)]] <- list(comp = cm$comp, meta = cm$meta, grid = cgrid)
  }

  # --- burn severity ---------------------------------------------------------------------------
  inside <- rasterize_mask(list(geom), grid)
  others_ll <- fires_in_area(bbox_ll, pre_window(rec$ignition, 365)[1], summer_window(last)[2])
  others <- lapply(others_ll, function(x) g_buffer(to_albers(x), 200, quad_segs = 16L))
  burnt_any <- rasterize_mask(c(others, list(g_buffer(geom, 500, quad_segs = 16L))), grid)
  control <- !burnt_any

  nbr_pre <- nbr(base$comp)
  nbr_s1 <- nbr(s1$comp)
  dnbr_raw <- nbr_pre - nbr_s1
  ok_ctrl <- control & is.finite(dnbr_raw)
  offset <- if (sum(ok_ctrl) > 1000) nanmedian(dnbr_raw[ok_ctrl]) else 0
  dnbr <- dnbr_raw - offset
  cls <- classify_dnbr(dnbr)
  px <- grid_pixel_ha(grid)
  valid_in <- inside & is.finite(dnbr)
  burnt <- is.finite(dnbr) & dnbr >= SEVERITY$lo[2]
  stats <- list(
    site = key, fire_id = as.integer(fire_id), fire_name = rec$name, region = region,
    ignition = format(rec$ignition), declared_out = format(end), cause = rec$cause,
    official_area_ha = round(rec$area_ha, 1),
    boundary_area_ha = round(sum(inside) * px, 1),
    valid_inside_pct = round(100 * sum(valid_in) / max(sum(inside), 1), 1),
    dnbr_offset = round(offset, 4),
    burnt_inside_ha = round(sum(burnt & valid_in) * px, 1),
    burnt_inside_pct = round(100 * sum(burnt & valid_in) / max(sum(valid_in), 1), 1),
    burnt_outside_ha = round(sum(burnt & control) * px, 1),
    control_ha = round(sum(ok_ctrl) * px, 1),
    pre_dates = paste(pre$meta$dates, collapse = ";"), baseline = base_name,
    baseline_dates = paste(base$meta$dates, collapse = ";"),
    immediate_dates = paste(imm$meta$dates, collapse = ";"),
    severity_basis = if (initial) "initial (just after the fire)" else sprintf("summer %d", first),
    severity_summer = first, severity_dates = paste(s1$meta$dates, collapse = ";"),
    res_m = grid$res)
  for (i in seq_len(nrow(SEVERITY))) {
    n <- sum(cls == i - 1L & valid_in)
    stats[[paste0("ha_", slugify(SEVERITY$name[i]))]] <- round(n * px, 1)
    stats[[paste0("pct_", slugify(SEVERITY$name[i]))]] <- round(100 * n / max(sum(valid_in), 1), 1)
  }
  cat("  reflectance ranges (pre): ", paste(vapply(names(BANDS), function(b)
    sprintf("%s %.3f-%.3f", b, pctl(pre$comp[[b]], 0.01), pctl(pre$comp[[b]], 0.99)), ""),
    collapse = ", "), "\n", sep = "")
  cat(sprintf("  severity baseline: %s\n", base_name))
  cat(sprintf("  NBR pre %.2f to %.2f; dNBR (offset %+.3f removed) %.2f to %.2f\n",
              pctl(nbr_pre, 0.01), pctl(nbr_pre, 0.99), offset, pctl(dnbr, 0.01), pctl(dnbr, 0.99)))
  cat(sprintf(paste0("  boundary %s ha (official %s); valid %s%%; burnt (dNBR >= 0.10) %s ha = %s%%",
                     " of the boundary; dNBR >= 0.10 outside any fire %s ha of %s ha unburnt land\n"),
              format(round(stats$boundary_area_ha), big.mark = ","),
              format(round(stats$official_area_ha), big.mark = ","), stats$valid_inside_pct,
              format(round(stats$burnt_inside_ha), big.mark = ","), stats$burnt_inside_pct,
              format(round(stats$burnt_outside_ha), big.mark = ","),
              format(round(stats$control_ha), big.mark = ",")))

  # --- recovery --------------------------------------------------------------------------------
  sub_grid <- function(m, f, g) m[seq(1, by = f, length.out = g$height), seq(1, by = f, length.out = g$width)]
  rows <- list()
  for (y in sort(as.integer(names(summers)))) {
    s <- summers[[as.character(y)]]
    g <- s$grid
    f <- if (g$res == grid$res) 1L else COARSE
    cls_g <- sub_grid(cls, f, g)
    in_g <- sub_grid(inside, f, g)
    ctrl_g <- sub_grid(control, f, g)
    d <- sub_grid(nbr_pre, f, g) - nbr(s$comp)
    nd <- ndvi(s$comp)
    ctrl_drop <- if (sum(is.finite(d[ctrl_g])) > 100) nanmedian(d[ctrl_g]) else 0
    row <- list(site = key, summer = y, years_since = y - as.integer(format(rec$ignition, "%Y")),
                dates = paste(s$meta$dates, collapse = ";"), control_drop = round(ctrl_drop, 4),
                control_ndvi = round(nanmedian(nd[ctrl_g]), 4))
    for (i in seq_len(nrow(SEVERITY))) {
      k <- slugify(SEVERITY$name[i])
      sel <- in_g & cls_g == i - 1L & is.finite(d)
      if (sum(sel) < 50) {
        row[[paste0("deficit_", k)]] <- NA
        row[[paste0("ndvi_", k)]] <- NA
      } else {
        row[[paste0("deficit_", k)]] <- round(nanmedian(d[sel]) - ctrl_drop, 4)
        row[[paste0("ndvi_", k)]] <- round(nanmedian(nd[sel]), 4)
      }
    }
    rows[[length(rows) + 1]] <- row
  }
  for (k in c("high", "moderate-high")) {
    a <- rows[[1]][[paste0("deficit_", k)]]
    b <- rows[[length(rows)]][[paste0("deficit_", k)]]
    stats[[sprintf("recovered_%s_pct", k)]] <-
      if (!is.na(a) && !is.na(b) && a > 0) round(100 * (1 - b / a), 1) else NA
  }
  stats$latest_summer <- max(as.integer(names(summers)))

  write_rows(file.path(out_dir, "summary.csv"), list(stats))
  write_rows(file.path(out_dir, "recovery.csv"), rows)
  cat(sprintf("  wrote %s: summary.csv, recovery.csv\n", sub(paste0(ROOT, "/"), "", out_dir, fixed = TRUE)))
  invisible(stats)
}

write_rows <- function(path, rows) {
  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, check.names = FALSE)))
  utils::write.csv(df, path, row.names = FALSE, na = "")
}

# ---------------------------------------------------------------------------------------------
# Compare with tas_fires.py
# ---------------------------------------------------------------------------------------------

compare_with_python <- function(key) {
  cat(sprintf("\n== %s: R (gdalraster) vs Python (osgeo.gdal) ==\n", key))
  for (f in c("summary.csv", "recovery.csv")) {
    pp <- file.path(ROOT, "outputs", "tasmania", key, f)
    rp <- file.path(OUTPUTS, key, f)
    if (!file.exists(pp) || !file.exists(rp)) { cat("  ", f, ": missing\n"); next }
    p <- utils::read.csv(pp, check.names = FALSE, colClasses = "character")
    r <- utils::read.csv(rp, check.names = FALSE, colClasses = "character")
    cols <- union(names(p), names(r))
    cat(sprintf("  %s: %d columns x %d rows\n", f, length(cols), nrow(p)))
    n_same <- 0
    for (cl in cols) {
      if (!cl %in% names(p) || !cl %in% names(r)) { cat(sprintf("    %-28s only in %s\n", cl,
        if (cl %in% names(p)) "Python" else "R")); next }
      for (i in seq_len(max(nrow(p), nrow(r)))) {
        a <- p[[cl]][i]; b <- r[[cl]][i]
        na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
        same <- if (!is.na(na) && !is.na(nb)) isTRUE(all.equal(na, nb)) else identical(a, b)
        if (same) n_same <- n_same + 1 else
          cat(sprintf("    %-28s row %d: Python %s | R %s\n", cl, i, a, b))
      }
    }
    cat(sprintf("    %d values identical\n", n_same))
  }
  pdir <- file.path(ROOT, "downloads", "tasmania", key)
  for (rf in sort(list.files(file.path(DOWNLOADS, key), "\\.tif$", full.names = TRUE))) {
    pf <- file.path(pdir, basename(rf))
    if (!file.exists(pf)) { cat("  ", basename(rf), ": no Python composite\n"); next }
    a <- read_composite(pf)$comp
    b <- read_composite(rf)$comp
    msg <- vapply(names(a), function(n) {
      x <- a[[n]]; y <- b[[n]]
      nanpat <- sum(is.na(x) != is.na(y))
      ok <- !is.na(x) & !is.na(y)
      sprintf("%s max|diff| %.3g, NaN mismatch %d", n, if (any(ok)) max(abs(x[ok] - y[ok])) else 0, nanpat)
    }, "")
    cat(sprintf("  %-22s %s\n", basename(rf), paste(msg, collapse = "; ")))
  }
}

# ---------------------------------------------------------------------------------------------

main <- function(args = commandArgs(TRUE)) {
  sites <- args[which(args == "--site") + 1]
  if (!length(sites) || "--list" %in% args) {
    cat("Built-in sites:", paste(names(SITES), collapse = ", "), "\n")
    return(invisible(0))
  }
  for (key in sites) {
    s <- SITES[[key]]
    if (is.null(s)) stop("Unknown site: ", key)
    if (!"--compare-only" %in% args) {
      analyse_fire(key, s$fire_id, s$region, refresh = "--refresh" %in% args)
    }
    if (any(c("--compare", "--compare-only") %in% args)) compare_with_python(key)
  }
  invisible(0)
}

if (sys.nframe() == 0L) main()
