## make_bioregion_viewer.R
##
## Builds Outputs/bioregion_viewer.html: an interactive Leaflet viewer that
## lets you switch between the three GDM-filled iNaturalist bioregion layers
## (K=194/53/15) and their matching official National Ecological Framework
## boundaries (ecoregions/ecoprovinces/ecozones), using the same
## Hungarian-algorithm colour matching as bioregionalization_recolor_match.R
## so corresponding regions share colours across layers.
##
## Reuses the cell-level assignments already saved by
## bioregionalization_inat2025_gdm_fill.R and _multiK.R (no recompute of the
## community matrix / NPP extraction / GDM steps).
##
## Output: Outputs/bioregion_viewer.html

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(data.table)
  library(jsonlite); library(clue)
})

ECO_FILE <- "Data/nef_ca_ter_ecoregion_v2_2.geojson"
ECOZONE_NAMES_FILE <- "Data/EckertJBiog/Code and Data/Data/Ecoframework Data/ecozone.data.NEW.RDS"
OUT_DIR  <- "Outputs"
sf_use_s2(FALSE)

pal_fn <- colorRampPalette(c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                              "#80BA5A","#E68310","#008695","#CF1C90","#f97b72",
                              "#4b4b8f","#A5AA99"))

match_official_colors <- function(cluster_filled, official_id, pal) {
  cluster_filled <- as.character(cluster_filled)
  official_id    <- as.character(official_id)
  ok  <- !is.na(official_id) & !is.na(cluster_filled)
  tab <- table(cluster_filled[ok], official_id[ok])
  clusters  <- rownames(tab)
  officials <- colnames(tab)
  n <- max(length(clusters), length(officials))
  cost <- matrix(0, n, n)
  cost[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  assign <- clue::solve_LSAP(cost, maximum = TRUE)
  cluster_pal <- setNames(pal, sort(as.integer(unique(cluster_filled))))
  official_color <- setNames(rep(NA_character_, length(officials)), officials)
  for (i in seq_len(length(clusters))) {
    j <- assign[i]
    if (j <= length(officials)) official_color[officials[j]] <- cluster_pal[[clusters[i]]]
  }
  official_color[is.na(official_color)] <- "#BBBBBB"
  official_color
}

geo_str <- function(sf_obj) {
  tmp <- tempfile(fileext = ".geojson")
  st_write(sf_obj, tmp, quiet = TRUE, delete_dsn = TRUE)
  s <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  unlink(tmp)
  s
}

# ═══ Official boundaries: load + simplify once ════════════════════════════════
cat("Loading & simplifying official boundaries...\n")
eco_raw <- st_read(ECO_FILE, quiet = TRUE)
eco_raw <- st_make_valid(eco_raw)
eco_raw <- st_simplify(eco_raw, dTolerance = 0.01)
eco_raw <- st_make_valid(eco_raw)

eco_diss <- function(field) {
  eco_raw |> group_by(across(all_of(field))) |> summarise(.groups = "drop") |> st_make_valid()
}

ecoregion_off   <- eco_diss("ECOREGION_ID")
ecozone_off     <- eco_diss("ECOZONE_ID")
ecoprovince_off <- eco_diss("ECOPROVINCE_ID")

# ecoregion names already present per-feature in eco_raw; recover one per group
eco_names <- eco_raw |> st_drop_geometry() |>
  distinct(ECOREGION_ID, ECOREGION_NAME_EN) |> filter(!is.na(ECOREGION_ID))
ecoregion_off <- left_join(ecoregion_off, eco_names, by = "ECOREGION_ID")

# ecozone names: nearest-centroid match against the Eckert framework's named
# 15-ecozone table (same underlying NEF classification used elsewhere in
# this project for the 194 ecoregions)
eck_zones <- readRDS(ECOZONE_NAMES_FILE)
zone_cent <- st_coordinates(st_centroid(ecozone_off))
zone_name <- sapply(seq_len(nrow(ecozone_off)), function(i) {
  d <- sqrt((eck_zones$Longitude - zone_cent[i,1])^2 + (eck_zones$Latitude - zone_cent[i,2])^2)
  eck_zones$Ecozone[which.min(d)]
})
ecozone_off$ECOZONE_NAME <- zone_name

ecoprovince_off$ECOPROVINCE_NAME <- sprintf("Ecoprovince %.1f", ecoprovince_off$ECOPROVINCE_ID)

cat(sprintf("  ecoregions: %d | ecozones: %d | ecoprovinces: %d\n",
            nrow(ecoregion_off), nrow(ecozone_off), nrow(ecoprovince_off)))

# ═══ Per-target: colour match + assemble geojson ═══════════════════════════════
build_target <- function(cells_geojson, cells_csv, official_sf, official_field,
                          official_name_field, K, label) {
  cat(sprintf("── %s (K=%d) ──\n", label, K))
  cells <- st_read(cells_geojson, quiet = TRUE)
  csv <- fread(cells_csv)
  if (!"official_id" %in% names(csv)) setnames(csv, official_field, "official_id")
  csv[, official_id := as.character(official_id)]

  pal <- pal_fn(K)
  cluster_ids_sorted <- sort(as.integer(unique(csv$cluster_filled)))
  cluster_color <- setNames(pal, cluster_ids_sorted)
  official_color <- match_official_colors(csv$cluster_filled, csv$official_id, pal)

  cells$color <- cluster_color[as.character(cells$cluster_filled)]
  cells_min <- cells[, c("cell_id","cluster_filled","data_status","npp_mean","color")]

  off <- official_sf
  off_id_chr <- as.character(st_drop_geometry(off)[[official_field]])
  off$color <- official_color[off_id_chr]
  off$color[is.na(off$color)] <- "#BBBBBB"
  off$region_name <- st_drop_geometry(off)[[official_name_field]]
  off_min <- off[, c("region_name","color")]
  names(off_min)[1] <- "region_name"

  list(
    filled_geojson   = geo_str(cells_min),
    official_geojson = geo_str(off_min),
    K = K,
    n_sampled = sum(csv$data_status == "sampled"),
    n_filled  = nrow(csv)
  )
}

ecoregion_res <- build_target(
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.geojson"),
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.csv"),
  ecoregion_off, "ECOREGION_ID", "ECOREGION_NAME_EN", 194, "ecoregion")

ecozone_res <- build_target(
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecozone_k15.geojson"),
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecozone_k15.csv"),
  ecozone_off, "ECOZONE_ID", "ECOZONE_NAME", 15, "ecozone")

ecoprovince_res <- build_target(
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecoprovince_k53.geojson"),
  file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecoprovince_k53.csv"),
  ecoprovince_off, "ECOPROVINCE_ID", "ECOPROVINCE_NAME", 53, "ecoprovince")

# ═══ Assemble HTML ══════════════════════════════════════════════════════════════
cat("Assembling HTML...\n")
meta <- list(
  ecoregion   = list(K = ecoregion_res$K,   n_sampled = ecoregion_res$n_sampled,   n_filled = ecoregion_res$n_filled,   label = "Ecoregions"),
  ecozone     = list(K = ecozone_res$K,     n_sampled = ecozone_res$n_sampled,     n_filled = ecozone_res$n_filled,     label = "Ecozones"),
  ecoprovince = list(K = ecoprovince_res$K, n_sampled = ecoprovince_res$n_sampled, n_filled = ecoprovince_res$n_filled, label = "Ecoprovinces")
)
meta_json <- toJSON(meta, auto_unbox = TRUE)

tmpl <- paste(readLines("bioregion_viewer_template.html", warn = FALSE), collapse = "\n")
html <- tmpl
html <- sub("__META_JSON__",                 meta_json,                       html, fixed = TRUE)
html <- sub("__GEO_ECOREGION_FILLED__",      ecoregion_res$filled_geojson,    html, fixed = TRUE)
html <- sub("__GEO_ECOREGION_OFFICIAL__",    ecoregion_res$official_geojson,  html, fixed = TRUE)
html <- sub("__GEO_ECOZONE_FILLED__",        ecozone_res$filled_geojson,      html, fixed = TRUE)
html <- sub("__GEO_ECOZONE_OFFICIAL__",      ecozone_res$official_geojson,    html, fixed = TRUE)
html <- sub("__GEO_ECOPROVINCE_FILLED__",    ecoprovince_res$filled_geojson,  html, fixed = TRUE)
html <- sub("__GEO_ECOPROVINCE_OFFICIAL__",  ecoprovince_res$official_geojson, html, fixed = TRUE)

out_path <- file.path(OUT_DIR, "bioregion_viewer.html")
writeLines(html, out_path)
cat(sprintf("\nSaved: %s (%.1f MB)\n", out_path, file.size(out_path) / 1024^2))
