## bioregionalization_inat2025_gdm_fill_multiK.R
##
## Same GDM gap-filling approach as bioregionalization_inat2025_gdm_fill.R,
## but targeting the two COARSER levels of the National Ecological Framework
## hierarchy that are already present in nef_ca_ter_ecoregion_v2_2.geojson:
##   - ECOZONE_ID:     15 terrestrial ecozones of Canada
##   - ECOPROVINCE_ID: 53 ecoprovinces
## (194 ECOREGION_ID was handled in bioregionalization_inat2025_gdm_fill.R.)
##
## The expensive steps (community matrix, NPP zonal extraction, GDM fit +
## transform) do not depend on K, so they are run ONCE; only the cluster cut
## (cutree), nearest-centroid extrapolation, and mapping are repeated per
## target level.
##
## Outputs (Outputs/):
##   bioregions_inat2025_gdm_filled_cells_ecozone_k15.csv / .geojson
##   bioregions_inat2025_gdm_filled_cells_ecoprovince_k53.csv / .geojson
##   bioregions_inat2025_gdm_filled_map_ecozone_k15.png
##   bioregions_inat2025_gdm_filled_map_ecoprovince_k53.png

suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(data.table); library(sf)
  library(vegan); library(gdm); library(terra)
  library(ggplot2); library(patchwork)
})

set.seed(1)

IN_PARQUET <- "Data/iNaturalist/inaturalist-canada-dec2025_smaller.parquet"
ECO_FILE   <- "Data/nef_ca_ter_ecoregion_v2_2.geojson"
MODIS_DIR  <- "Data/MODIS8day"
OUT_DIR    <- "Outputs"

CELL_SIZE   <- 50000
MIN_RECORDS <- 20
MIN_TAXA    <- 8
ANIMAL_TAXA <- c("Aves","Mammalia","Insecta","Arachnida","Mollusca",
                 "Amphibia","Reptilia","Actinopterygii")
CANADA_AEA  <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

TARGETS <- list(
  list(field = "ECOZONE_ID",     label = "ecozone",     tag = "ecozone_k15"),
  list(field = "ECOPROVINCE_ID", label = "ecoprovince", tag = "ecoprovince_k53")
)

ari_fn <- function(cl, ec) {
  ok <- !is.na(ec) & !is.na(cl)
  tab <- table(cl[ok], ec[ok])
  n <- sum(tab); a <- rowSums(tab); b <- colSums(tab)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab)); sum_a <- sum(comb2(a)); sum_b <- sum(comb2(b))
  expected <- sum_a * sum_b / comb2(n); maxidx <- (sum_a + sum_b) / 2
  (sum_ij - expected) / (maxidx - expected)
}

# ═══ 1. Rebuild 2025 iNat community matrix (same as before) ═══════════════════
cat("── [1/6] Rebuilding 2025 iNat community matrix ─────────────\n")
eco <- st_read(ECO_FILE, quiet = TRUE)
TARGETS[[1]]$K <- length(unique(eco$ECOZONE_ID))
TARGETS[[2]]$K <- length(unique(eco$ECOPROVINCE_ID))
cat(sprintf("  Targets: %s\n", paste(sprintf("%s (K=%d)", sapply(TARGETS, `[[`, "label"),
                                              sapply(TARGETS, `[[`, "K")), collapse = ", ")))

ds <- open_dataset(IN_PARQUET)
d <- ds |>
  filter(year == "2025", quality_grade != "casual",
         (iconic_taxon_name == "Plantae") | (iconic_taxon_name %in% ANIMAL_TAXA)) |>
  select(longitude, latitude, scientific_name, iconic_taxon_name) |>
  collect()
setDT(d)
d <- d[!is.na(longitude) & !is.na(latitude) & !is.na(scientific_name) & scientific_name != ""]

d[, n_words := lengths(strsplit(scientific_name, "\\s+"))]
d[, w1 := sub("^(\\S+).*", "\\1", scientific_name)]
d[, w2 := ifelse(n_words >= 2, sub("^\\S+\\s+(\\S+).*", "\\1", scientific_name), NA_character_)]
BAD_W2 <- c("sp","sp.","spp","spp.","cf","cf.","aff","aff.","complex","group","x","hybrid")

plants <- d[iconic_taxon_name == "Plantae"]
plants[, taxon := paste0("Plant genus: ", w1)]
animals <- d[iconic_taxon_name %in% ANIMAL_TAXA & n_words >= 2 &
             grepl("^[a-z]+\\.?$", w2) & !(tolower(w2) %in% BAD_W2)]
animals[, taxon := paste0("Animal sp: ", w1, " ", w2)]

comb <- rbind(plants[, .(longitude, latitude, taxon)], animals[, .(longitude, latitude, taxon)])
cat(sprintf("  Combined: %s records, %s taxa\n", format(nrow(comb), big.mark=","), format(uniqueN(comb$taxon), big.mark=",")))

pts <- st_as_sf(comb, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
xy <- st_coordinates(st_transform(pts, CANADA_AEA))
comb[, `:=`(x = xy[,1], y = xy[,2])]
comb[, `:=`(col = floor(x / CELL_SIZE), row = floor(y / CELL_SIZE))]
comb[, cell_id := paste(col, row, sep = "_")]

cell_stats <- comb[, .(n_records = .N, n_taxa = uniqueN(taxon)), by = .(cell_id, col, row)]
well <- cell_stats[n_records >= MIN_RECORDS & n_taxa >= MIN_TAXA]
cat(sprintf("  %s well-sampled cells (of %s with any data)\n", format(nrow(well), big.mark=","), format(nrow(cell_stats), big.mark=",")))

comb_well <- comb[cell_id %in% well$cell_id]
mat_dt <- dcast(comb_well, cell_id ~ taxon, fun.aggregate = length, value.var = "taxon", fill = 0)
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- mat_dt$cell_id
mat <- mat[well$cell_id, ]
cat(sprintf("  Community matrix: %d cells x %d taxa\n", nrow(mat), ncol(mat)))

mat_hel <- decostand(mat, method = "hellinger")
hc <- hclust(dist(mat_hel, method = "euclidean"), method = "ward.D2")
cat("  Ward dendrogram built (K-independent; cutree happens per target below)\n")

# ═══ 2. Full 50km land grid over the ecoregion extent ═════════════════════════
cat("── [2/6] Building full land grid ────────────────────────────\n")
eco_p <- st_transform(eco, CANADA_AEA)
eco_union <- st_union(st_make_valid(eco_p))
bb <- st_bbox(eco_p)
col_range <- floor(bb["xmin"] / CELL_SIZE):floor(bb["xmax"] / CELL_SIZE)
row_range <- floor(bb["ymin"] / CELL_SIZE):floor(bb["ymax"] / CELL_SIZE)
full_grid <- CJ(col = col_range, row = row_range)
full_grid[, cell_id := paste(col, row, sep = "_")]

make_sq <- function(col, row) {
  xmin <- col * CELL_SIZE; xmax <- xmin + CELL_SIZE
  ymin <- row * CELL_SIZE; ymax <- ymin + CELL_SIZE
  st_polygon(list(matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin),
                          ncol = 2, byrow = TRUE)))
}
geoms <- mapply(make_sq, full_grid$col, full_grid$row, SIMPLIFY = FALSE)
full_sf <- st_sf(full_grid, geometry = st_sfc(geoms, crs = CANADA_AEA))
on_land <- lengths(st_intersects(full_sf, eco_union)) > 0
full_sf <- full_sf[on_land, ]
full_sf$cell_num <- seq_len(nrow(full_sf))
cat(sprintf("  %s land cells (50 km) across the full ecoregion extent\n", format(nrow(full_sf), big.mark=",")))

cent_xy <- st_coordinates(st_centroid(full_sf))
full_sf$x <- cent_xy[,1]; full_sf$y <- cent_xy[,2]

# ═══ 3. Environmental layer: 2025 NPP (PsnNet) mean & SD per cell ════════════
cat("── [3/6] Extracting 2025 NPP per grid cell ──────────────────\n")
npp_files <- list.files(MODIS_DIR, pattern = "PsnNet.*doy2025.*[.]tif$", full.names = TRUE)
cat(sprintf("  %d 2025 PsnNet (8-day NPP) files\n", length(npp_files)))

template <- rast(npp_files[1])
grid_v <- vect(st_transform(full_sf[, "cell_num"], crs(template)))
template_agg <- aggregate(template, fact = 4)
zone <- rasterize(grid_v, template_agg, field = "cell_num")

npp_cols <- vector("list", length(npp_files))
for (i in seq_along(npp_files)) {
  r  <- rast(npp_files[i])
  rc <- clamp(r, lower = -3, upper = 3, values = FALSE)
  ra <- aggregate(rc, fact = 4, fun = "mean", na.rm = TRUE)
  zm <- zonal(ra, zone, fun = "mean", na.rm = TRUE)
  names(zm) <- c("cell_num", "psn")
  npp_cols[[i]] <- as.data.table(zm)
  if (i %% 10 == 0) cat(sprintf("    %d/%d files\n", i, length(npp_files)))
}
npp_long <- rbindlist(npp_cols, idcol = "file_i")
npp_stats <- npp_long[, .(npp_mean = mean(psn, na.rm = TRUE),
                           npp_sd   = sd(psn, na.rm = TRUE)), by = cell_num]
full_sf <- merge(full_sf, npp_stats, by = "cell_num", all.x = TRUE)
cat(sprintf("  NPP extracted for %d / %d cells (%d NA)\n",
            sum(!is.na(full_sf$npp_mean)), nrow(full_sf), sum(is.na(full_sf$npp_mean))))

# ═══ 4. Fit GDM once (K-independent) ═══════════════════════════════════════════
cat("── [4/6] Fitting GDM ───────────────────────────────────────\n")
well_env <- as.data.table(st_drop_geometry(full_sf))[cell_id %in% well$cell_id]
well <- merge(well, well_env[, .(cell_id, cell_num, x, y, npp_mean, npp_sd)], by = "cell_id")
well <- well[!is.na(npp_mean) & !is.na(npp_sd)]
mat_fit <- mat[well$cell_id, ]
cat(sprintf("  %d well-sampled cells with valid NPP retained for GDM fitting\n", nrow(well)))

bray <- as.matrix(vegdist(mat_fit, method = "bray"))
bioDataTab <- data.frame(site = well$cell_num, bray)
predData_fit <- data.frame(site = well$cell_num, X = well$x, Y = well$y,
                            npp_mean = well$npp_mean, npp_sd = well$npp_sd)

sp_table <- formatsitepair(bioDataTab, bioFormat = 3, siteColumn = "site",
                            XColumn = "X", YColumn = "Y", predData = predData_fit)
mod <- gdm(sp_table, geo = TRUE)
cat(sprintf("  GDM deviance explained: %.1f%%\n", mod$explained))

full_valid <- full_sf[!is.na(full_sf$npp_mean) & !is.na(full_sf$npp_sd), ]
trans_all <- gdm.transform(mod, data = data.frame(
  X = full_valid$x, Y = full_valid$y,
  npp_mean = full_valid$npp_mean, npp_sd = full_valid$npp_sd))
trans_all <- as.data.table(trans_all)
trans_mat <- as.matrix(trans_all)
trans_well_idx <- match(well$cell_num, full_valid$cell_num)

# ── Attach official hierarchy IDs (ecozone/ecoprovince/ecoregion) per cell ────
full_valid_ll <- st_transform(full_valid, 4326)
cent_ll <- st_centroid(full_valid_ll)
eco_match <- suppressWarnings(st_join(cent_ll,
  eco[, c("ECOREGION_ID","ECOZONE_ID","ECOPROVINCE_ID")], join = st_intersects))
full_valid_ll$ECOREGION_ID     <- eco_match$ECOREGION_ID
full_valid_ll$ECOZONE_ID       <- eco_match$ECOZONE_ID
full_valid_ll$ECOPROVINCE_ID   <- eco_match$ECOPROVINCE_ID

cat(sprintf("  %s cells ready for extrapolation (%s originally well-sampled)\n",
            format(nrow(full_valid), big.mark=","), format(nrow(well), big.mark=",")))

# ═══ 5. Per-target: cutree, nearest-centroid extrapolation, map ═══════════════
pal <- colorRampPalette(c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                           "#80BA5A","#E68310","#008695","#CF1C90","#f97b72",
                           "#4b4b8f","#A5AA99"))

for (t in TARGETS) {
  cat(sprintf("── [5/6] Target: %s (K=%d) ──────────────────────────\n", t$label, t$K))
  clust_k <- cutree(hc, k = t$K)
  well_cluster <- clust_k[match(well$cell_id, names(clust_k))]

  trans_well <- trans_mat[trans_well_idx, , drop = FALSE]
  centroids <- rowsum(trans_well, group = well_cluster) / as.vector(table(well_cluster))
  nearest <- apply(trans_mat, 1, function(v) {
    dd <- sqrt(rowSums(sweep(centroids, 2, v)^2))
    as.integer(rownames(centroids)[which.min(dd)])
  })

  res <- full_valid_ll
  res$cluster_filled <- nearest
  res$data_status <- ifelse(res$cell_id %in% well$cell_id, "sampled", "gap-filled")
  res$official_id <- st_drop_geometry(res)[[t$field]]

  ari_filled  <- ari_fn(res$cluster_filled, res$official_id)
  ari_sampled <- ari_fn(res$cluster_filled[res$data_status == "sampled"],
                         res$official_id[res$data_status == "sampled"])
  ward_sampled_id <- res$official_id[match(well$cell_id, res$cell_id)]
  ari_ward <- ari_fn(well_cluster, ward_sampled_id)
  cat(sprintf("  ARI vs. official %ss — Ward (sampled): %.3f | GDM nearest-centroid (sampled): %.3f | GDM-filled (all): %.3f\n",
              t$label, ari_ward, ari_sampled, ari_filled))

  fwrite(as.data.table(st_drop_geometry(res))[, .(
    cell_id, cluster_filled, data_status, npp_mean, npp_sd, official_id
  )], file.path(OUT_DIR, sprintf("bioregions_inat2025_gdm_filled_cells_%s.csv", t$tag)))
  st_write(res[, c("cell_id","cluster_filled","data_status","npp_mean","npp_sd","official_id")],
            file.path(OUT_DIR, sprintf("bioregions_inat2025_gdm_filled_cells_%s.geojson", t$tag)),
            delete_dsn = TRUE, quiet = TRUE)

  sampled_only <- res[res$data_status == "sampled", ]
  p1 <- ggplot(sampled_only) +
    geom_sf(aes(fill = factor(cluster_filled)), color = NA) +
    scale_fill_manual(values = pal(t$K), guide = "none") +
    labs(title = sprintf("iNat 2025 bioregions (sparse, K=%d)", t$K),
         subtitle = sprintf("%s well-sampled cells", format(nrow(sampled_only), big.mark=","))) +
    theme_void(base_size = 10)

  p2 <- ggplot(res) +
    geom_sf(aes(fill = factor(cluster_filled)), color = NA) +
    scale_fill_manual(values = pal(t$K), guide = "none") +
    labs(title = sprintf("GDM-filled iNat 2025 bioregions (K=%d)", t$K),
         subtitle = sprintf("%s cells; NPP + geographic distance GDM (%.0f%% dev. explained)",
                             format(nrow(res), big.mark=","), mod$explained)) +
    theme_void(base_size = 10)

  official_diss <- eco_p |> group_by(across(all_of(t$field))) |>
    summarise(.groups = "drop") |> st_transform(4326)
  p3 <- ggplot(official_diss) +
    geom_sf(aes(fill = factor(.data[[t$field]])), color = "white", linewidth = 0.1) +
    scale_fill_manual(values = pal(t$K), guide = "none") +
    labs(title = sprintf("Official %ss (K=%d)", t$label, t$K),
         subtitle = "National Ecological Framework for Canada") +
    theme_void(base_size = 10)

  combined_plot <- p1 + p2 + p3 +
    plot_annotation(caption = sprintf(
      "ARI vs. official %ss — Ward (sampled): %.3f | GDM nearest-centroid (sampled): %.3f | GDM-filled (all cells): %.3f",
      t$label, ari_ward, ari_sampled, ari_filled))

  out_png <- file.path(OUT_DIR, sprintf("bioregions_inat2025_gdm_filled_map_%s.png", t$tag))
  ggsave(out_png, combined_plot, width = 20, height = 7.5, dpi = 200, bg = "white")
  cat(sprintf("  Saved: %s\n", out_png))
}

cat("\n── [6/6] Done. ─────────────────────────────────────────────\n")
