## bioregionalization_inat2025_gdm_fill.R
##
## Extends bioregionalization_inat2025.R: instead of leaving grid cells with
## insufficient iNaturalist sampling blank, fill them in by extrapolating the
## community-composition pattern using a Generalized Dissimilarity Model
## (GDM; Ferrier et al. 2007) fitted on the well-sampled cells, driven by
## environmental layers available in Data/ (MODIS MOD17A2H PsnNet, a NPP
## proxy) plus geographic distance.
##
## Steps:
##   1. Rebuild the 2025 plant-genus + animal-species community matrix and
##      Ward/Hellinger clustering (same as bioregionalization_inat2025.R) to
##      get 194 reference bioregions for the well-sampled cells.
##   2. Build the FULL 50 km equal-area land grid over the extent of the 194
##      official ecoregions (not just the occurrence bounding box) -> this
##      includes many cells with zero/insufficient iNat records.
##   3. Extract mean & SD of 2025 NPP (PsnNet) per grid cell (all cells) via
##      zonal statistics against the MODIS sinusoidal grid.
##   4. Fit a GDM: Bray-Curtis compositional dissimilarity (well-sampled
##      cells only) ~ f(geographic distance, NPP mean, NPP seasonality).
##   5. Use the fitted GDM to transform every grid cell's environment
##      (including blank cells) into GDM space, where Euclidean distance
##      approximates predicted compositional dissimilarity.
##   6. Assign every grid cell to its nearest cluster centroid in GDM space
##      (nearest-centroid extrapolation) -> gap-filled bioregion map.
##
## Outputs (Outputs/):
##   bioregions_inat2025_gdm_filled_cells.geojson / .csv
##   bioregions_inat2025_gdm_filled_map.png  (sparse vs. GDM-filled vs. official)

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

# ═══ 1. Rebuild well-sampled community matrix + Ward clusters ════════════════
cat("── [1/6] Rebuilding 2025 iNat community matrix ─────────────\n")
eco <- st_read(ECO_FILE, quiet = TRUE)
K <- length(unique(eco$ECOREGION_ID))
cat(sprintf("  Target K = %d (official ecoregions)\n", K))

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
mat <- mat[well$cell_id, ]   # ensure row order matches `well`
cat(sprintf("  Community matrix: %d cells x %d taxa\n", nrow(mat), ncol(mat)))

mat_hel <- decostand(mat, method = "hellinger")
hc <- hclust(dist(mat_hel, method = "euclidean"), method = "ward.D2")
clust <- cutree(hc, k = K)
well[, cluster := clust[match(cell_id, names(clust))]]
cat(sprintf("  Realised clusters: %d\n", uniqueN(clust)))

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
cat("── [3/6] Extracting 2025 NPP per grid cell (this takes a few minutes) ──\n")
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

# ═══ 4. Fit GDM on well-sampled cells ══════════════════════════════════════════
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
cat("  Predictor importance (summed I-spline coefficients):\n")
print(round(mod$coefficients, 4))

# ═══ 5. Transform all grid cells into GDM space ════════════════════════════════
cat("── [5/6] Extrapolating to all grid cells via GDM transform ──\n")
full_valid <- full_sf[!is.na(full_sf$npp_mean) & !is.na(full_sf$npp_sd), ]
trans_all <- gdm.transform(mod, data = data.frame(
  X = full_valid$x, Y = full_valid$y,
  npp_mean = full_valid$npp_mean, npp_sd = full_valid$npp_sd))
trans_all <- as.data.table(trans_all)

trans_well <- trans_all[match(well$cell_num, full_valid$cell_num)]
centroids <- trans_well[, lapply(.SD, mean), by = well$cluster]
setnames(centroids, "well", "cluster")
cent_mat <- as.matrix(centroids[, -1])
rownames(cent_mat) <- centroids$cluster

trans_mat <- as.matrix(trans_all)
nearest <- apply(trans_mat, 1, function(v) {
  dd <- sqrt(rowSums(sweep(cent_mat, 2, v)^2))
  as.integer(names(which.min(dd)))
})
full_valid$cluster_filled <- nearest
full_valid$data_status <- ifelse(full_valid$cell_id %in% well$cell_id, "sampled", "gap-filled")
cat(sprintf("  %s cells assigned (%s originally sampled, %s gap-filled)\n",
            format(nrow(full_valid), big.mark=","),
            format(sum(full_valid$data_status == "sampled"), big.mark=","),
            format(sum(full_valid$data_status == "gap-filled"), big.mark=",")))

full_valid <- st_transform(full_valid, 4326)
cent_ll <- st_centroid(full_valid)
eco_match <- suppressWarnings(st_join(cent_ll, eco[, c("ECOREGION_ID","ECOREGION_NAME_EN")], join = st_intersects))
full_valid$ECOREGION_ID <- eco_match$ECOREGION_ID
full_valid$ECOREGION_NAME_EN <- eco_match$ECOREGION_NAME_EN

ari_fn <- function(cl, ec) {
  ok <- !is.na(ec)
  tab <- table(cl[ok], ec[ok])
  n <- sum(tab); a <- rowSums(tab); b <- colSums(tab)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab)); sum_a <- sum(comb2(a)); sum_b <- sum(comb2(b))
  expected <- sum_a * sum_b / comb2(n); maxidx <- (sum_a + sum_b) / 2
  (sum_ij - expected) / (maxidx - expected)
}
ari_filled <- ari_fn(full_valid$cluster_filled, full_valid$ECOREGION_ID)
ari_gdm_sampled <- ari_fn(full_valid$cluster_filled[full_valid$data_status == "sampled"],
                           full_valid$ECOREGION_ID[full_valid$data_status == "sampled"])
well_eco_id <- full_valid$ECOREGION_ID[match(well$cell_id, full_valid$cell_id)]
ari_ward_sampled <- ari_fn(well$cluster, well_eco_id)
cat(sprintf("  ARI (original Ward clusters, sampled cells only):    %.3f\n", ari_ward_sampled))
cat(sprintf("  ARI (GDM nearest-centroid, sampled cells only):      %.3f\n", ari_gdm_sampled))
cat(sprintf("  ARI (GDM-filled, ALL cells incl. gap-filled):        %.3f\n", ari_filled))

# ═══ 6. Save outputs & render map ══════════════════════════════════════════════
cat("── [6/6] Saving outputs & rendering map ─────────────────────\n")
fwrite(as.data.table(st_drop_geometry(full_valid))[, .(
  cell_id, cluster_filled, data_status, npp_mean, npp_sd,
  ECOREGION_ID, ECOREGION_NAME_EN
)], file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.csv"))

st_write(full_valid[, c("cell_id","cluster_filled","data_status","npp_mean","npp_sd",
                         "ECOREGION_ID","ECOREGION_NAME_EN")],
          file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.geojson"),
          delete_dsn = TRUE, quiet = TRUE)

pal <- colorRampPalette(c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                           "#80BA5A","#E68310","#008695","#CF1C90","#f97b72",
                           "#4b4b8f","#A5AA99"))(K)

sampled_only_sf <- full_valid[full_valid$data_status == "sampled", ]
p1 <- ggplot(sampled_only_sf) +
  geom_sf(aes(fill = factor(cluster_filled)), color = NA) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "iNat 2025 bioregions (sparse, well-sampled cells only)",
       subtitle = sprintf("%s cells", format(nrow(sampled_only_sf), big.mark=","))) +
  theme_void(base_size = 10)

p2 <- ggplot(full_valid) +
  geom_sf(aes(fill = factor(cluster_filled)), color = NA) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "GDM-filled iNat 2025 bioregions",
       subtitle = sprintf("%s cells; NPP + geographic distance GDM (%.0f%% dev. explained)",
                           format(nrow(full_valid), big.mark=","), mod$explained)) +
  theme_void(base_size = 10)

eco_ll <- st_transform(eco, 4326)
p3 <- ggplot(eco_ll) +
  geom_sf(aes(fill = factor(ECOREGION_ID)), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = colorRampPalette(pal)(length(unique(eco_ll$ECOREGION_ID))), guide = "none") +
  labs(title = "Official terrestrial ecoregions", subtitle = "National Ecological Framework for Canada") +
  theme_void(base_size = 10)

combined_plot <- p1 + p2 + p3 +
  plot_annotation(caption = sprintf(
    "ARI vs. official ecoregions — Ward clusters (sampled cells): %.3f | GDM nearest-centroid (sampled cells): %.3f | GDM-filled (all cells): %.3f",
    ari_ward_sampled, ari_gdm_sampled, ari_filled))

ggsave(file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_map.png"), combined_plot,
       width = 20, height = 7.5, dpi = 200, bg = "white")

cat(sprintf("\nDone.\n  Cells:   %s\n  Map:     %s\n",
            file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.geojson"),
            file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_map.png")))
