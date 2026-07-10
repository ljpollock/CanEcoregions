## bioregionalization_inat2025.R
##
## Data-driven bioregionalization of Canada from 2025 iNaturalist occurrence
## records, using plant genera and animal species composition.
##
## Rationale for taxonomic resolution: plants are aggregated to GENUS because
## iNaturalist plant IDs are frequently coarse/uncertain at species level
## (especially for grasses, sedges, willows), while most animal groups
## (birds, mammals, herps, fish, many insects/arachnids/molluscs) are
## reliably identified to species and species is the ecologically
## meaningful unit for them.
##
## Method:
##   1. Pull 2025 iNaturalist records (Canada), research+needs_id grade only
##      (casual grade over-represents cultivated plants / captive animals).
##   2. Plant genus + animal species occurrence counts per grid cell
##      (50 km equal-area cells, Canada Albers Equal-Area Conic).
##   3. Keep cells with >=20 records and >=8 taxa ("well-sampled").
##   4. Hellinger-transform the site x taxon count matrix (Legendre &
##      Gallagher 2001) and cluster with Ward's method (Euclidean distance
##      on Hellinger-transformed data = Hellinger distance -> a proper
##      metric, appropriate for Ward clustering unlike raw Bray-Curtis).
##   5. Cut the dendrogram at k = 194 clusters to match the number of
##      terrestrial ecoregions in the National Ecological Framework for
##      Canada (Data/nef_ca_ter_ecoregion_v2_2.geojson / the Eckert et al.
##      framework used elsewhere in this project), so the emergent
##      "iNat bioregions" are directly comparable to the official regions.
##   6. Compare data-driven clusters to official ecoregions (majority
##      overlap + Adjusted Rand Index) and summarise characteristic taxa
##      per cluster.
##
## Outputs (Outputs/):
##   bioregions_inat2025_cells.csv   - per-cell assignment + stats
##   bioregions_inat2025_cells.geojson - cell polygons w/ cluster + best-match ecoregion
##   bioregions_inat2025_summary.csv - per-cluster summary + top indicator taxa
##   bioregions_inat2025_map.png     - iNat bioregions vs. official ecoregions

suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(data.table); library(sf)
  library(vegan); library(ggplot2); library(patchwork)
})

set.seed(1)

IN_PARQUET <- "Data/iNaturalist/inaturalist-canada-dec2025_smaller.parquet"
ECO_FILE   <- "Data/nef_ca_ter_ecoregion_v2_2.geojson"
OUT_DIR    <- "Outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

CELL_SIZE   <- 50000    # metres (50 km equal-area grid)
MIN_RECORDS <- 20       # min occurrence records for a cell to be "well-sampled"
MIN_TAXA    <- 8        # min distinct taxa (genera+species) for a cell
ANIMAL_TAXA <- c("Aves","Mammalia","Insecta","Arachnida","Mollusca",
                 "Amphibia","Reptilia","Actinopterygii")
CANADA_AEA  <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# ── 1. Number of target clusters = number of official ecoregions ──────────────
cat("── Loading official ecoregions ─────────────────────────────\n")
eco <- st_read(ECO_FILE, quiet = TRUE)
K <- length(unique(eco$ECOREGION_ID))
cat(sprintf("  %d official ecoregions (%s) -> target K = %d\n",
            K, basename(ECO_FILE), K))

# ── 2. Load & filter 2025 iNaturalist records ──────────────────────────────────
cat("── Loading 2025 iNaturalist records ────────────────────────\n")
ds <- open_dataset(IN_PARQUET)
d <- ds |>
  filter(year == "2025", quality_grade != "casual",
         (iconic_taxon_name == "Plantae") | (iconic_taxon_name %in% ANIMAL_TAXA)) |>
  select(longitude, latitude, scientific_name, iconic_taxon_name) |>
  collect()
setDT(d)
d <- d[!is.na(longitude) & !is.na(latitude) & !is.na(scientific_name) & scientific_name != ""]
cat(sprintf("  %s candidate records (plants + animals, non-casual)\n", format(nrow(d), big.mark=",")))

d[, n_words := lengths(strsplit(scientific_name, "\\s+"))]
d[, w1 := sub("^(\\S+).*", "\\1", scientific_name)]
d[, w2 := ifelse(n_words >= 2, sub("^\\S+\\s+(\\S+).*", "\\1", scientific_name), NA_character_)]
BAD_W2 <- c("sp","sp.","spp","spp.","cf","cf.","aff","aff.","complex","group","x","hybrid")

plants <- d[iconic_taxon_name == "Plantae"]
plants[, `:=`(taxon = paste0("Plant genus: ", w1), unit = w1, rank = "genus")]

animals <- d[iconic_taxon_name %in% ANIMAL_TAXA & n_words >= 2 &
             grepl("^[a-z]+\\.?$", w2) & !(tolower(w2) %in% BAD_W2)]
animals[, `:=`(taxon = paste0("Animal sp: ", w1, " ", w2), unit = paste(w1, w2), rank = "species")]

cat(sprintf("  Plant records:  %s  (%s genera)\n",
            format(nrow(plants), big.mark=","), format(uniqueN(plants$unit), big.mark=",")))
cat(sprintf("  Animal records: %s  (%s species; %s dropped as subspecies/unresolved)\n",
            format(nrow(animals), big.mark=","), format(uniqueN(animals$unit), big.mark=","),
            format(nrow(d[iconic_taxon_name %in% ANIMAL_TAXA]) - nrow(animals), big.mark=",")))

comb <- rbind(plants[, .(longitude, latitude, taxon, kingdom = "Plant")],
              animals[, .(longitude, latitude, taxon, kingdom = "Animal")])
cat(sprintf("  Combined: %s records, %s taxa (genera+species)\n",
            format(nrow(comb), big.mark=","), format(uniqueN(comb$taxon), big.mark=",")))

# ── 3. Assign to equal-area grid ───────────────────────────────────────────────
cat("── Building equal-area grid ────────────────────────────────\n")
pts <- st_as_sf(comb, coords = c("longitude","latitude"), crs = 4326, remove = FALSE)
pts_p <- st_transform(pts, CANADA_AEA)
xy <- st_coordinates(pts_p)
comb[, `:=`(x = xy[,1], y = xy[,2])]
comb[, `:=`(col = floor(x / CELL_SIZE), row = floor(y / CELL_SIZE))]
comb[, cell_id := paste(col, row, sep = "_")]

cell_stats <- comb[, .(n_records = .N, n_taxa = uniqueN(taxon),
                        n_plant_genera = uniqueN(taxon[kingdom == "Plant"]),
                        n_animal_species = uniqueN(taxon[kingdom == "Animal"])),
                    by = .(cell_id, col, row)]
well <- cell_stats[n_records >= MIN_RECORDS & n_taxa >= MIN_TAXA]
cat(sprintf("  %s grid cells (%d km) total, %s well-sampled (>=%d records, >=%d taxa)\n",
            format(nrow(cell_stats), big.mark=","), CELL_SIZE/1000,
            format(nrow(well), big.mark=","), MIN_RECORDS, MIN_TAXA))

# ── 4. Site x taxon matrix (well-sampled cells only) ───────────────────────────
cat("── Building community matrix ───────────────────────────────\n")
comb_well <- comb[cell_id %in% well$cell_id]
mat_dt <- dcast(comb_well, cell_id ~ taxon, fun.aggregate = length, value.var = "taxon", fill = 0)
cell_ids <- mat_dt$cell_id
mat <- as.matrix(mat_dt[, -1, with = FALSE])
rownames(mat) <- cell_ids
cat(sprintf("  Matrix: %d cells x %d taxa\n", nrow(mat), ncol(mat)))

# Hellinger transform (Legendre & Gallagher 2001): sqrt of within-row relative
# abundance. Euclidean distance on this = Hellinger distance, a proper metric
# well suited to Ward's minimum-variance clustering (unlike raw Bray-Curtis).
mat_hel <- decostand(mat, method = "hellinger")

# ── 5. Cluster into K bioregions ───────────────────────────────────────────────
cat(sprintf("── Clustering into K = %d bioregions ───────────────────────\n", K))
d_hel <- dist(mat_hel, method = "euclidean")
hc <- hclust(d_hel, method = "ward.D2")
clust <- cutree(hc, k = K)
well[, cluster := clust[match(cell_id, names(clust))]]
cat(sprintf("  Realised clusters: %d (some target clusters may be empty if K > distinct groups)\n",
            uniqueN(clust)))
cluster_sizes <- table(clust)
cat(sprintf("  Cluster size range: %d - %d cells (median %d)\n",
            min(cluster_sizes), max(cluster_sizes), stats::median(cluster_sizes)))

# ── 6. Build grid-cell polygons & compare to official ecoregions ──────────────
cat("── Building polygons & comparing to official ecoregions ────\n")
well[, `:=`(xmin = col * CELL_SIZE, xmax = (col + 1) * CELL_SIZE,
            ymin = row * CELL_SIZE, ymax = (row + 1) * CELL_SIZE)]

make_sq <- function(xmin, xmax, ymin, ymax) {
  st_polygon(list(matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin),
                          ncol = 2, byrow = TRUE)))
}
geoms <- mapply(make_sq, well$xmin, well$xmax, well$ymin, well$ymax, SIMPLIFY = FALSE)
cells_sf <- st_sf(well, geometry = st_sfc(geoms, crs = CANADA_AEA))
cells_sf <- st_transform(cells_sf, 4326)
cells_sf$cell_area_km2 <- as.numeric(st_area(st_transform(cells_sf, CANADA_AEA))) / 1e6

# Assign each cell to whichever official ecoregion contains its centroid
cent <- st_centroid(cells_sf)
eco_match <- suppressWarnings(st_join(cent, eco[, c("ECOREGION_ID","ECOREGION_NAME_EN")], join = st_intersects))
cells_sf$ECOREGION_ID <- eco_match$ECOREGION_ID
cells_sf$ECOREGION_NAME_EN <- eco_match$ECOREGION_NAME_EN

# Adjusted Rand Index between iNat clusters and official ecoregions
ari <- {
  ok <- !is.na(cells_sf$ECOREGION_ID)
  tab <- table(cells_sf$cluster[ok], cells_sf$ECOREGION_ID[ok])
  n <- sum(tab)
  a <- rowSums(tab); b <- colSums(tab)
  comb2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(comb2(tab))
  sum_a  <- sum(comb2(a))
  sum_b  <- sum(comb2(b))
  expected <- sum_a * sum_b / comb2(n)
  maxidx   <- (sum_a + sum_b) / 2
  (sum_ij - expected) / (maxidx - expected)
}
cat(sprintf("  Adjusted Rand Index (iNat clusters vs. official ecoregions): %.3f\n", ari))
cat("  (0 = no better than random grouping, 1 = identical partitions)\n")

# ── 7. Characteristic taxa per cluster (simple IndVal: specificity x fidelity) ─
cat("── Computing characteristic taxa per cluster ───────────────\n")
clust_of_cell <- setNames(well$cluster, well$cell_id)
grp <- clust_of_cell[rownames(mat)]
col_sums_by_grp <- rowsum(mat, group = grp)          # K x taxa: total count per group
grp_totals <- rowSums(col_sums_by_grp)                 # total records per group
specificity <- sweep(col_sums_by_grp, 1, grp_totals, "/")   # A: share of group's records
specificity <- sweep(specificity, 2, colSums(specificity), "/") # normalise across groups -> specificity to THIS group
pa <- (mat > 0) * 1
fidelity <- rowsum(pa, group = grp)
fidelity <- sweep(fidelity, 1, table(grp)[rownames(fidelity)], "/")  # B: fraction of group's cells with taxon
indval <- specificity * fidelity

top_taxa <- rbindlist(lapply(rownames(indval), function(g) {
  v <- indval[g, ]
  v <- v[fidelity[g, ] >= 0.3]                     # present in >=30% of the cluster's cells
  if (length(v) == 0) return(NULL)
  v <- sort(v, decreasing = TRUE)[seq_len(min(5, length(v)))]
  data.table(cluster = as.integer(g), taxon = names(v), indval = round(unname(v), 3))
}))
top_taxa_txt <- top_taxa[, .(top_taxa = paste(sprintf("%s (%.2f)", sub("^(Plant genus: |Animal sp: )", "", taxon), indval),
                                               collapse = "; ")), by = cluster]

# ── 8. Per-cluster summary ─────────────────────────────────────────────────────
clust_summary <- well[, .(n_cells = .N, n_records = sum(n_records)), by = cluster]

cent_xy <- st_coordinates(st_centroid(cells_sf))
cent_dt <- as.data.table(cent_xy)[, cluster := cells_sf$cluster]
centroid_by_clust <- cent_dt[, .(centroid_lon = mean(X), centroid_lat = mean(Y)), by = cluster]
clust_summary <- merge(clust_summary, centroid_by_clust, by = "cluster")

comb_well[, cluster := clust_of_cell[cell_id]]
tax_by_clust <- comb_well[, .(
  n_plant_genera = uniqueN(taxon[kingdom == "Plant"]),
  n_animal_species = uniqueN(taxon[kingdom == "Animal"])
), by = cluster]
clust_summary <- merge(clust_summary, tax_by_clust, by = "cluster")

# Majority-overlap official ecoregion for each cluster
best_eco <- as.data.table(st_drop_geometry(cells_sf))[!is.na(ECOREGION_NAME_EN),
  .N, by = .(cluster, ECOREGION_NAME_EN)][order(cluster, -N)][, .SD[1], by = cluster]
setnames(best_eco, c("ECOREGION_NAME_EN","N"), c("best_match_ecoregion","best_match_n_cells"))

clust_summary <- merge(clust_summary, best_eco, by = "cluster", all.x = TRUE)
clust_summary <- merge(clust_summary, top_taxa_txt, by = "cluster", all.x = TRUE)
setorder(clust_summary, -n_cells)
clust_summary[, area_km2 := n_cells * (CELL_SIZE/1000)^2]

# ── 9. Save outputs ─────────────────────────────────────────────────────────────
cat("── Saving outputs ──────────────────────────────────────────\n")
fwrite(as.data.table(st_drop_geometry(cells_sf))[, .(
  cell_id, cluster, n_records, n_taxa, n_plant_genera, n_animal_species,
  cell_area_km2, ECOREGION_ID, ECOREGION_NAME_EN
)], file.path(OUT_DIR, "bioregions_inat2025_cells.csv"))

st_write(cells_sf[, c("cell_id","cluster","n_records","n_taxa","n_plant_genera",
                       "n_animal_species","ECOREGION_ID","ECOREGION_NAME_EN")],
          file.path(OUT_DIR, "bioregions_inat2025_cells.geojson"), delete_dsn = TRUE, quiet = TRUE)

fwrite(clust_summary, file.path(OUT_DIR, "bioregions_inat2025_summary.csv"))

# ── 10. Comparison map ──────────────────────────────────────────────────────────
cat("── Rendering comparison map ────────────────────────────────\n")
pal <- colorRampPalette(c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                           "#80BA5A","#E68310","#008695","#CF1C90","#f97b72",
                           "#4b4b8f","#A5AA99"))(K)
cells_sf$cluster_f <- factor(cells_sf$cluster)

p1 <- ggplot(cells_sf) +
  geom_sf(aes(fill = cluster_f), color = NA) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = sprintf("iNaturalist 2025 bioregions (data-driven, K=%d)", K),
       subtitle = sprintf("%s well-sampled 50 km cells; plant genera + animal species composition",
                           format(nrow(cells_sf), big.mark=","))) +
  theme_void(base_size = 10)

eco_ll <- st_transform(eco, 4326)
p2 <- ggplot(eco_ll) +
  geom_sf(aes(fill = factor(ECOREGION_ID)), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = colorRampPalette(pal)(length(unique(eco_ll$ECOREGION_ID))), guide = "none") +
  labs(title = sprintf("Official terrestrial ecoregions (K=%d)", K),
       subtitle = "National Ecological Framework for Canada") +
  theme_void(base_size = 10)

combined_plot <- p1 + p2 +
  plot_annotation(caption = sprintf("Adjusted Rand Index vs. official ecoregions: %.3f", ari))

ggsave(file.path(OUT_DIR, "bioregions_inat2025_map.png"), combined_plot,
       width = 14, height = 7, dpi = 200, bg = "white")

cat(sprintf("\nDone.\n  Cells:   %s\n  Summary: %s\n  Map:     %s\n  ARI:     %.3f\n",
            file.path(OUT_DIR, "bioregions_inat2025_cells.geojson"),
            file.path(OUT_DIR, "bioregions_inat2025_summary.csv"),
            file.path(OUT_DIR, "bioregions_inat2025_map.png"), ari))
