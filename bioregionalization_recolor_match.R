## bioregionalization_recolor_match.R
##
## Re-renders the sparse / GDM-filled / official comparison maps for all
## three target levels (ecoregion K=194, ecozone K=15, ecoprovince K=53)
## with a shared colour scheme: each official region is coloured to match
## whichever predicted (GDM-filled) cluster overlaps it most, found via
## optimal one-to-one bipartite matching (Hungarian algorithm, clue::solve_LSAP)
## on the cluster x official-region cell-overlap contingency table. This makes
## visual comparison between "estimated" and "official" panels much easier —
## big, well-recovered regions end up the same colour in both panels.
##
## Reuses the cell-level assignments already saved by
## bioregionalization_inat2025_gdm_fill.R and
## bioregionalization_inat2025_gdm_fill_multiK.R (no need to redo the
## expensive community-matrix / NPP-extraction / GDM steps).
##
## Outputs (Outputs/):
##   bioregions_inat2025_gdm_filled_map_ecoregion_k194_matched.png
##   bioregions_inat2025_gdm_filled_map_ecozone_k15_matched.png
##   bioregions_inat2025_gdm_filled_map_ecoprovince_k53_matched.png

suppressPackageStartupMessages({
  library(sf); library(data.table); library(dplyr)
  library(ggplot2); library(patchwork); library(clue)
})

ECO_FILE <- "Data/nef_ca_ter_ecoregion_v2_2.geojson"
OUT_DIR  <- "Outputs"
eco <- st_read(ECO_FILE, quiet = TRUE)

pal_fn <- colorRampPalette(c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74",
                              "#80BA5A","#E68310","#008695","#CF1C90","#f97b72",
                              "#4b4b8f","#A5AA99"))

## Optimal cluster -> official-region colour matching (maximises total
## overlapping cells across all pairs). Returns a named vector: official
## region ID (as character) -> hex colour.
match_official_colors <- function(cluster_filled, official_id, pal) {
  cluster_filled <- as.character(cluster_filled)
  official_id    <- as.character(official_id)
  ok  <- !is.na(official_id) & !is.na(cluster_filled)
  tab <- table(cluster_filled[ok], official_id[ok])   # rows = clusters, cols = official regions
  clusters  <- rownames(tab)
  officials <- colnames(tab)
  n <- max(length(clusters), length(officials))
  cost <- matrix(0, n, n)
  cost[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  assign <- clue::solve_LSAP(cost, maximum = TRUE)   # assign[i] = matched column for row i

  cluster_pal <- setNames(pal, sort(as.integer(unique(cluster_filled))))
  official_color <- setNames(rep(NA_character_, length(officials)), officials)
  for (i in seq_len(length(clusters))) {
    j <- assign[i]
    if (j <= length(officials)) {
      official_color[officials[j]] <- cluster_pal[[clusters[i]]]
    }
  }
  # any official regions left unmatched (shouldn't normally happen when
  # n_clusters >= n_officials) get a neutral fallback grey
  official_color[is.na(official_color)] <- "#BBBBBB"
  official_color
}

render_target <- function(cells_csv, official_field, K, label, tag, geom_field_in_csv = "official_id") {
  cat(sprintf("── %s (K=%d) ──────────────────────────\n", label, K))
  cells <- fread(cells_csv)
  setnames(cells, geom_field_in_csv, "official_id", skip_absent = TRUE)
  cells[, official_id := as.character(official_id)]

  pal <- pal_fn(K)
  official_colors <- match_official_colors(cells$cluster_filled, cells$official_id, pal)

  # Rebuild sf for the filled cells (regenerate 50km squares from cell_id "col_row")
  CELL_SIZE  <- 50000
  CANADA_AEA <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  cells[, c("col","row") := tstrsplit(cell_id, "_", type.convert = TRUE)]
  make_sq <- function(col, row) {
    xmin <- col * CELL_SIZE; xmax <- xmin + CELL_SIZE
    ymin <- row * CELL_SIZE; ymax <- ymin + CELL_SIZE
    st_polygon(list(matrix(c(xmin,ymin, xmax,ymin, xmax,ymax, xmin,ymax, xmin,ymin),
                            ncol = 2, byrow = TRUE)))
  }
  geoms <- mapply(make_sq, cells$col, cells$row, SIMPLIFY = FALSE)
  cells_sf <- st_sf(cells, geometry = st_sfc(geoms, crs = CANADA_AEA))
  cells_sf <- st_transform(cells_sf, 4326)
  cells_sf$cluster_col <- pal[match(cells_sf$cluster_filled, sort(unique(cells$cluster_filled)))]

  sampled_only <- cells_sf[cells_sf$data_status == "sampled", ]
  p1 <- ggplot(sampled_only) +
    geom_sf(aes(fill = cluster_col), color = NA) +
    scale_fill_identity() +
    labs(title = sprintf("iNat 2025 bioregions (sparse, K=%d)", K),
         subtitle = sprintf("%s well-sampled cells", format(nrow(sampled_only), big.mark=","))) +
    theme_void(base_size = 10)

  p2 <- ggplot(cells_sf) +
    geom_sf(aes(fill = cluster_col), color = NA) +
    scale_fill_identity() +
    labs(title = sprintf("GDM-filled iNat 2025 bioregions (K=%d)", K),
         subtitle = sprintf("%s cells", format(nrow(cells_sf), big.mark=","))) +
    theme_void(base_size = 10)

  official_diss <- eco |> group_by(across(all_of(official_field))) |>
    summarise(.groups = "drop") |> st_transform(4326)
  official_diss$id_chr <- as.character(st_drop_geometry(official_diss)[[official_field]])
  official_diss$fill_col <- official_colors[official_diss$id_chr]
  official_diss$fill_col[is.na(official_diss$fill_col)] <- "#BBBBBB"

  p3 <- ggplot(official_diss) +
    geom_sf(aes(fill = fill_col), color = "white", linewidth = 0.1) +
    scale_fill_identity() +
    labs(title = sprintf("Official %ss (K=%d)", label, K),
         subtitle = "coloured to match best-overlapping GDM-filled cluster") +
    theme_void(base_size = 10)

  combined_plot <- p1 + p2 + p3 +
    plot_annotation(caption = sprintf(
      "Colours matched via optimal bipartite assignment (Hungarian algorithm) on cluster x %s cell overlap.",
      label))

  out_png <- file.path(OUT_DIR, sprintf("bioregions_inat2025_gdm_filled_map_%s_matched.png", tag))
  ggsave(out_png, combined_plot, width = 20, height = 7.5, dpi = 200, bg = "white")
  cat(sprintf("  Saved: %s\n", out_png))
}

render_target(file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells.csv"),
              "ECOREGION_ID", 194, "ecoregion", "ecoregion_k194",
              geom_field_in_csv = "ECOREGION_ID")

render_target(file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecozone_k15.csv"),
              "ECOZONE_ID", 15, "ecozone", "ecozone_k15")

render_target(file.path(OUT_DIR, "bioregions_inat2025_gdm_filled_cells_ecoprovince_k53.csv"),
              "ECOPROVINCE_ID", 53, "ecoprovince", "ecoprovince_k53")

cat("\nDone.\n")
