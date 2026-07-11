# Net change in suitable habitat extent (area-weighted)
library(terra)
library(fs)
library(dplyr)
library(tidyr)

# Paths
base_dir <- "./results_glossa"          # per-species GLOSSA projection outputs
out_dir  <- "./output/net_change"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Species, scenarios, decades (proj_1 .. proj_7)
species_list <- c(
  "Balaenoptera_acutorostrata", "Balaenoptera_borealis", "Balaenoptera_edeni",
  "Balaenoptera_musculus",      "Balaenoptera_physalus", "Delphinus_delphis",
  "Globicephala_macrorhynchus", "Grampus_griseus",       "Megaptera_novaeangliae",
  "Orcinus_orca",               "Physeter_macrocephalus","Pseudorca_crassidens",
  "Stenella_coeruleoalba",      "Stenella_frontalis",    "Tursiops_truncatus",
  "Ziphius_cavirostris"
)
ssps    <- c("ssp119", "ssp585")
decades <- c(2030, 2040, 2050, 2060, 2070, 2080, 2090)

# Per-species TSS-max thresholds defining suitable habitat
thresholds <- c(
  Balaenoptera_acutorostrata = 0.444134530113127,
  Balaenoptera_borealis      = 0.472731661493837,
  Balaenoptera_edeni         = 0.521948297923397,
  Balaenoptera_musculus      = 0.452278246331992,
  Balaenoptera_physalus      = 0.527919510923481,
  Delphinus_delphis          = 0.515739226839754,
  Globicephala_macrorhynchus = 0.492305050733093,
  Grampus_griseus            = 0.488290915917543,
  Megaptera_novaeangliae     = 0.509966205966828,
  Orcinus_orca               = 0.534303946162485,
  Physeter_macrocephalus     = 0.511006364633733,
  Pseudorca_crassidens       = 0.514510216988369,
  Stenella_coeruleoalba      = 0.496225572445452,
  Stenella_frontalis         = 0.504466217535132,
  Tursiops_truncatus         = 0.485043726590147,
  Ziphius_cavirostris        = 0.503008611615610
)

# Raster path helpers (GLOSSA output layout)
get_baseline <- function(sp_dir) {
  f <- dir_ls(file.path(sp_dir, "sh/proj/proj_baseline/mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  if (length(f) == 0) NA_character_ else f[1]
}
get_future <- function(sp_dir, ssp, k) {
  f <- dir_ls(file.path(sp_dir, "sh/proj", paste0("proj_", ssp), "mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  hit <- f[grepl(paste0("_", k, "_mean\\.tif$"), basename(f))]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Suitable habitat area (km^2): true geodesic cell area summed over suitable cells.
# Area weighting (not pixel count) corrects for cells shrinking towards the pole,
# which matters when suitable habitat shifts latitudinally.
area_suitable <- function(r, thr, area_r) {
  v <- values(r); a <- values(area_r)
  sum(a[!is.na(v) & v >= thr], na.rm = TRUE)
}

# Per-decade net change for each species x scenario
rows <- list()
for (sp in species_list) {
  sp_dir <- file.path(base_dir, sp)
  thr    <- thresholds[[sp]]
  bf <- get_baseline(sp_dir); if (is.na(bf)) next
  r_base    <- clamp(rast(bf), 0, 1)
  area_r    <- cellSize(r_base, unit = "km")          # latitude-corrected cell area
  base_area <- area_suitable(r_base, thr, area_r)

  for (ssp in ssps) {
    for (k in seq_along(decades)) {
      ff <- get_future(sp_dir, ssp, k); if (is.na(ff)) next
      r_fut    <- clamp(rast(ff), 0, 1)
      fut_area <- area_suitable(r_fut, thr, area_r)
      rows[[paste(sp, ssp, k)]] <- data.frame(
        species           = sp,
        ssp               = ssp,
        decade            = decades[k],
        baseline_area_km2 = base_area,
        future_area_km2   = fut_area,
        net_change        = (fut_area - base_area) / base_area   # +gain / -loss (proportion)
      )
    }
  }
}
traj <- bind_rows(rows)
write.csv(traj, file.path(out_dir, "net_change_perdecade.csv"), row.names = FALSE)

# End-century (proj_7 = 2090) net change, one row per species, wide by scenario
net_end <- traj %>%
  filter(decade == 2090) %>%
  select(species, ssp, net_change) %>%
  pivot_wider(names_from = ssp, values_from = net_change, names_prefix = "net_change_")
write.csv(net_end, file.path(out_dir, "net_change_endcentury.csv"), row.names = FALSE)
