# Habitat-based refugia: continuous in-situ index and refugial retention
library(terra)
library(fs)
library(dplyr)

# Paths
base_dir <- "./results_glossa"          # per-species GLOSSA projection outputs
out_dir  <- "./output/refugia"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Species and scenarios
species_list <- c(
  "Balaenoptera_acutorostrata", "Balaenoptera_borealis", "Balaenoptera_edeni",
  "Balaenoptera_musculus",      "Balaenoptera_physalus", "Delphinus_delphis",
  "Globicephala_macrorhynchus", "Grampus_griseus",       "Megaptera_novaeangliae",
  "Orcinus_orca",               "Physeter_macrocephalus","Pseudorca_crassidens",
  "Stenella_coeruleoalba",      "Stenella_frontalis",    "Tursiops_truncatus",
  "Ziphius_cavirostris"
)
ssps <- c("ssp119", "ssp585")

# Per-species TSS-max thresholds defining baseline suitable habitat
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

# Raster helpers (GLOSSA output layout); future = end-century surface (proj_7)
get_baseline <- function(sp_dir) {
  f <- dir_ls(file.path(sp_dir, "sh/proj/proj_baseline/mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  if (length(f) == 0) NULL else clamp(rast(f[1]), 0, 1)
}
get_future <- function(sp_dir, ssp) {
  f <- dir_ls(file.path(sp_dir, "sh/proj", paste0("proj_", ssp), "mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  f <- f[grepl("_7_mean\\.tif$", basename(f))]
  if (length(f) == 0) NULL else clamp(rast(f[1]), 0, 1)
}

summary_rows <- list()
for (sp in species_list) {
  sp_dir <- file.path(base_dir, sp)
  thr    <- thresholds[[sp]]
  r_base <- get_baseline(sp_dir); if (is.null(r_base)) next

  # Refugia index: baseline suitability x end-century suitability, per scenario
  insitu <- list()
  for (ssp in ssps) {
    r_fut <- get_future(sp_dir, ssp); if (is.null(r_fut)) next
    r_fut <- clamp(resample(r_fut, r_base, method = "bilinear"), 0, 1)
    insitu[[ssp]] <- r_base * r_fut
  }
  if (length(insitu) < length(ssps)) next

  # Summed over the two scenarios (bounded 0-2): high where habitat is suitable
  # now AND remains suitable under future climate. One raster per species.
  refugia_index <- Reduce(`+`, insitu)
  writeRaster(refugia_index, file.path(out_dir, paste0(sp, "_refugia_index.tif")), overwrite = TRUE)

  # Restrict to baseline suitable habitat so the metric describes persistence
  # within the species' occupied range rather than across the full study extent.
  in_range   <- r_base >= thr
  index_vals <- values(mask(refugia_index, in_range, maskvalues = FALSE), na.rm = TRUE)

  # Refugial retention: index relative to a no-change benchmark (baseline x
  # baseline summed over the two scenarios = 2 * baseline^2). 100% = suitability
  # fully retained within baseline habitat; <100% = erosion.
  benchmark_vals <- values(mask(2 * r_base^2, in_range, maskvalues = FALSE), na.rm = TRUE)

  per_ssp <- sapply(insitu, function(x)
    mean(values(mask(x, in_range, maskvalues = FALSE), na.rm = TRUE)))

  summary_rows[[sp]] <- data.frame(
    species                = sp,
    n_range_cells          = length(index_vals),
    mean_refugia_index     = round(mean(index_vals), 4),
    p90_refugia_index      = round(quantile(index_vals, 0.90), 4),
    mean_ssp119            = round(per_ssp[["ssp119"]], 4),
    mean_ssp585            = round(per_ssp[["ssp585"]], 4),
    refugial_retention_pct = round(100 * mean(index_vals) / mean(benchmark_vals), 1)
  )
}

write.csv(bind_rows(summary_rows),
          file.path(out_dir, "refugia_summary.csv"), row.names = FALSE)
