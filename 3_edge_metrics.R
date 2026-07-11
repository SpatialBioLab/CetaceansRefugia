# Latitudinal edge shifts of suitable habitat (area-weighted)
library(terra)
library(fs)
library(dplyr)

# Paths
base_dir <- "./results_glossa"          # per-species GLOSSA projection outputs
out_dir  <- "./output/edge_metrics"
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
KM_PER_DEG_LAT <- 111.32

# Trailing / median / leading edge as area-weighted suitability quantiles
edges <- c(trailing = 0.05, median = 0.50, leading = 0.95)

species_label <- function(sp) {
  s <- gsub("_", " ", sp)
  if (s == "Balaenoptera edeni") "Balaenoptera edeni brydei" else s
}

# Raster helpers (GLOSSA output layout)
get_baseline <- function(sp_dir) {
  f <- dir_ls(file.path(sp_dir, "sh/proj/proj_baseline/mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  if (length(f) == 0) NULL else clamp(rast(f[1]), 0, 1)
}
get_future <- function(sp_dir, ssp, k) {
  f <- dir_ls(file.path(sp_dir, "sh/proj", paste0("proj_", ssp), "mean"), glob = "*.tif", fail = FALSE)
  f <- f[!grepl("_cat", basename(f))]
  hit <- f[grepl(paste0("_", k, "_mean\\.tif$"), basename(f))]
  if (length(hit) == 0) NULL else clamp(rast(hit[1]), 0, 1)
}

# Latitude at which the cumulative area-weighted suitability mass reaches q.
# Weighting by cell area avoids overweighting smaller high-latitude cells.
quantile_lat <- function(r, area_r, q) {
  df <- as.data.frame(c(r, area_r), xy = TRUE, na.rm = TRUE)
  names(df)[3:4] <- c("suit", "area")
  df <- df[df$suit > 0, ]; if (nrow(df) == 0) return(NA_real_)
  prof <- df %>%
    mutate(mass = suit * area) %>%
    group_by(y) %>% summarise(w = sum(mass), .groups = "drop") %>%
    arrange(y) %>% mutate(cumw = cumsum(w) / sum(w))
  if (q <= prof$cumw[1]) return(prof$y[1])
  hi <- which(prof$cumw >= q)[1]; lo <- hi - 1
  if (prof$cumw[hi] == prof$cumw[lo]) return(prof$y[hi])
  frac <- (q - prof$cumw[lo]) / (prof$cumw[hi] - prof$cumw[lo])
  prof$y[lo] + frac * (prof$y[hi] - prof$y[lo])
}

# Per-decade edge shifts (km) relative to baseline, for both scenarios
rows <- list()
for (sp in species_list) {
  rb <- get_baseline(file.path(base_dir, sp)); if (is.null(rb)) next
  area_r   <- cellSize(rb, unit = "km")
  base_lat <- sapply(edges, function(q) quantile_lat(rb, area_r, q))

  for (ssp in ssps) {
    for (k in seq_along(decades)) {
      rf <- get_future(file.path(base_dir, sp), ssp, k); if (is.null(rf)) next
      rf      <- clamp(resample(rf, rb, method = "bilinear"), 0, 1)
      fut_lat <- sapply(edges, function(q) quantile_lat(rf, area_r, q))
      rows[[paste(sp, ssp, k)]] <- data.frame(
        species           = species_label(sp),
        ssp               = ssp,
        decade            = decades[k],
        trailing_shift_km = (fut_lat["trailing"] - base_lat["trailing"]) * KM_PER_DEG_LAT,
        median_shift_km   = (fut_lat["median"]   - base_lat["median"])   * KM_PER_DEG_LAT,
        leading_shift_km  = (fut_lat["leading"]  - base_lat["leading"])  * KM_PER_DEG_LAT
      )
    }
  }
}
traj <- bind_rows(rows)

# Anchor the baseline (2020) as zero shift so rates are measured from baseline
anchor <- traj %>% distinct(species, ssp) %>%
  mutate(decade = 2020, trailing_shift_km = 0, median_shift_km = 0, leading_shift_km = 0)
traj <- bind_rows(anchor, traj) %>% arrange(species, ssp, decade)
write.csv(traj, file.path(out_dir, "edge_shifts_perdecade.csv"), row.names = FALSE)

# Shift rate (km/decade) as the OLS slope of edge shift over time; + = poleward
fit_rate <- function(g, col) {
  x <- (g$decade - 2020) / 10; y <- g[[col]]
  data.frame(
    rate_km_per_decade = round(coef(lm(y ~ x))[2], 1),
    traj_range_km      = round(max(y) - min(y)),
    fit_r              = round(suppressWarnings(cor(x, y)), 2)
  )
}
rate_rows <- list()
for (sp in unique(traj$species)) {
  for (ssp in ssps) {
    g <- traj %>% filter(species == sp, ssp == !!ssp) %>% arrange(decade)
    r <- data.frame(species = sp, ssp = ssp)
    for (e in names(edges)) {
      f <- fit_rate(g, paste0(e, "_shift_km")); names(f) <- paste0(e, "_", names(f))
      r <- bind_cols(r, f)
    }
    rate_rows[[paste(sp, ssp)]] <- r
  }
}
write.csv(bind_rows(rate_rows), file.path(out_dir, "edge_shift_rates.csv"), row.names = FALSE)
