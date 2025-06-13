### refugia areas for each species

library(terra)
library(Kendall)

disp_s19 <- list.files(pattern = "disp_species_model_ssp119.tif")
disp_s85 <- list.files(pattern = "disp_species_model_ssp585.tif")

r_disp_s19 <- rast(disp_s19)
r_disp_s85 <- rast(disp_s85)
names(r_disp_s19) <- c("2040", "2050", "2060", "2070", "2080", "2090", "2100")
names(r_disp_s85) <- c("2040", "2050", "2060", "2070", "2080", "2090", "2100")

nodisp_s19 <- list.files(pattern = "nodisp_species_model_ssp119.tif")
nodisp_s85 <- list.files(pattern = "nodisp_species_model_sspp585.tif")

r_nodisp_s19 <- rast(nodisp_s19)
r_nodisp_s85 <- rast(nodisp_s85)
names(r_nodisp_s19) <- c("2040", "2050", "2060", "2070", "2080", "2090", "2100")
names(r_nodisp_s85) <- c("2040", "2050", "2060", "2070", "2080", "2090", "2100")


# MANN-KENDALL
# MK_disp_s19
MK_disp_s19 <- function(r_disp_s19, type=c("trend","pval","both")) {
  layers = nlyr(r_disp_s19)
  ncell = ncell(r_disp_s19)
  ncol = ncol(r_disp_s19)
  nrow = nrow(r_disp_s19)
  crs = crs(r_disp_s19)
  extent = ext(r_disp_s19)
  mtrx <- as.matrix(r_disp_s19)
  
  # Initialize an empty matrix to store results
  empt <- matrix(nrow = ncell, ncol = 2)
  
  # Start processing each pixel
  if (type == "trend") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
      }}
    
    # Create a SpatRaster to store trend values
    trend <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(trend) <- empt[, 1]
    trend} 
  else if (type == "pval") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
      }}
    
    # Create a SpatRaster to store p-values
    pval <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(pval) <- empt[, 1]
    pval}
  else if (type == "both") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else {
        tryCatch({
          empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
          empt[i, 2] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
        }, error = function(e) {
          empt[i, 1] <- NA
          empt[i, 2] <- NA
        })}}
    
    tr <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    pv <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(tr) <- empt[, 1]
    values(pv) <- empt[, 2]
    
    brk <- c(tr, pv)
    names(brk) <- c("trend", "p.value")
    brk}}

# trend slope and p-value
trend_disp_s19 <- MK_disp_s19(r_disp_s19, type = "trend")
pval_disp_s19 <- MK_disp_s19(r_disp_s19, type = "pval")


# MK_disp_s85
MK_disp_s85 <- function(r_disp_s85, type=c("trend","pval","both")) {
  layers = nlyr(r_disp_s85)
  ncell = ncell(r_disp_s85)
  ncol = ncol(r_disp_s85)
  nrow = nrow(r_disp_s85)
  crs = crs(r_disp_s85)
  extent = ext(r_disp_s85)
  mtrx <- as.matrix(r_disp_s85)
  
  # Initialize an empty matrix to store results
  empt <- matrix(nrow = ncell, ncol = 2)
  
  # Start processing each pixel
  if (type == "trend") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
      }}
    
    # Create a SpatRaster to store trend values
    trend <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(trend) <- empt[, 1]
    trend} 
  else if (type == "pval") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
      }}
    
    # Create a SpatRaster to store p-values
    pval <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(pval) <- empt[, 1]
    pval}
  else if (type == "both") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else {
        tryCatch({
          empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
          empt[i, 2] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
        }, error = function(e) {
          empt[i, 1] <- NA
          empt[i, 2] <- NA
        })}}
    
    tr <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    pv <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(tr) <- empt[, 1]
    values(pv) <- empt[, 2]
    
    brk <- c(tr, pv)
    names(brk) <- c("trend", "p.value")
    brk}}

# trend slope and p-value
trend_disp_s85 <- MK_disp_s85(r_disp_s85, type = "trend")
pval_disp_s85 <- MK_disp_s85(r_disp_s85, type = "pval")


# MK_nodisp_s19
MK_nodisp_s19 <- function(r_nodisp_s19, type=c("trend","pval","both")) {
  layers = nlyr(r_nodisp_s19)
  ncell = ncell(r_nodisp_s19)
  ncol = ncol(r_nodisp_s19)
  nrow = nrow(r_nodisp_s19)
  crs = crs(r_nodisp_s19)
  extent = ext(r_nodisp_s19)
  mtrx <- as.matrix(r_nodisp_s19)
  
  # Initialize an empty matrix to store results
  empt <- matrix(nrow = ncell, ncol = 2)
  
  # Start processing each pixel
  if (type == "trend") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
      }}
    
    # Create a SpatRaster to store trend values
    trend <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(trend) <- empt[, 1]
    trend} 
  else if (type == "pval") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
      }}
    
    # Create a SpatRaster to store p-values
    pval <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(pval) <- empt[, 1]
    pval}
  else if (type == "both") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else {
        tryCatch({
          empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
          empt[i, 2] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
        }, error = function(e) {
          empt[i, 1] <- NA
          empt[i, 2] <- NA
        })}}
    
    tr <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    pv <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(tr) <- empt[, 1]
    values(pv) <- empt[, 2]
    
    brk <- c(tr, pv)
    names(brk) <- c("trend", "p.value")
    brk}}

# trend slope and p-value
trend_nodisp_s19 <- MK_nodisp_s19(r_nodisp_s19, type = "trend")
pval_nodisp_s19 <- MK_nodisp_s19(r_nodisp_s19, type = "pval")


# MK_nodisp_s85
MK_nodisp_s85 <- function(r_nodisp_s85, type=c("trend","pval","both")) {
  layers = nlyr(r_nodisp_s85)
  ncell = ncell(r_nodisp_s85)
  ncol = ncol(r_nodisp_s85)
  nrow = nrow(r_nodisp_s85)
  crs = crs(r_nodisp_s85)
  extent = ext(r_nodisp_s85)
  mtrx <- as.matrix(r_nodisp_s85)
  
  # Initialize an empty matrix to store results
  empt <- matrix(nrow = ncell, ncol = 2)
  
  # Start processing each pixel
  if (type == "trend") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
      }}
    
    # Create a SpatRaster to store trend values
    trend <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(trend) <- empt[, 1]
    trend} 
  else if (type == "pval") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
      } else {
        empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
      }}
    
    # Create a SpatRaster to store p-values
    pval <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(pval) <- empt[, 1]
    pval}
  else if (type == "both") {
    for (i in 1:ncell) {
      if (all(is.na(mtrx[i,]))) { 
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else if (sum(!is.na(mtrx[i,])) < 4) {
        empt[i, 1] <- NA
        empt[i, 2] <- NA
      } else {
        tryCatch({
          empt[i, 1] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$tau)
          empt[i, 2] <- as.numeric(MannKendall(as.numeric(mtrx[i,]))$sl)
        }, error = function(e) {
          empt[i, 1] <- NA
          empt[i, 2] <- NA
        })}}
    
    tr <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    pv <- rast(nrows = nrow, ncols = ncol, crs = crs, ext = extent)
    values(tr) <- empt[, 1]
    values(pv) <- empt[, 2]
    
    brk <- c(tr, pv)
    names(brk) <- c("trend", "p.value")
    brk}}

# trend slope and p-value
trend_nodisp_s85 <- MK_nodisp_s85(r_nodisp_s85, type = "trend")
pval_nodisp_s85 <- MK_nodisp_s85(r_nodisp_s85, type = "pval")



# Stable areas: Non-significant p-values
stable1 <- (pval_disp_s19 >= 0.05)
stable2 <- (pval_disp_s85 >= 0.05)
stable3 <- (pval_nodisp_s19 >= 0.05)
stable4 <- (pval_nodisp_s85 >= 0.05)

# Refuge = Increasing + Stable areas
refuge1 <- (trend_disp_s19 > 0 & pval_disp_s19 < 0.05) | stable1
refuge2 <- (trend_disp_s85 > 0 & pval_disp_s85 < 0.05) | stable2
refuge3 <- (trend_nodisp_s19 > 0 & pval_nodisp_s19 < 0.05) | stable3
refuge4 <- (trend_nodisp_s85 > 0 & pval_nodisp_s85 < 0.05)

r2 <- resample(refuge2, refuge1, method = "near")
r3 <- resample(refuge3, refuge1, method = "near")
r4 <- resample(refuge4, refuge1, method = "near")

# Combine across scenarios
refuge_count <- refuge1 + r2 + r3 + r4

pa <- rast("path_to_binary_raster_of_baseline.tif")

pa <- resample(pa, refuge_count, method = "near")
high_suitability <- pa >= 1 
refuge_filtered <- refuge_count * high_suitability  # Mask refugia using high suitability

categories_filtered <- classify(refuge_filtered, rcl = matrix(c(
  0, 0, 5,  # No refuge (Category 5)
  1, 1, 4,  # Refuge under 1 scenario (Category 4)
  2, 2, 3,  # Refuge under 2 scenarios (Category 3)
  3, 3, 2,  # Refuge under 3 scenarios (Category 2)
  4, 4, 1   # Refuge under all 4 scenarios (Category 1)
), ncol = 3, byrow = TRUE))

plot(categories_filtered)
