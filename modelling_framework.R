# MODELLING FRAMEWORK

# --- LOAD NECESSARY LIBRARIES ---
library(terra)            
library(sf)             
library(dismo)            
library(biomod2)          
library(dplyr)            
library(blockCV)          
library(ggplot2)          
library(tools)         


# --- GLOBAL SETUP AND ENVIRONMENTAL DATA PREPARATION ---

# Load baseline environmental rasters (2010-2020)
n1 <- list.files(path = "path_to_baseline", pattern = "2010-2020.tif$", full.names = TRUE)
rasters <- rast(n1)
ref <- rasters[[1]]
# Align and stack rasters
aligned_rasters <- lapply(rasters, function(r) project(r, ref, method = "bilinear"))
aligned_stack <- rast(aligned_rasters)

# Crop to study extent
e <- ext(c(-35,10,0,70))
bioc <- crop(aligned_stack, e)

# Load kernel raster (sampling bias)
kernel_terra <- rast("path_to_kernel.tif")

# Crop and resample environmental data to match kernel raster extent and resolution
bioc_crop <- crop(bioc, ext(kernel_terra))
bioc_resampled <- resample(bioc_crop, kernel_terra, method = "bilinear")
bioc_mask <- mask(bioc_resampled, kernel_terra)

# Resample sampling bias layer to match environmental predictors
bias_raster_resampled <- resample(kernel_terra, bioc_mask[[1]], method = "bilinear")
names(bias_raster_resampled) <- "sampling_bias"

# Combine environmental predictors with sampling bias 
env_predictors <- c(bioc_mask, bias_raster_resampled)
env_stack_bias_for_bg <- raster::stack(env_predictors)
NAvalue(env_stack_bias_for_bg) <- -9999

# without sampling bias for modelling
keep_layers_for_modeling <- which(!(names(env_stack_bias_for_bg) %in% c("sampling_bias")))
env_stack_for_modeling <- subset(env_stack_bias_for_bg, keep_layers_for_modeling)
bioc_r_raster <- raster::stack(bioc)
NAvalue(bioc_r_raster) <- -9999


# --- LOAD FUTURE ENVIRONMENTAL DATA ---

# ssp 119 
preds_2040_s19 <- rast(list.files(path = "path_to_ssp119", pattern="2030-2040.tif$", full.names = TRUE)) # first layer, repeat for the remaining years


# Ensure future rasters have same names as baseline predictors and crop them

future_preds_raw <- list(
  ssp119_2040 = preds_2040_s19, ssp119_2050 = preds_2050_s19, ssp119_2060 = preds_2060_s19,
  ssp119_2070 = preds_2070_s19, ssp119_2080 = preds_2080_s19, ssp119_2090 = preds_2090_s19,
  ssp119_2100 = preds_2100_s19, ssp585_2040 = preds_2040_s85, ssp585_2050 = preds_2050_s85,    ssp585_2060 = preds_2060_s85, ssp585_2070 = preds_2070_s85, ssp585_2080 = preds_2080_s85,    ssp585_2090 = preds_2090_s85, ssp585_2100 = preds_2100_s85)

future_preds_raster <- list()
for (key in names(future_preds_raw)) {
  temp_rast <- future_preds_raw[[key]]
  names(temp_rast) <- names(bioc_r)
  temp_stack <- raster::stack(temp_rast) 
  NAvalue(temp_stack) <- -9999
  future_preds_raster[[key]] <- temp_stack
}

# Define scenario names and years for looping
scenarios <- c("ssp119", "ssp585")
years <- c("2040", "2050", "2060", "2070", "2080", "2090", "2100")


# --- LOAD AND PREPARE SPECIES DATA ---
sp_df <- read.csv("path_to_csv") %>%
  select(lon, lat, Species)
all_species_names <- unique(sp_df$Species)


# --- MAIN LOOP FOR EACH SPECIES ---
base_output_dir <- "./Output_SDM_Results" # Base directory for all results

for (species_name in all_species_names) {
  cat(paste("\n--- Starting modeling for species:", species_name, "---\n"))
  
  # Create species-specific output directory structure
  species_dir <- file.path(base_output_dir, species_name)
  
  tryCatch({
    sp_sub <- sp_df %>% filter(Species == species_name)
    
    # Create presence points matrix 
    presences <- sp_sub %>% select(lon, lat) %>% as.matrix()
    
    # Generate background points
    # Using the sampling bias raster to influence background point selection
    set.seed(123 + which(all_species_names == species_name)) 
    bg_points <- dismo::randomPoints(mask = env_stack_bias_for_bg[["sampling_bias"]], n = 10000)
    background <- bg_points
    
    
    # --- SPATIAL BLOCK CROSS-VALIDATION FOLDS ---
  
    # Convert presence and background points to sf
    pres_sf <- st_as_sf(data.frame(lon=presences[,1], lat=presences[,2]), coords = c("lon", "lat"), crs=4326)
    bg_sf <- st_as_sf(data.frame(lon=background[,1], lat=background[,2]), coords = c("lon", "lat"), crs=4326)
    
    # Combine for folds generation
    all_pts <- rbind(
      cbind(data.frame(st_coordinates(pres_sf)), Species_PA=1),
      cbind(data.frame(st_coordinates(bg_sf)), Species_PA=0)
    )
    
    all_sf <- st_as_sf(all_pts, coords = c("X", "Y"), crs = 4326)
    sf_filt <- all_sf %>% filter(Species_PA == 1)
    
    # Calculate spatial autocorrelation range for block size
      SAC_sp <- cv_spatial_autocor(
        x = sf_filt,
        column = "Species_PA",
        num_sample = min(5000, nrow(sf_filt)),
        deg_to_metre = 111325)
      block_size <- SAC_sp$range
    
    # Generate spatial blocks for cross-validation
    set.seed(42 + which(all_species_names == species_name)) 
    folds <- cv_spatial(x = all_sf,
                        column = "Species_PA", 
                        r = env_stack_for_modeling,
                        size = block_size,
                        k = 3, 
                        selection = "systematic",
                        plot = FALSE, 
                        deg_to_metre = 111325,
                        progress = FALSE) 
    
    # Add fold IDs to presences and background points
    all_sf$foldID <- folds$folds_ids
    
    # Combine presences and background with PA = 0 for background
    coords <- sf::st_coordinates(all_sf)
    
    # Add coordinates to the sf object
    all_sf$lon <- coords[, 1]
    all_sf$lat <- coords[, 2]
    
    pa_data <- data.frame(
      lon = all_sf$lon,
      lat = all_sf$lat,
      PA = all_sf$Species_PA,
      fold = all_sf$foldID
    )
    
    # Extract environmental values from raster
    env_values <- raster::extract(env_stack_for_modeling, pa_data[, c("lon", "lat")])
    
    # Combine environmental values with PA data
    pa_env_data <- cbind(pa_data, env_values)
    pa_env_data <- pa_env_data[complete.cases(pa_env_data), ]
    
    
    # --- BIOMOD2 DATA PREPARATION ---
    myRespName <- species_name
    myResp <- pa_env_data$PA
    myRespXY <- pa_env_data[, c("lon", "lat")]
    myExpl <- env_stack_for_modeling
    
    # To enable spatial CV in biomod2, we use the 'DataSplitTable' argument
    n_folds <- length(unique(pa_env_data$fold))
    n_points <- nrow(pa_env_data)
    
    # Initialize data split matrix
    DataSplitTable <- matrix(NA, nrow = n_points, ncol = n_folds)
    
    for(i in 1:n_folds) {
      DataSplitTable[, i] <- ifelse(pa_env_data$fold == i, FALSE, TRUE) 
    }
    colnames(DataSplitTable) <- paste0("_allData_RUN", 1:n_folds)
    
    # Add a column for "_allData_allRun" 
    num_observations <- nrow(pa_env_data)
    all_data_for_training_col <- rep(TRUE, num_observations)
    DataSplitTable_final <- cbind(DataSplitTable, '_allData_allRun' = all_data_for_training_col)
    
    # Prepare data object
    biomodData <- biomod2::BIOMOD_FormatingData(
      resp.var = myResp,
      expl.var = myExpl,
      resp.xy = myRespXY,
      resp.name = myRespName,
      PA.nb.rep = 1,
      PA.nb.absences = 10000
    )
    
    # --- MODELING OPTIONS ---
    
    user.MAXNET <- list('_allData_allRun' = list(
      regmult = 3,
      classes = "lh" # Linear, Hinge features
    ))
    
    # Combine them into a modeling options object
    biomodOptions <- bm_ModelingOptions(
      data.type = "binary",
      models = "MAXNET", 
      strategy = "user.defined",
      bm.format = biomodData, 
      user.val = list(
        GBM.binary.gbm.gbm = user.GBM,
        MAXNET.binary.maxnet.maxnet = user.MAXNET),
      user.base = "bigboss")
    
    # --- RUN MODELS ---
    biomodModelOut <- BIOMOD_Modeling(
      bm.format      = biomodData,
      models         = "MAXNET",
      OPT.user       = biomodOptions,
      CV.strategy    = "user.defined",
      CV.user.table  = DataSplitTable_final, 
      var.import     = 3, 
      metric.eval    = c("TSS","ROC", "BOYCE"), 
      modeling.id    = paste0(myRespName, "_modelling")
    )
    
    # --- SAVE MODEL OBJECT ---
    saveRDS(biomodModelOut, file.path(species_dir, "models", paste0("biomodModelOut_", species_name, ".rds")))
    
    # --- MODEL EVALUATION RESULTS ---
    models_scores <- get_evaluations(biomodModelOut)
    write.csv(models_scores, file.path(species_dir, "performance", paste0("performance_metrics_", species_name, ".csv")), row.names = FALSE)
    
    # Summarize mean and SD of performance metrics
    mean_sd_metrics <- models_scores %>%
      select(algo, metric.eval, calibration, validation) %>%
      group_by(algo, metric.eval) %>%
      summarise(
        Mean_Calibration = mean(calibration, na.rm = TRUE),
        SD_Calibration = sd(calibration, na.rm = TRUE),
        Mean_Validation = mean(validation, na.rm = TRUE),
        SD_Validation = sd(validation, na.rm = TRUE),
        .groups = 'drop'
      )

    # --- VARIABLE IMPORTANCE ---
    var_importance <- get_variables_importance(biomodModelOut)
    
    # --- RESPONSE CURVES ---
    mods_built <- get_built_models(biomodModelOut, run = 'allRun') 
    
    # Generate response curve data
    table_rc <- bm_PlotResponseCurves(
      bm.out = biomodModelOut,
      models.chosen = get_built_models(biomodModelOut),
      do.bivariate = FALSE,
      fixed.var = "mean",
      show.variables.responses = "all",
      save.file = "no" 
    )
    df_rc <- table_rc$tab
  
    df_summary_max <- df_rc %>%
      group_by(expl.name, expl.val) %>%
      summarise(
        mean = mean(pred.val),
        lower = quantile(pred.val, 0.025),
        upper = quantile(pred.val, 0.975),
        .groups = "drop")
    
    p1 <- ggplot(df_summary_max, aes(x = expl.val, y = mean)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80", alpha = 0.5) +
      geom_line(color = "black", linewidth = 0.7) +
      facet_wrap(~ expl.name, scales = "free_x") +
      theme_minimal() +
      labs(y = "Predicted suitability", x = "Environmental variable")

    
    # --- SPATIAL PREDICTION (BASELINE) ---
    proj_folder_baseline <- file.path(species_dir, "projections", "baseline")
    dir.create(proj_folder_baseline, recursive = TRUE, showWarnings = FALSE)
    
    biomodProj_baseline <- BIOMOD_Projection(
      bm.mod = biomodModelOut,
      proj.name = paste0(species_name, "_baseline_proj"),
      new.env = bioc_r_raster, 
      models.chosen = mods_built, 
      metric.binary = "TSS",
      build.clamping.mask = TRUE,
      output.format = ".tif")
  
    proj_raster_stack_baseline <- get_predictions(biomodProj_baseline)
  
    
    # --- SPATIAL PREDICTION (FUTURE) ---
    for (scenario in scenarios) {
      for (year in years) {
        
        key <- paste0(scenario, "_", year)
        future_env_raster_current <- future_preds_raster[[key]]
        names(future_env_raster_current) <- names(bioc_r_raster)
        
        # Define output folder for this specific future projection
        proj_folder_future_scenario <- file.path(species_dir, "projections", "future", scenario)
        dir.create(proj_folder_future_scenario, recursive = TRUE, showWarnings = FALSE)
        
        # Run BIOMOD Projection for future
        proj_obj_future <- BIOMOD_Projection(
          bm.mod = biomodModelOut,
          new.env = future_env_raster_current,
          proj.name = paste0(species_name, "_", key, "_proj"), 
          models.chosen = mods_built, 
          binary.meth = NULL,
          filtered.meth = NULL,
          compress = 'xz',
          output.format = '.tif'
        )
        
        proj_stack_future <- get_predictions(proj_obj_future)
      }
    } 
  } 


