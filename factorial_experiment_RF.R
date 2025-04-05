# Random Forest Factor Contribution Analysis for NDVI
# This script uses random forest models to analyze the contribution of different environmental factors to NDVI

# Install required packages if needed
# install.packages(c( "raster", "rgdal", "readxl",  "randomForest"))

# Load required libraries
library(raster)
library(rgdal)
library(readxl)
library(randomForest)

# Set paths to input data
ndvi_path <- "/home/you_ch/NDVI/max_grass_0.1"
precipitation_path <- "/home/you_ch/ERA5/PPT/grass"
temperature_path <- "/home/you_ch/ERA5/TA/grass"
radiation_path <- "/home/you_ch/ERA5/Rad/grass"
nitrogen_path <- "/home/you_ch/ERA5/N/grass"
cropland_path <- "/home/you_ch/landuse/crop_0.1/grass"
rangeland_path <- "/home/you_ch/landuse/range_0.1/grass"
co2_excel_path <- "/home/you_ch/ERA5/CO2.xlsx"

# Function to read and stack raster files
read_variable_files <- function(path, pattern) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  stack <- stack(files)
  return(as.array(stack))
}

# Read NDVI data
ndvi_files <- list.files(ndvi_path, pattern = "cropped_MAX_NDVI_global_0.1_.*\\.tif$", full.names = TRUE)
ndvi_stack <- stack(ndvi_files)
n_years <- nlayers(ndvi_stack)
ndvi_array <- as.array(ndvi_stack)

# Read environmental variables
precipitation_array <- read_variable_files(precipitation_path, "cropped_EAR5_PPT_.*\\.tiff$")
temperature_array <- read_variable_files(temperature_path, "cropped_EAR5_TA_.*\\.tiff$")
radiation_array <- read_variable_files(radiation_path, "cropped_EAR5_Rad_.*\\.tiff$")
nitrogen_array <- read_variable_files(nitrogen_path, "cropped_trendy_N_.*\\.tiff$")
cropland_array <- read_variable_files(cropland_path, "cropped_crop_0.1_.*\\.tif$")
rangeland_array <- read_variable_files(rangeland_path, "cropped_range_0.1_.*\\.tif$")

# Read and process CO2 data
co2_data <- read_excel(co2_excel_path)
co2_values <- co2_data$CO2
co2_array <- array(rep(co2_values, each = ncell(ndvi_stack)), 
                   dim = c(dim(ndvi_stack)[1], dim(ndvi_stack)[2], n_years))

# Initialize arrays to store predictions for all scenarios
RF_ALL <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))
RF_CLI <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))
RF_CO2 <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))
RF_N <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))
RF_Crop <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))
RF_Range <- array(NA, dim = c(dim(ndvi_array)[1], dim(ndvi_array)[2], 40))

# Main analysis loop
for (i in 1:dim(ndvi_array)[1]) {
  for (j in 1:dim(ndvi_array)[2]) {
    print(paste("Processing pixel - Row:", i, "Column:", j, 
                "Total rows:", dim(ndvi_array)[1], "Total columns:", dim(ndvi_array)[2]))
    
    ndvi_values <- ndvi_array[i, j, ]
    
    # Skip if no valid NDVI data
    if (all(is.na(ndvi_values)) || all(ndvi_values == 0) || length(unique(ndvi_values)) == 1) next
    
    # Extract environmental variables
    crop_values <- cropland_array[i, j, ]
    range_values <- rangeland_array[i, j, ]
    co2_values <- co2_array[i, j, ]
    n_values <- nitrogen_array[i, j, ]
    ta_values <- temperature_array[i, j, ]
    rad_values <- radiation_array[i, j, ]
    ppt_values <- precipitation_array[i, j, ]
    years <- 1982:2021
    
    # Create data frame
    var_values <- data.frame(
      Crop = crop_values,
      CO2 = co2_values,
      N = n_values,
      TA = ta_values,
      Rad = rad_values,
      PPT = ppt_values,
      Range = range_values,
      Year = years
    )
    
    # Remove rows with missing values
    cordata <- data.frame(NDVImax = ndvi_values, var_values)
    cordata <- cordata[complete.cases(cordata), ]
    
    # Only proceed if sufficient data (>=20 observations)
    if (nrow(cordata) >= 20) {
      set.seed(123) # For reproducibility
      
      # Define parameter grid for tuning (only mtry in this case)
      tune_grid <- expand.grid(
        mtry = c(1, 2, 3, 4, 5)
      )
      
      # Set up cross-validation
      train_control <- trainControl(method = "cv", number = 5, 
                                    verboseIter = FALSE, 
                                    savePredictions = "all")
      
      # Train random forest model with automatic parameter tuning
      rf_model <- train(NDVImax ~ Crop + CO2 + N + TA + Rad + PPT + Range,
                        data = cordata,
                        method = "rf",
                        trControl = train_control,
                        ntree = 500,
                        nodesize = 5,
                        maxnodes = 20,
                        tuneGrid = tune_grid, 
                        importance = TRUE)
      
      # Make predictions and store results
      predictions <- predict(rf_model, newdata = cordata)
      
      # Get valid year indices
      valid_years <- cordata$Year
      
      # Store predictions for each scenario
      for (k in 1:length(valid_years)) {
        year_index <- which(years == valid_years[k])
        
        # Full model predictions
        RF_ALL[i, j, year_index] <- predictions[k]
        
        # Climate-fixed scenario (PPT, Rad, TA fixed to first year)
        cli_data <- cordata
        cli_data$PPT <- ppt_values[1]
        cli_data$Rad <- rad_values[1]
        cli_data$TA <- ta_values[1]
        RF_CLI[i, j, year_index] <- predict(rf_model, newdata = cli_data)[k]
        
        # CO2-fixed scenario
        co2_data <- cordata
        co2_data$CO2 <- co2_values[1]
        RF_CO2[i, j, year_index] <- predict(rf_model, newdata = co2_data)[k]
        
        # Nitrogen-fixed scenario
        n_data <- cordata
        n_data$N <- n_values[1]
        RF_N[i, j, year_index] <- predict(rf_model, newdata = n_data)[k]
        
        # Cropland-fixed scenario
        crop_data <- cordata
        crop_data$Crop <- crop_values[1]
        RF_Crop[i, j, year_index] <- predict(rf_model, newdata = crop_data)[k]
        
        # Rangeland-fixed scenario
        range_data <- cordata
        range_data$Range <- range_values[1]
        RF_Range[i, j, year_index] <- predict(rf_model, newdata = range_data)[k]
      }
    }
  }
}

# Save results as raster files
years <- 1982:2021
output_folder <- "/home/you_ch/RF/factor experiment/"

# Read reference raster for projection information
refer_ras <- raster("/home/you_ch/NDVI/max_grass_0.1/cropped_MAX_NDVI_global_0.1_1982.tif")

# Function to create properly formatted raster
create_raster <- function(data) {
  ras <- raster(nrow = 1730, ncol = 3600, 
                xmn = -180, xmx = 180, 
                ymn = -89.9, ymx = 83.1, 
                crs = CRS("+proj=longlat +datum=WGS84"))
  values(ras) <- as.matrix(data)
  return(ras)
}

# Save yearly predictions for each scenario
for (year in years) {
  year_idx <- which(years == year)
  
  # Extract data for current year
  year_rf_all <- RF_ALL[,,year_idx]
  year_rf_cli <- RF_CLI[,,year_idx]
  year_rf_co2 <- RF_CO2[,,year_idx]
  year_rf_n <- RF_N[,,year_idx]
  year_rf_crop <- RF_Crop[,,year_idx]
  year_rf_range <- RF_Range[,,year_idx]
  
  # Create rasters
  rf_all_raster <- create_raster(year_rf_all)
  rf_cli_raster <- create_raster(year_rf_cli)
  rf_co2_raster <- create_raster(year_rf_co2)
  rf_n_raster <- create_raster(year_rf_n)
  rf_crop_raster <- create_raster(year_rf_crop)
  rf_range_raster <- create_raster(year_rf_range)
  
  # Construct filenames
  file_suffix <- paste0("_", year, ".tif")
  
  # Save rasters
  writeRaster(rf_all_raster, 
              filename = paste0(output_folder, "RF_ALL_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  writeRaster(rf_cli_raster, 
              filename = paste0(output_folder, "RF_CLI_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  writeRaster(rf_co2_raster, 
              filename = paste0(output_folder, "RF_CO2_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  writeRaster(rf_n_raster, 
              filename = paste0(output_folder, "RF_N_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  writeRaster(rf_crop_raster, 
              filename = paste0(output_folder, "RF_Crop_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  writeRaster(rf_range_raster, 
              filename = paste0(output_folder, "RF_Range_3", file_suffix), 
              format = "GTiff", NAflag = -9999, overwrite = TRUE)
  
  print(paste("Year", year, "saved"))
}

print("Analysis complete! All results saved.")