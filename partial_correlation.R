# Partial Correlation Analysis for NDVI Drivers
# This script calculates partial correlations between NDVI and environmental variables
# while controlling for temporal trends and other variables

# Install required packages if needed
# install.packages(c( "raster", "rgdal", "readxl", "ppcor"))

# Load required libraries
library(raster)
library(rgdal)
library(readxl)
library(ppcor)

# Set paths to input data
ndvi_path <- "/data/NDVI"
precipitation_path <- "/data/PPT"
temperature_path <- "/data/TA"
radiation_path <- "/data/Rad"
nitrogen_path <- "/data/N"
cropland_path <- "/data/Crop"
rangeland_path <- "/data/Range"
co2_excel_path <- "/data/CO2/CO2.xlsx"

# Output directory for results
output_dir <- "data/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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

# Initialize output rasters
ta_pcor_raster <- ppt_pcor_raster <- rad_pcor_raster <- co2_pcor_raster <- 
  n_pcor_raster <- crop_pcor_raster <- range_pcor_raster <- raster(ndvi_stack[[1]])

p_ta_pcor_raster <- p_ppt_pcor_raster <- p_rad_pcor_raster <- p_co2_pcor_raster <- 
  p_n_pcor_raster <- p_crop_pcor_raster <- p_range_pcor_raster <- raster(ndvi_stack[[1]])



# Main analysis loop
for (i in 1:dim(ndvi_array)[1]) {
  for (j in 1:dim(ndvi_array)[2]) {
    # Print progress
    if (j == 1) {
      print(paste("Processing row", i, "of", dim(ndvi_array)[1], 
                  "(", round(i/dim(ndvi_array)[1]*100, 1), "%)"))
    }
    
    # Extract NDVI time series
    ndvi_values <- ndvi_array[i, j, ]
    
    # Skip if no data or constant values
    if (all(is.na(ndvi_values)) || all(ndvi_values == 0) || 
        length(unique(ndvi_values)) == 1) next
    
    # Extract environmental variables
    env_data <- data.frame(
      NDVI = ndvi_values,
      Crop = cropland_array[i, j, ],
      Range = rangeland_array[i, j, ],
      CO2 = co2_array[i, j, ],
      N = nitrogen_array[i, j, ],
      TA = temperature_array[i, j, ],
      Rad = radiation_array[i, j, ],
      PPT = precipitation_array[i, j, ],
      Year = 1982:2021
    )
    
    # Remove rows with missing values
    env_data <- env_data[complete.cases(env_data), ]
    
    # Skip if insufficient data
    if (nrow(env_data) < 20) next
    
    # Detrend variables
    detrend_var <- function(y, x) {
      model <- lm(y ~ x)
      return(resid(model))
    }
    
    # Create detrended data frame
    detrended <- data.frame(NDVI = detrend_var(env_data$NDVI, env_data$Year))
    
    # Detrend each environmental variable
    for (var in c("Crop", "Range", "CO2", "N", "TA", "Rad", "PPT")) {
      if (length(unique(env_data[[var]])) > 1) {  # Only process if variable varies
        detrended[[var]] <- detrend_var(env_data[[var]], env_data$Year)
      }
    }
    
    # Remove near-constant variables
    sd_check <- sapply(detrended, function(x) sd(x, na.rm = TRUE) > threshold)
    detrended <- detrended[, sd_check]
    
    # Skip if only NDVI remains
    if (ncol(detrended) < 2) next
    
    # Calculate partial correlations
    pcor_result <- tryCatch(
      {
        pcor(detrended)
      },
      error = function(e) {
        return(NULL)
      }
    )
    
    # Skip if calculation failed
    if (is.null(pcor_result)) next
    
    # Store results in output rasters
    for (var in c("TA", "PPT", "Rad", "CO2", "N", "Crop", "Range")) {
      if (var %in% rownames(pcor_result$estimate)) {
        # Get raster objects dynamically
        est_raster <- get(paste0(tolower(var), "_pcor_raster"))
        pval_raster <- get(paste0("p_", tolower(var), "_pcor_raster"))
        
        # Update values
        est_raster[i, j] <- pcor_result$estimate[var, "NDVI"]
        pval_raster[i, j] <- pcor_result$p.value[var, "NDVI"]
        
        # Save back to environment
        assign(paste0(tolower(var), "_pcor_raster"), est_raster)
        assign(paste0("p_", tolower(var), "_pcor_raster"), pval_raster)
      }
    }
  }
}

# Save all output rasters
variables <- c("ta", "ppt", "rad", "co2", "n", "crop", "range")

for (var in variables) {
  # Save estimate raster
  writeRaster(
    get(paste0(var, "_pcor_raster")),
    filename = file.path(output_dir, paste0(var, "_pcor_raster.tif")),
    format = "GTiff",
    overwrite = TRUE
  )
  
  # Save p-value raster
  writeRaster(
    get(paste0("p_", var, "_pcor_raster")),
    filename = file.path(output_dir, paste0("p_", var, "_pcor_raster.tif")),
    format = "GTiff",
    overwrite = TRUE
  )
}

print("Analysis complete! All results saved.")

