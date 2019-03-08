# This script houses the functions used in "correlate/correlate_meanders.R"


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(doMC); doMC::registerDoMC(cores = 50)


# Files -------------------------------------------------------------------

# Bounding boxes
source("setup/regionDefinition.R")

# The AVISO NetCDF files
AVISO_files <- dir(path = "../data/AVISO", pattern = "CMEMS", full.names = T)

# The MHW result files
MHW_files <- dir(path = "../data/MHW", pattern = "MHW.calc.", full.names = T)

# The OISST lon index
load("../tikoraluk/metadata/lon_OISST.RData")


# Subset AVISO ------------------------------------------------------------

# testers...
# file_name <- AVISO_files[1]
# region <- "AC"
# coords <- bbox$AC
# var_id <- "ugos"

# Ease function for extracting AVISO variables
AVISO_var <- function(file_name, var_id, coords){
  # Define lon/lat ranges for extraction
  lat_range <- seq(coords[1]-0.125, coords[2]+0.125, 0.25)
  lon_range <- seq(coords[3]-0.125, coords[4]+0.125, 0.25)
  nc <- nc_open(as.character(file_name))
  lat_sub <- which(nc$dim$latitude$vals %in% lat_range)
  lon_sub <- which(nc$dim$longitude$vals %in% lon_range)
  # Subset data
  res <- ncvar_get(nc, varid = var_id, 
                   start = c(lon_sub[1], lat_sub[1], 1),
                   count = c(length(lon_range), length(lat_range), -1))
  dimnames(res) <- list(lon = lon_range,
                        lat = lat_range,
                        t = nc$dim$time$vals)
  nc_close(nc)
  # Melt and finish
  res <- as.data.frame(reshape2::melt(res, value.name = var_id), row.names = NULL) %>% 
    mutate(t = as.Date(t, origin = "1950-01-01"))
  return(res)
}

# Load AVSIO anomaly data and subset accordingly
AVISO_sub_load <- function(file_name, coords){
  # Extract the four desired variables
  ugos <- AVISO_var(file_name, "ugos", coords)
  vgos <- AVISO_var(file_name, "vgos", coords)
  ugosa <- AVISO_var(file_name, "ugosa", coords)
  vgosa <- AVISO_var(file_name, "vgosa", coords)
  # Steek'im baiby
  res <- left_join(ugos, vgos, by = c("lon", "lat", "t")) %>% 
    left_join(ugosa, by = c("lon", "lat", "t")) %>% 
    left_join(vgosa, by = c("lon", "lat", "t"))
  return(res)
}

# Function for combining and saving the subsetted AVISO data
AVISO_sub_save <- function(region){
  coords <- bbox[colnames(bbox) == region][1:4,]
  AVISO_sub <- plyr::ldply(AVISO_files,
                                .fun = AVISO_sub_load, 
                                .parallel = TRUE, 
                                coords = coords)
  save(AVISO_sub, file = paste0("../data/WBC/AVISO_sub_",region,".Rdata"))
}


# Subset MHW --------------------------------------------------------------

# testers...
# file_name <- MHW_files_sub[1]

# Load AVSIO anomaly data and subset accordingly
MHW_sub_load <- function(file_name, lat_range){
  load(file_name)
  MHW_res_sub <- MHW_res %>% 
    filter(lat %in% lat_range)
  rm(MHW_res)
  return(MHW_res_sub)
}

# Function for combining and saving the subsetted AVISO data
MHW_sub_save <- function(region){
  # Define lon/lat ranges for extraction from sub files
  coords <- bbox[colnames(bbox) == region][1:4,]
  lat_range <- seq(coords[1]-0.125, coords[2]+0.125, 0.25)
  lon_range <- seq(coords[3]-0.125, coords[4]+0.125, 0.25)
  MHW_files_sub <- MHW_files[which(lon_OISST %in% lon_range)]
  # Grab the data
  MHW_sub <- plyr::ldply(MHW_files_sub,
                         .fun = MHW_sub_load, 
                         .parallel = TRUE,
                         lat_range = lat_range)
  save(MHW_sub, file = paste0("../data/WBC/MHW_sub_",region,".Rdata"))
}


# Calculate EKE -----------------------------------------------------------

ke.fun <- function(df) {
  out <- kdf %>%
    na.omit() %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))) %>%
    dplyr::filter(eke > -999) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  return(out)
}
