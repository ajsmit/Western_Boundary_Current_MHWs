# This script houses the functions used in "correlate/correlate_meanders.R"


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(maps)
# library(furrr)
# library(multidplyr, lib.loc = "../R-packages/")
library(ncdf4)
library(RcppRoll)
library(data.table)
library(ggridges)
doMC::registerDoMC(cores = 50)
# tikoraluk_cluster <- multidplyr::create_cluster(50)
# future::plan(multiprocess)


# Files -------------------------------------------------------------------

# Bounding boxes
source("setup/regionDefinition.R")

# The AVISO NetCDF files
AVISO_files <- dir(path = "../data/AVISO", pattern = "CMEMS", full.names = T)

# The MHW result files
MHW_files <- dir(path = "../data/MHW", pattern = "MHW.calc.", full.names = T)

# The MHW region files
MHW_clim_files <- dir(path = "../data/WBC", pattern = "MHW_clim", full.names = T)

# The OISST lon index
load("../tikoraluk/metadata/lon_OISST.RData")


# Subset AVISO and calculate KE -------------------------------------------

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

# Load AVISO anomaly data and subset accordingly
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
  rm(ugos, vgos, ugosa, vgosa); gc()
  return(res)
}

# Calculate EKE and MKE
# This is designed to run on a dataframe that has already been grouped by lon/lat
AVISO_ke_calc <- function(df) {
  res <- df %>%
    na.omit() %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    # dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))) %>%
    dplyr::filter(eke > -999) %>%
    dplyr::group_by(t) %>% 
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup()
  return(res)
}

# Function for combining and saving the subsetted AVISO KE data
AVISO_KE_save <- function(region){
  
  # Determine the coordinates
  coords <- bbox[colnames(bbox) == region][1:4,]
  
  # Subset the AVISO data
  AVISO_sub <- plyr::ldply(AVISO_files,
                           .fun = AVISO_sub_load, 
                           .parallel = TRUE, 
                           coords = coords)
  
  # Calculate KE and save
  AVISO_KE <- plyr::ddply(AVISO_sub, .variables = c("lon", "lat"),
                          .fun = AVISO_ke_calc, .parallel = TRUE)
  rm(AVISO_sub); gc()
  fwrite(AVISO_KE, file = paste0("../data/WBC/AVISO_KE_",region,".csv"), nThread = 10)
  rm(AVISO_KE); gc()
}


# Subset MHW --------------------------------------------------------------

# testers...
# file_name <- MHW_files_sub[1]

# Wrapper function to load and subset MHW files by longitude slice
MHW_sub_load <- function(file_name, lat_range){
  load(file_name)
  MHW_res_sub <- MHW_res %>% 
    filter(lat %in% lat_range)
  rm(MHW_res); gc()
  return(MHW_res_sub)
}

# Function for combining and saving the subsetted MHW data
MHW_sub_save <- function(region){
  
  # Define lon/lat ranges for extraction from sub files
  coords <- bbox[colnames(bbox) == region][1:4,]
  lat_range <- seq(coords[1]-0.125, coords[2]+0.125, 0.25)
  lon_range <- seq(coords[3]-0.125, coords[4]+0.125, 0.25)
  MHW_files_sub <- MHW_files[which(lon_OISST %in% lon_range)]
  
  # Grab the data and save
  MHW_sub <- plyr::ldply(MHW_files_sub,
                         .fun = MHW_sub_load, 
                         .parallel = TRUE,
                         lat_range = lat_range)
  
  # Create an event data.frame
  MHW_event <- MHW_sub %>% 
    select(-cat) %>% 
    unnest(event) %>% 
    filter(row_number() %% 2 == 0) %>% 
    unnest(event) %>% 
    as.data.table()
  cols <- names(MHW_event)[11:24]
  MHW_event[,(cols) := round(.SD, 3), .SDcols = cols]
  
  # Create a clim data.frame
  MHW_clim <- MHW_sub %>% 
    select(-cat) %>%
    unnest(event) %>% 
    filter(row_number() %% 2 == 1) %>%
    unnest(event) %>%
    select(lon, lat, t, temp, seas, thresh, event, event_no) %>% 
    as.data.table()
  cols <- names(MHW_clim)[5:6]
  MHW_clim[,(cols) := round(.SD, 3), .SDcols = cols]

  # Save
  rm(MHW_sub); gc()
  fwrite(MHW_event, file = paste0("../data/WBC/MHW_event_",region,".csv"), nThread = 10)
  fwrite(MHW_clim, file = paste0("../data/WBC/MHW_clim_",region,".csv"), nThread = 10)
  rm(MHW_event, MHW_clim); gc()
}


# Create MKE percentile masks ---------------------------------------------

# testers...
# region <- "AC"
# AVISO <- AVISO_KE
# MHW <- MHW_event

# Function for creating MKE and max intensity masks
# It expects to be given the AVISO and MHW data as dataframes
mke_masks <- function(AVISO, MHW){
  
  # Find pixels with 90th percentile for max int.
  max_90 <- MHW %>% 
    group_by(lon, lat) %>%
    summarise(intensity_max = max(intensity_max, na.rm = T)) %>% 
    mutate(max_90 = quantile(intensity_max, probs = 0.9)) %>% 
    filter(intensity_max >= max_90)
    
  # Find MKE percentiles
  mke_perc <- AVISO %>%
    group_by(lon, lat) %>%
    summarise(mke = mean(mke, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(mke_90 = quantile(mke, probs = 0.9),
           mke_75 = quantile(mke, probs = 0.75),
           mke_50 = quantile(mke, probs = 0.5))
  
  # Screen out pixels accordingly
  mke_90 <- mke_perc %>%
    filter(mke >= mke_90) %>% 
    select(-mke_75, -mke_50)
  mke_75 <- mke_perc %>% 
    filter(mke >= mke_75, mke < mke_90) %>% 
    select(-mke_75)
  mke_50 <- mke_perc %>% 
    filter(mke >= mke_50, mke < mke_90) %>% 
    select(-mke_50)
  
  # Combine and save
  masks <- list(mke_90 = mke_90,
                mke_75 = mke_75,
                mke_50 = mke_50,
                max_90 = max_90)
  # rm(AVISO_KE)
  saveRDS(masks, file = paste0("masks/masks_",region,".Rds"))
  rm(masks); gc()
  return()
}


# Create eddy trajectory masks --------------------------------------------

# Wrapper function to load a day of AVISO data to get the active grid
# cells inside each bounding box
bbox_cells <- function(region){
  AVISO_sub_load(AVISO_files[1],
                 coords = bbox[colnames(bbox) == region][1:4,]) %>% 
    filter(t == "1993-01-01") %>% 
    select(lon, lat) %>% 
    # Create pixel sides
    mutate(lon = case_when(lon > 180 ~ lon-360,
                           lon <= 180 ~ lon),
           bottom = lat-0.125,
           top = lat+0.125,
           left = lon-0.125,
           right = lon+0.125)
}

# Function that looks at the distance between two points and decides if
# they fall within the radius of the eddy on that day
# testers...
# df <- slice(eddies_sub, 1) %>%
# # unnest() %>%
# select(-lon2, -lat2, -time)
# df <- test %>%
# select(-lon2, -lat2, -time)
dist_calc <- function(df, cells){
  if(nrow(df) != 1) stop("df should be one row exactly")
  # Create a rough, very conservative, radius value for further filtering of pixels
  # before calculating distance of eddies from nearby pixels
  rough_radius <- ceiling(df$speed_radius/50)
  
  # Filter pixels by rough radius
  cells_sub <- cells %>% 
    filter(lon <= df$lon+rough_radius, 
           lon >= df$lon-rough_radius,
           lat <= df$lat+rough_radius,
           lat >= df$lat-rough_radius)
  
  if(nrow(cells_sub) > 0){
    # If there are pixels to mask, find their distances from the eddy centre
    cells_dist <- cells_sub %>% 
      mutate(bottom_left = round(geosphere::distGeo(cbind(.$left, .$bottom), c(df$lon, df$lat))/1000),
             bottom_right = round(geosphere::distGeo(cbind(.$right, .$bottom), c(df$lon, df$lat))/1000),
             top_left = round(geosphere::distGeo(cbind(.$left, .$top), c(df$lon, df$lat))/1000),
             top_right = round(geosphere::distGeo(cbind(.$right, .$top), c(df$lon, df$lat))/1000),
             centre = round(geosphere::distGeo(cbind(.$lon, .$lat), c(df$lon, df$lat))/1000))
    
    # Filter out only pixels with some part covered by the eddy
    cells_filter <- cells_dist %>% 
      filter(bottom_left <= df$speed_radius |
               bottom_right <= df$speed_radius |
               top_left <= df$speed_radius |
               top_right <= df$speed_radius |
               centre <= df$speed_radius) %>% 
      select(lon, lat) %>% 
      dplyr::rename(lon_mask = lon, lat_mask = lat)

    if(nrow(cells_filter) == 0){
      cells_filter <- data.frame(lon_mask = NA, lat_mask = NA)
    }
    # Combine and exit
    # cells_res <- cbind(df, cells_filter)
    # return(cells_res)
    
  } else {
    cells_filter <- data.frame(lon_mask = NA, lat_mask = NA)
  }
  
  # Combine and exit
  cells_res <- cbind(df, cells_filter)
  return(cells_res)
  # return()
}

# This intermediary function helps lighten the load on ddply() in eddy_cells()
# df <- slice(eddies_sub, 1:100) %>% 
  # select(-time)
eddies_inter <- function(df, cells){
  res <- df %>% 
    mutate(lon2 = lon, 
           lat2 = lat) %>% 
    group_by(lon2, lat2) %>% 
    nest() %>% 
    mutate(eddy_res = map(data, dist_calc, cells = cells)) %>% 
    select(eddy_res) %>% 
    unnest()
  return(res)
}

# This function creates the daily index of pixels with eddies in them
# testers...
# region <- "BC"
# eddies_here <- eddies
eddy_cells <- function(region, eddies_here){
  
  cells <- bbox_cells(region)
  
  # The +-4 is to allow for a larger bbox due to the size of the largest eddy
  # We still want the extent of the larger eddies to be considered even when
  # the centre of the eddy is no longer within the bounding box
  eddies_sub <- eddies_here %>% 
    filter(lon <= max(cells$right)+4,
           lon >= min(cells$left)-4,
           lat <= max(cells$top)+4,
           lat >= min(cells$bottom)-4) %>% 
    mutate(lon2 = lon,
           lat2 = lat)
  
  # test <- eddies_sub[103,]
  # test <- eddies_sub[123,]
  
  # Find the pixels masks by each eddy on each day of the dataset
  # doMC::registerDoMC(cores = 50)
  # system.time(
  eddies_res <- plyr::ddply(eddies_sub, .variables = c("time"),
                            # .fun = mean, x = "speed_radius", .parallel = T)
                            .fun = eddies_inter, cells = cells, .parallel = TRUE) %>% 
    na.omit()
  # ) # ~ 363 seconds for EAC
  fwrite(eddies_res, file = paste0("../data/WBC/AVISO_eddy_mask_",region,".csv"), nThread = 10)
  # return(eddies_res)
  rm(eddies_sub, eddies_res); gc()
}

# Function for creating the eddy masks
# This loads the eddy data but then calls eddy_cells() to do the heavy lifting
# testers...
# region <- "BC"
eddy_masks <- function(){
  
  # Load eddy data
  nc <- nc_open("correlate/eddy_trajectory_2.0exp_19930101_20180118.nc")
  eddies <- tibble(lat = round(ncvar_get(nc, varid = "latitude"), 4), # observation longitude
                   lon = round(ncvar_get(nc, varid = "longitude"), 4), # observation latitude
                   # days since 1950-01-01 00:00:00 UTC:
                   time = as.Date(ncvar_get(nc, varid = "time"), origin = "1950-01-01"),
                   # magnitude of the height difference between (Obs) cm the extremum
                   # of SLA within the eddy andthe SLA around the contour defining the
                   # eddy perimeter:
                   amplitude = round(ncvar_get(nc, varid = "amplitude"), 3),
                   # flow orientation of eddies -1 is Cyclonic and 1 is Anticyclonic:
                   cyclonic_type = ncvar_get(nc, varid = "cyclonic_type"),
                   # observation sequence number, days from eddy start:
                   observation_number = ncvar_get(nc, varid = "observation_number"),
                   # flag indicating if the value is interpolated between two observations
                   # or not (0: observed, 1: interpolated):
                   observed_flag = ncvar_get(nc, varid = "observed_flag"),
                   # average speed of the contour defining the radius scale L:
                   speed_average = ncvar_get(nc, varid = "speed_average"),
                   # radius of a circle whose area is equal to that enclosed by the
                   # contour of maximum circum-average speed:
                   speed_radius = ncvar_get(nc, varid = "speed_radius"),
                   # eddy identification number:
                   track = ncvar_get(nc, varid = "track"))
  nc_close(nc); rm(nc)
  
  # correct lon 
  # NB: Needs to be corrected twice due to the way the tracks were assembled
  eddies <- eddies %>% 
    mutate(lon_cor = case_when(lon > 180 ~ lon-360,
                               lon <= 180 ~ lon),
           lon = case_when(lon_cor > 180 ~ lon_cor-360,
                           lon_cor <= 180 ~ lon_cor)) %>% 
    select(-lon_cor)
  # max(eddies$lon)
  # max(eddies$lon, na.rm = T)
  
  # Calculate the eddy masks for each bbox and save
  print(paste0("Began creating eddy masks for AC at ",Sys.time()))
  eddy_cells("AC", eddies)
  print(paste0("Began creating eddy masks for BC at ",Sys.time()))
  eddy_cells("BC", eddies)
  print(paste0("Began creating eddy masks for EAC at ",Sys.time()))
  eddy_cells("EAC", eddies)
  print(paste0("Began creating eddy masks for GS at ",Sys.time()))
  eddy_cells("GS", eddies)
  print(paste0("Began creating eddy masks for KC at ",Sys.time()))
  eddy_cells("KC", eddies)
  
  # Clean up and exit
  rm(eddies); gc()
}


# Mask AVISO and MHW data -------------------------------------------------

# Wrapper function for screening out pixels within an eddy masks
# This function is designed to be fed a dataframe already grouped by date with plyr::ddply
# df <- AVISO_MHW_50 %>%
# filter(t == "1995-01-01")
# eddy_mask <- eddy_masks_region
# merge_eddy <- function(df, eddy_mask){
#   
#   # Filter only the one date used for the calculation
#   # Also filter out pixels with multiple eddies present
#   # selecting for the largest of the eddies
#   eddy_mask_sub <- filter(eddy_mask, time == as.character(df$t[1])) %>% 
#     unite(lon_mask, lat_mask, col = "coord_index", sep = "", remove = F) %>% 
#     dplyr::rename(lon_eddy = lon, lat_eddy = lat) %>% 
#     group_by(coord_index) %>% 
#     filter(speed_radius == max(speed_radius, na.rm = T))
#   
#   # Combine and exit
#   res <- df %>% 
#     left_join(eddy_mask_sub, by = "coord_index") %>% 
#     mutate_if(is.numeric, round, 3)
#   return(res)
# }

# Create the masked dataframes
# region <- "BC"
mask_region <- function(region){
  
  # Load AVISO, MHW, and MKE masks
  AVISO_KE <- fread(paste0("~/data/WBC/AVISO_KE_",region,".csv"), nThread = 30) %>%
    mutate_if(is.numeric, round, 3)
  AVISO_eddy <- fread(paste0("~/data/WBC/AVISO_eddy_mask_",region,".csv"), nThread = 30) %>% 
    # unite(lon_mask, lat_mask, col = "coord_index", sep = "", remove = F) %>% 
    dplyr::rename(lon_eddy = lon, lat_eddy = lat,
                  lon = lon_mask, lat = lat_mask,
                  t = time) %>%
    mutate(lon = ifelse(lon < 0, lon+360, lon),
           lon_eddy = ifelse(lon_eddy < 0, lon_eddy+360, lon_eddy)) %>% 
    mutate_if(is.numeric, round, 3)
  MHW_clim <- fread(paste0("~/data/WBC/MHW_clim_",region,".csv"), nThread = 30) %>%
    filter(event == TRUE) %>%
    mutate(intensity = temp-thresh) %>%
    mutate_if(is.numeric, round, 3)
  mke_masks <- readRDS(paste0("masks/masks_",region,".Rds"))
  
  # First mask the MHW data as these will be left_join() to AVISO data
  print("Filtering MHW data")
  MHW_50 <- left_join(mke_masks$mke_50[,c("lon", "lat")], MHW_clim, by = c("lon", "lat"))
  MHW_max <- left_join(mke_masks$max_90[,c("lon", "lat")], MHW_clim, by = c("lon", "lat"))
  rm(MHW_clim); gc()
  
  # Then mask the AVISO data
  print("Filtering AVISO data")
  AVISO_50 <- left_join(mke_masks$mke_50[,c("lon", "lat")], AVISO_KE, by = c("lon", "lat"))
  AVISO_max <- left_join(mke_masks$max_90[,c("lon", "lat")], AVISO_KE, by = c("lon", "lat"))
  rm(AVISO_KE); gc()
  
  # Join MHW and AVISO data
  print("Joining AVISO and MHW data")
  AVISO_MHW_50 <- left_join(AVISO_50, MHW_50, by = c("lon", "lat", "t"))
  AVISO_MHW_max <- left_join(AVISO_max, MHW_max, by = c("lon", "lat", "t"))
  rm(AVISO_50, MHW_50, AVISO_max, MHW_max); gc()
  
  # Then mask the eddy data
  print("Joining eddy data")
  AVISO_MHW_eddy_50 <- left_join(AVISO_MHW_50, AVISO_eddy, by = c("lon", "lat", "t")) #%>% 
    # group_by(lon, lat, t) %>% 
    # NB: Filter out pixels with multiple eddies present
      # selecting for the largest of the eddies
  # Not currently doing this as it is too slow and not terribly important
    # filter(speed_radius == max(speed_radius, na.rm = T)) %>% 
    # ungroup()
  AVISO_MHW_eddy_max <- left_join(AVISO_MHW_max, AVISO_eddy, by = c("lon", "lat", "t")) #%>% 
    # group_by(lon, lat, t) %>% 
    # filter(speed_radius == max(speed_radius, na.rm = T)) %>% 
    # ungroup()
  rm(AVISO_MHW_50, AVISO_MHW_max); gc()
  
  # Save and clean up
  fwrite(AVISO_MHW_eddy_50, file = paste0("~/data/WBC/AVISO_MHW_eddy_50_",region,".csv"), nThread = 10)
  fwrite(AVISO_MHW_eddy_max, file = paste0("~/data/WBC/AVISO_MHW_eddy_max_",region,".csv"), nThread = 10)
  rm(AVISO_MHW_eddy_50, AVISO_MHW_eddy_max); gc()
}


# Calculate cooccurrence and correlation ----------------------------------

# Wrapper function for screening out multiple things
# testers...
# df <- eddy_merge_50
# mke_val <- mke_masks_region$mke_90$mke_90[1]
count_mke_eddy <- function(df, mke_val){
  
  # Count base days + days with MHWs
  days_base <- df %>% 
    # NB: This is where month grouping would be added to investigate that pattern
    # An earlier version of the code showed no clear monthly pattern
    # So I've stopped calculating it here
    # mutate(month = as.character(lubridate::month(t, label = T))) %>%
    group_by(lon, lat) %>%
    mutate(mke_base_days = n()) %>% 
    group_by(lon, lat, mke_base_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_base_days = n())
  
  # Count days at/above 90th perc. MKE + days with MHWs
  days_90 <- df %>% 
    filter(mke >= mke_val) %>%
    group_by(lon, lat) %>%
    mutate(mke_90_days = n()) %>% 
    group_by(lon, lat, mke_90_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_90_days = n())
  
  # Count days with eddies + days with MHWs
  days_eddy <- df %>% 
    filter(!is.na(amplitude)) %>%
    group_by(lon, lat) %>%
    mutate(eddy_days = n()) %>% 
    group_by(lon, lat, eddy_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_eddy_days = n())
  
  # Count days with cyclonic eddies + days with MHWs
  days_cyclone <- df %>% 
    filter(cyclonic_type == -1) %>%
    group_by(lon, lat) %>%
    mutate(cyclone_days = n()) %>% 
    group_by(lon, lat, cyclone_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_cyclone_days = n())
  
  # Count days with anticyclonic eddies + days with MHWs
  days_anticyclone <- df %>% 
    filter(cyclonic_type == 1) %>%
    group_by(lon, lat) %>%
    mutate(anticyclone_days = n()) %>% 
    group_by(lon, lat, anticyclone_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_anticyclone_days = n())
  
  # Count days without eddies + days with MHWs
  days_no_eddy <- df %>% 
    filter(is.na(amplitude)) %>%
    group_by(lon, lat) %>%
    mutate(no_eddy_days = n()) %>% 
    group_by(lon, lat, no_eddy_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_no_eddy_days = n())
  
  # Count days without eddies and in the 90th perc. mke + days with MHWs
  days_no_eddy_90 <- df %>% 
    filter(is.na(amplitude),
           mke >= mke_val) %>%
    group_by(lon, lat) %>%
    mutate(no_eddy_90_days = n()) %>% 
    group_by(lon, lat, no_eddy_90_days) %>% 
    filter(event == TRUE) %>%
    summarise(MHW_no_eddy_90_days = n())
  
  # Combine and exit
  res <- left_join(days_base, days_90, by = c("lon", "lat")) %>% 
    left_join(days_eddy, by = c("lon", "lat")) %>% 
    left_join(days_cyclone, by = c("lon", "lat")) %>% 
    left_join(days_anticyclone, by = c("lon", "lat")) %>% 
    left_join(days_no_eddy, by = c("lon", "lat")) %>% 
    left_join(days_no_eddy_90, by = c("lon", "lat"))
  return(res)
}

# Wrapper function for finding co-occurrence rates
# df <- count_50
cooc_90 <- function(df){
  
  # Calculate co-occurrence by month
  # res_month <- df %>% 
  #   mutate(month = as.character(lubridate::month(t, label = T))) %>% 
  #   group_by(lon, lat, month, mke_base_days, MHW_base_days, mke_90_days, MHW_90_days) %>% 
  #   summarise(cooc_base = MHW_base_days[1]/mke_base_days[1],
  #             cooc_90 = MHW_90_days[1]/mke_90_days[1],
  #             diff_MHW = MHW_90_days[1]/MHW_base_days[1])
  
  # Calculate total co-occurrence and proportions of occurrence
  res_total <- df %>% 
    # mutate(month = "total") %>% 
    group_by(lon, lat) %>% 
    # Co-occurrence values
    summarise(cooc_base = MHW_base_days/mke_base_days,
              cooc_90 = MHW_90_days/mke_90_days,
              cooc_eddy = MHW_eddy_days/eddy_days,
              cooc_cyclone = MHW_cyclone_days/cyclone_days,
              cooc_anticyclone = MHW_anticyclone_days/anticyclone_days,
              cooc_no_eddy = MHW_no_eddy_days/no_eddy_days,
              cooc_no_eddy_90 = MHW_no_eddy_90_days/no_eddy_90_days,
              # Proportion values
              prop_90 = mke_90_days/mke_base_days,
              prop_90_MHW = MHW_90_days/MHW_base_days,
              prop_eddy = eddy_days/mke_base_days,
              prop_eddy_MHW = MHW_eddy_days/MHW_base_days,
              prop_cyclone = cyclone_days/mke_base_days,
              prop_cyclone_MHW = MHW_cyclone_days/MHW_base_days,
              prop_anticyclone = anticyclone_days/mke_base_days,
              prop_anticyclone_MHW = MHW_anticyclone_days/MHW_base_days,
              prop_cyclone_anticyclone = cyclone_days/anticyclone_days,
              prop_no_eddy = no_eddy_days/mke_base_days,
              prop_no_eddy_MHW = MHW_no_eddy_days/MHW_base_days,
              prop_no_eddy_90 = no_eddy_90_days/mke_base_days,
              prop_no_eddy_90_MHW = MHW_no_eddy_90_days/MHW_base_days)
  
  # Combine and exit
  return(res_total)
  # res <- rbind(res_month, res_total)
  # return(res)
}

# Wrapper function for calculating correlations between metrics
corr_calc <- function(df){
  
  # Calculate correlations by month
  # res_month <- df %>% 
  #   mutate(month = as.character(lubridate::month(t, label = T))) %>% 
  #   group_by(lon, lat, month) %>% 
  #   summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))
  
  # Calculate total correlations
  res_int_mke <- df %>% 
    # mutate(month = "total") %>% 
    group_by(lon, lat) %>%
    summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))
  
  res_int_rad <- df %>% 
    group_by(lon, lat) %>%
    summarise(rad_intensity_r = cor(intensity, speed_radius, use = "pairwise.complete.obs"))
  
  # Combine and exit
  res <- left_join(res_int_mke, res_int_rad, by = c("lon", "lat"))
  return(res)
}

# Run all the calculations etc.
# testers...
# region <- "BC"
meander_co_calc <- function(region){
  
  # Load masked AVISO and MHW data as well as masks
  AVISO_MHW_eddy_50 <- fread(paste0("../data/WBC/AVISO_MHW_eddy_50_",region,".csv"), nThread = 30)
  AVISO_MHW_eddy_max <- fread(paste0("../data/WBC/AVISO_MHW_eddy_max_",region,".csv"), nThread = 30)
  mke_masks_region <- readRDS(paste0("masks/masks_",region,".Rds"))
  
  # Prep the data for further calculations
  # Screen out days when MKE below 90th perc. and days with no MHWs
  # Also screen days with/out eddies present
  print("Screening by MKE and eddies")
  count_50 <- count_mke_eddy(AVISO_MHW_eddy_50, mke_masks_region$mke_90$mke_90[1])
  count_max <- count_mke_eddy(AVISO_MHW_eddy_max, mke_masks_region$mke_90$mke_90[1])
  
  # Calculate rates of co-occurrence between types of days and MHWs
  print("Calculating co-occurrence")
  cooc_50 <- cooc_90(count_50)
  cooc_max <- cooc_90(count_max)
  
  # Calculate correlations
  print("Calculating correlations")
  corr_50 <- corr_calc(AVISO_MHW_eddy_50)
  corr_max <- corr_calc(AVISO_MHW_eddy_max)

  # Merge, save, and clean up
  meander_res <- list(cooc_50 = cooc_50,
                      cooc_max = cooc_max,
                      corr_50 = corr_50,
                      corr_max = corr_max)
  saveRDS(meander_res, file = paste0("correlate/meander_res_",region,".Rds"))
  rm(AVISO_MHW_50, AVISO_MHW_max, mke_masks_region, eddy_masks_region, 
     eddy_merge_50, eddy_merge_max, count_50, count_max, cooc_50, cooc_max, 
     corr_50, corr_max, meander_res); gc()
}


# Visualise results -------------------------------------------------------

# Quick list capable prep function for plotting consistency
list_prep <- function(df){
  # if("cooc_flat" %in% colnames(df)){
    res <- df %>% 
      # ungroup() %>% 
      # select(lon, lat) %>% 
      gather(key = metric, value = val, -lon, -lat)
  # } else {
  #   res <- df %>% 
  #     dplyr::rename(val = mke_intensity_r) %>% 
  #     mutate(metric = "mke_intensity_r")
  # }
  # Not working for "total"...
  # res$month <- factor(res$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
  #                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "total"))
  # res$month[is.na(res$month)] <- "total"
  res <- res %>% 
    ungroup() %>% 
    mutate(lon = ifelse(lon > 180, lon-360, lon))
  return(res)
}

# Create a single faceted plot
# df <- filter(meander_prep$cooc_50)
# df <- filter(meander_prep$corr_50)
# scale_type <- "cooc"
# scale_type <- "prop"
# scale_type <- "cor"
# coords <- coords
plot_res <- function(df, region, scale_type, coords){
  
  # Screen out desired data
  if(scale_type == "cooc"){
    df_prep <- df %>%
      filter(!grepl("prop", metric))
  } else if(scale_type == "prop") {
    df_prep <- df %>%
      filter(!grepl("cooc", metric))
  } else {
    df_prep <- df
  }
  
  # Create the base figure
  fig_res <- ggplot(df_prep, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = val), na.rm = T) +
    borders(fill = "grey80", colour = "black") +
    coord_equal(xlim = c(coords[3:4]), ylim = c(coords[1:2])) +
    labs(x = NULL, y = NULL) +
    ggtitle(region) +
    # scale_fill_viridis_c(scale_label,  option = "A", na.value = NA) +
    # facet_wrap(~metric, ncol = 1) +
    facet_wrap(~metric) +
    theme(legend.position = c(0.7, 0.15),
          legend.key.width = unit(((coords[2]-coords[1])/30), units = "cm"))
  if(scale_type == "cooc"){
    # fig_scale <- fig_res + scale_fill_viridis_c(scale_label,  option = "A", na.value = NA)
    fig_scale <- fig_res + 
      scale_fill_gradientn(guide = guide_colorbar(title = "Co-ocurrence",
                                                  direction = "horizontal",
                                                  title.position = "top"),
                           colors = viridis::viridis(n = 9, option = "A"), limits = c(0, 1))
  } else if(scale_type == "prop"){
    fig_scale <- fig_res + 
      scale_fill_gradientn(guide = guide_colourbar(title = "Proportion",
                                                   direction = "horizontal",
                                                   title.position = "top"),
                           colors = viridis::viridis(n = 9, option = "B"), limits = c(0, 1))
  }  else {
    fig_scale <- fig_res + 
      scale_fill_gradient2(guide = guide_colourbar(title = "Correlation (r)",
                                                   direction = "horizontal",
                                                   title.position = "top"),low = "blue", high = "red") +
      theme(legend.position = "bottom")
  }
  # fig_scale
  return(fig_scale)
}

# This function visualises the co-occurrence and proportion results
# as ridgeplots to give a more numeric representation of the results
# df <- filter(meander_prep$cooc_50)
# df <- filter(meander_prep$corr_50)
# scale_type <- "cooc"
# scale_type <- "prop"
# scale_type <- "cor"
# region <- "BC"
ridge_res <- function(df, region, scale_type){
  # Screen out desired data
  if(scale_type == "cooc"){
    df_prep <- df %>%
      filter(!grepl("prop", metric))
  } else if(scale_type == "prop") {
    df_prep <- df %>%
      filter(!grepl("cooc", metric))
  } else {
    df_prep <- df
  }
  
  
  ggplot(df_prep, aes(x = val, y = metric)) +
    stat_density_ridges(aes(fill = factor(..quantile..)),
                        geom = "density_ridges_gradient", calc_ecdf = TRUE,
                        quantiles = 4, quantile_lines = TRUE) +
    viridis::scale_fill_viridis(discrete = TRUE, name = "Quartiles", alpha = 0.7) +
    ggtitle(region) +
    labs(y = NULL, x = scale_type) +
    scale_x_continuous(limits = c(0, 1))
  
  # ggplot(df_prep, aes(x = val, y = metric)) +
  #   geom_density_ridges(
  #     jittered_points = TRUE,
  #     position = position_points_jitter(width = 0, height = 0),
  #     point_shape = "|", point_size = 2, point_alpha = 0.7, alpha = 0) +
  #   stat_density_ridges(aes(fill = factor(..quantile..)),
  #                       geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
  #   scale_fill_manual(name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
  #                     labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
  #   ggtitle(region) +
  #   labs(y = NULL, x = scale_type) +
  #   scale_x_continuous(limits = c(0, 1))
    
  # ggplot(df_prep, aes(x = val, y = metric)) +
  #   ggridges::geom_density_ridges() +
  #   ggtitle(region) +
  #   labs(y = NULL, x = scale_type) +
  #   scale_x_continuous(limits = c(0, 1))
    # facet_wrap(~metric)
}


# This function creates all of the current visual outputs
# region <- "BC"
meander_vis <- function(region){
  
  # Determine coordinates and figure dimensions
  coords <- bbox[colnames(bbox) == region][1:4,]
  if(coords[3] > 180) coords[3] <- coords[3]-360
  if(coords[4] > 180) coords[4] <- coords[4]-360
  fig_height <- length(seq(coords[1], coords[2], 1))/3
  fig_width <- length(seq(coords[3], coords[4], 1))/4
  
  # Prep the data
  meander_res <- readRDS(paste0("correlate/meander_res_",region,".Rds"))
  meander_prep <- lapply(meander_res, list_prep)
  
  # Create the 50th perc. MKE co-occurrence and proportion figures
  cooc_50 <- plot_res(meander_prep$cooc_50, region, "cooc", coords)
  ggsave(cooc_50, filename = paste0("figures/",region,"_cooc_50.pdf"), 
         width = fig_width, height = fig_height)
  cooc_50_ridge <- ridge_res(meander_prep$cooc_50, region, "cooc")
  ggsave(cooc_50_ridge, filename = paste0("figures/",region,"_cooc_50_ridge.pdf"), 
         width = 6, height = 6)
  prop_50 <- plot_res(meander_prep$cooc_50, region, "prop", coords)
  ggsave(prop_50, filename = paste0("figures/",region,"_prop_50.pdf"), 
         width = fig_width, height = fig_height)
  prop_50_ridge <- ridge_res(meander_prep$cooc_50, region, "prop")
  ggsave(prop_50_ridge, filename = paste0("figures/",region,"_prop_50_ridge.pdf"), 
         width = 6, height = 6)
  
  # Create the total 90th perc. max int. co-occurrence figure
  cooc_max <- plot_res(meander_prep$cooc_max, region, "cooc", coords)
  ggsave(cooc_max, filename = paste0("figures/",region,"_cooc_max.pdf"), 
         width = fig_width, height = fig_height)
  cooc_max_ridge <- ridge_res(meander_prep$cooc_max, region, "cooc")
  ggsave(cooc_max_ridge, filename = paste0("figures/",region,"_cooc_max_ridge.pdf"), 
         width = 6, height = 6)
  prop_max <- plot_res(meander_prep$cooc_max, region, "prop", coords)
  ggsave(prop_max, filename = paste0("figures/",region,"_prop_max.pdf"), 
         width = fig_width, height = fig_height)
  prop_max_ridge <- ridge_res(meander_prep$cooc_max, region, "prop")
  ggsave(prop_max_ridge, filename = paste0("figures/",region,"_prop_max_ridge.pdf"), 
         width = 6, height = 6)
  
  # NB: There is no clear correlation pattern so I am not running the code below
  # They will run if they are uncommented
  # Create the total 50th perc. MKE correltaionfigure
  # corr_50_total <- plot_res(filter(meander_prep$corr_50, month == "total"), F,
  #                           "correlation (r)", "corr", coords)
  # corr_50_total
  # ggsave(corr_50_total, filename = paste0("figures/",region,"_corr_50_total.pdf"), 
  #        width = fig_width, height = fig_height)
  
  # Create the total 90th perc. max int. correlation figure
  # corr_max_total <- plot_res(filter(meander_prep$corr_max, month == "total"), F,
  #                            "correlation (r)", "corr", coords)
  # corr_max_total
  # ggsave(corr_max_total, filename = paste0("figures/",region,"_corr_max_total.pdf"), 
  #        width = fig_width, height = fig_height)
  
  # NB: There appears to be little in the way of any monthly trend so 
  # I am not running the monthly figures here
  # Simply un-comment them and they will run
  # Create the monthly 50th perc. MKE  co-occurence figure
  # cooc_50_month <- plot_res(filter(meander_prep$cooc_50, month != "total"), T,
  #                           "co-occurrence\nproportion", "A")
  # cooc_50_month
  # ggsave(cooc_50_month, filename = paste0("figures/",region,"_cooc_50_month.pdf"), 
  #        width = fig_width*1.5, height = fig_height*2)
  
  # Create the monthly 50th perc. MKE  co-occurence figure
  # cooc_max_month <- plot_res(filter(meander_prep$cooc_max, month != "total"), T,
  #                            "co-occurrence\nproportion", "A")
  # cooc_max_month
  # ggsave(cooc_max_month, filename = paste0("figures/",region,"_cooc_max_month.pdf"), 
  #        width = fig_width*1.5, height = fig_height*2)
}

