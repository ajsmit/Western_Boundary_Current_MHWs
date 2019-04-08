# This script houses the functions used in "correlate/correlate_meanders.R"


# Load libraries ----------------------------------------------------------

library(tidyverse)
# library(furrr)
# library(multidplyr, lib.loc = "../R-packages/")
library(ncdf4)
library(RcppRoll)
library(data.table)
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
  
  # Filter pixles by rough radius
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

# testers...
# region <- "EAC"
# df <- AVISO_KE
# mask_sub <- masks$mke_50
# AVISO <- AVISO_KE
# MHW <- MHW_50

# This function is designed to be passed ddply data grouped by c("lon", "lat")
MHW_clim_sub <- function(df, mask_sub){
  
  # Create coordinate indexes
  mask_sub$coord_index <- paste0(mask_sub$lon, mask_sub$lat)
  df$coord_index <- paste0(df$lon, df$lat)
  
  # Filter by event column if the data are within the mask
  if(df$coord_index %in% mask_sub$coord_index){
    res <- df %>% 
      select(-coord_index) %>% 
      filter(event == TRUE) %>%
      mutate(intensity = temp-thresh)
    return(res)
  }
}

# This function is designed to be passed ddply data grouped by c("lon", "lat")
AVISO_KE_sub <- function(AVISO, mask_sub){
  
  # Create coordinate indexes
  mask_sub$coord_index <- paste0(mask_sub$lon, mask_sub$lat)
  AVISO$coord_index <- paste0(AVISO$lon, AVISO$lat)
  
  # Filter and join as required
  if(AVISO$coord_index %in% mask_sub$coord_index){
    res <- AVISO
  return(res)
  }
}

# Create the masked dataframes
mask_region <- function(region){
  
  # Load AVISO, MHW, and masks
  AVISO_KE <- fread(paste0("../data/WBC/AVISO_KE_",region,".csv"), nThread = 10)
  MHW_clim <- fread(paste0("../data/WBC/MHW_clim_",region,".csv"), nThread = 10)
  masks <- readRDS(paste0("masks/masks_",region,".Rds"))
  
  # First mask the MHW data as these will be left_join() to AVISO data
  print("Filtering MHW data")
  MHW_50 <- plyr::ddply(MHW_clim, .variables = c("lon", "lat"), 
                        .fun = MHW_clim_sub, .parallel = T, mask_sub = masks$mke_50)
  MHW_max <- plyr::ddply(MHW_clim, .variables = c("lon", "lat"), 
                         .fun = MHW_clim_sub, .parallel = T, mask_sub = masks$max_90)
  rm(MHW_clim); gc()
  
  # Then mask the AVISO data
  print("Filtering AVISO data")
  AVISO_50 <- plyr::ddply(AVISO_KE, .variables = c("lon", "lat"), 
                          .fun = AVISO_KE_sub, .parallel = T, mask_sub = masks$mke_50)
  AVISO_max <- plyr::ddply(AVISO_KE, .variables = c("lon", "lat"), 
                           .fun = AVISO_KE_sub, .parallel = T, mask_sub = masks$max_90)
  rm(AVISO_KE); gc()
  
  # Join and save
  print("Joining data")
  AVISO_MHW_50 <- left_join(AVISO_50, MHW_50, by = c("lon", "lat", "t"))
  fwrite(AVISO_MHW_50, file = paste0("../data/WBC/AVISO_MHW_50_",region,".csv"), nThread = 10)
  AVISO_MHW_max <- left_join(AVISO_max, MHW_max, by = c("lon", "lat", "t"))
  fwrite(AVISO_MHW_max, file = paste0("../data/WBC/AVISO_MHW_max_",region,".csv"), nThread = 10)
  rm(AVISO_MHW_50, AVISO_MHW_max); gc()
}


# Calculate cooccurrence and correlation ----------------------------------

# testers...
# region <- "EAC"

# Wrapper function for screening out pixels within an eddy masks
# This function is designed to be fed a dataframe already group by date with plyr::ddply
# df <- AVISO_MHW_50 %>%
  # filter(t == "1995-01-01")
# eddy_mask <- eddy_masks_region
merge_eddy <- function(df, eddy_mask){
  
  # Filter only the one date used for the calculation
  # Also filter out pixels with multiple eddies present
    # selecting for the largest of the eddies
  eddy_mask_sub <- filter(eddy_mask, time == as.character(df$t[1])) %>% 
    unite(lon_mask, lat_mask, col = "coord_index", sep = "", remove = F) %>% 
    dplyr::rename(lon_eddy = lon, lat_eddy = lat) %>% 
    group_by(coord_index) %>% 
    filter(speed_radius == max(speed_radius, na.rm = T))
  
  # Combine and exit
  res <- df %>% 
    left_join(eddy_mask_sub, by = "coord_index") #%>%
    # group_by(lat, lon, t) %>% 
    # filter(speed_radius == max(speed_radius, na.rm = T)) %>% 
    # summarise_if(is.numeric, mean, na.rm = T) %>% 
    # summarise()
    # ungroup()
  res[is.na(res)] <- NA
  return(res)
}

# Wrapper function for screening out multiple things
# testers...
# df <- AVISO_MHW_50
# eddy_mask <- eddy_masks_region
screen_mke_eddy <- function(df, eddy_mask, mke_val){
  
  # Add eddy mask information into the AVISO dataframe
  # This is a very large task...
  df_eddy <- plyr::ddply(df, .variables = "t", .fun = merge_eddy,
                         eddy_mask = eddy_mask, .parallel = T)
  
  # Count base days and the occurrence of MHWs
  res_base <- df %>% 
    group_by(lon, lat) %>%
    mutate(mke_base_days = n()) %>%
    filter(event == TRUE) %>%
    # filter(mke >= mke_val) %>%
    mutate(MHW_base_days = n()) #%>% 
    # mutate(MHW_days = n())
  
  # Screen out days below 90th perc. MKE and count days with MHWs
  res_90 <- df %>% 
    group_by(lon, lat) %>%
    # mutate(mke_days = n()) %>%
    filter(mke >= mke_val) %>%
    mutate(mke_90_days = n()) %>% 
    filter(event == TRUE) %>%
    mutate(MHW_90_days = n())
  
  # Combine and exit
  res <- left_join(res_base, res_90, 
                   by = c("lon", "lat", "t", "mke", "eke", "coord_index",
                          "temp", "seas", "thresh", "event", "event_no", "intensity"))
  return(res)
}

# Wrapper function for finding co-occurrence rates
cooc_90 <- function(df){
  
  # Calculate co-occurrence by month
  res_month <- df %>% 
    mutate(month = as.character(lubridate::month(t, label = T))) %>% 
    group_by(lon, lat, month, mke_base_days, MHW_base_days, mke_90_days, MHW_90_days) %>% 
    summarise(cooc_base = MHW_base_days[1]/mke_base_days[1],
              cooc_90 = MHW_90_days[1]/mke_90_days[1],
              diff_MHW = MHW_90_days[1]/MHW_base_days[1])
  
  # Calculate total co-occurrence
  res_total <- df %>% 
    mutate(month = "total") %>% 
    group_by(lon, lat, month, mke_base_days, MHW_base_days, mke_90_days, MHW_90_days) %>% 
    summarise(cooc_base = MHW_base_days[1]/mke_base_days[1],
              cooc_90 = MHW_90_days[1]/mke_90_days[1],
              diff_MHW = MHW_90_days[1]/MHW_base_days[1])
  
  # Combine and exit
  res <- rbind(res_month, res_total)
  return(res)
}

corr_calc <- function(df){
  
  # Calculate crrelations by month
  res_month <- df %>% 
    mutate(month = as.character(lubridate::month(t, label = T))) %>% 
    group_by(lon, lat, month) %>% 
    summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))
  
  # Calculate total correlations
  res_total <- df %>% 
    mutate(month = "total") %>% 
    group_by(lon, lat, month) %>% 
    summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))
  
  # Combine and exit
  res <- rbind(res_month, res_total)
  return(res)
}

# Run all the calculations etc.
# testers...
# region <- "AC"
meander_co_calc <- function(region){
  
  # Load masked AVISO and MHW data as well as masks
  AVISO_MHW_50 <- fread(paste0("../data/WBC/AVISO_MHW_50_",region,".csv"), nThread = 10)
  AVISO_MHW_max <- fread(paste0("../data/WBC/AVISO_MHW_max_",region,".csv"), nThread = 10)
  mke_masks_region <- readRDS(paste0("masks/masks_",region,".Rds"))
  eddy_masks_region <- fread(paste0("~/data/WBC/AVISO_eddy_mask_",region,".csv"), nThread = 10)
  
  # Prep the data for further calculations
  # CScreen out days when MKE below 90th perc. and days with no MHWs
  # Also screen days with/out eddies present
  print("Screening by MKE and eddies")
  screen_50 <- screen_mke_eddy(AVISO_MHW_50, mke_masks_region$mke_90$mke_90[1])
  screen_max <- screen_mke_eddy(AVISO_MHW_max, mke_masks_region$mke_90$mke_90[1])
  
  # Calculate the proportion of days when MKE is in the 
  # 90th percentile and MHWs are occurring
  print("Calculating co-occurrence")
  cooc_50 <- cooc_90(screen_50)
  cooc_max <- cooc_90(screen_max)
  
  # Calculate correlations
  print("Calculating correlations")
  corr_50 <- corr_calc(calc_50_eddy)
  corr_max <- corr_calc(calc_max_eddy)

  # Merge, save, and clean up
  meander_res <- list(cooc_50 = cooc_50,
                      cooc_max = cooc_max,
                      corr_50 = corr_50,
                      corr_max = corr_max)
  saveRDS(meander_res, file = paste0("correlate/meander_res_",region,".Rds"))
  rm(AVISO_MHW_50, AVISO_MHW_max, calc_50_base, calc_max_base, calc_50_eddy, calc_max_eddy,
     cooc_50, cooc_max, corr_50, corr_max,meander_res); gc()
}


# Visualise results -------------------------------------------------------

# testers...
# region <- "GS"
# df <- meander_res$corr_50

# Quick list capable prep function for plotting consistency
list_prep <- function(df){
  if("cooc_flat" %in% colnames(df)){
    res <- df %>% 
      ungroup() %>% 
      select(lon, lat, month, cooc_flat, cooc_90) %>% 
      gather(key = metric, value = val, -lon, -lat, -month)
  } else {
    res <- df %>% 
      dplyr::rename(val = mke_intensity_r) %>% 
      mutate(metric = "mke_intensity_r")
  }
  # Not working for "total"...
  res$month <- factor(res$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "total"))
  res$month[is.na(res$month)] <- "total"
  res <- res %>% 
    ungroup() %>% 
    mutate(lon = ifelse(lon > 180, lon-360, lon),
           lat = ifelse(lat > 180, lat-360, lat))
  return(res)
}

# Create a single faceted plot
plot_res <- function(df, month_facet, scale_label, scale_type, coords){
  fig_res <- ggplot(df, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = val), na.rm = T) +
    borders(fill = "grey80", colour = "black") +
    coord_equal(xlim = c(coords[3:4]), ylim = c(coords[1:2])) +
    labs(x = NULL, y = NULL) +
    ggtitle(region) +
    # scale_fill_viridis_c(scale_label,  option = "A", na.value = NA) +
    # facet_wrap(~metric, ncol = 1) +
    theme(legend.position = "bottom")
  if(scale_type == "cooc"){
    # fig_scale <- fig_res + scale_fill_viridis_c(scale_label,  option = "A", na.value = NA)
    fig_scale <- fig_res + scale_fill_gradientn(colors = viridis::viridis(n = 9, option = "A"), limits = c(0, 1))
  } else {
    fig_scale <- fig_res + scale_fill_gradient2(low = "blue", high = "red")
  }
  if(month_facet){
    fig_res_facet <- fig_scale +
      facet_wrap(month~metric, ncol = 6)
  } else {
    fig_res_facet <- fig_scale +
      facet_wrap(~metric, ncol = 1)
  }
  return(fig_res_facet)
}

# This function creates all of the current visual outputs
meander_vis <- function(region){
  
  # Determine coordinates and figure dimensions
  coords <- bbox[colnames(bbox) == region][1:4,]
  if(coords[3] > 180) coords[3] <- coords[3]-360
  if(coords[4] > 180) coords[4] <- coords[4]-360
  fig_height <- length(seq(coords[1], coords[2], 1))/3
  fig_width <- length(seq(coords[3], coords[4], 1))/5
  
  # Prep the data
  meander_res <- readRDS(paste0("correlate/meander_res_",region,".Rds")) 
  meander_prep <- lapply(meander_res, list_prep)
  
  # Create the total 50th perc. MKE co-occurrence figure
  cooc_50_total <- plot_res(filter(meander_prep$cooc_50, month == "total"), F,
                            "co-occurrence\nproportion", "cooc", coords)
  # cooc_50_total
  ggsave(cooc_50_total, filename = paste0("figures/",region,"_cooc_50_total.pdf"), 
         width = fig_width, height = fig_height)
  
  # Create the total 90th perc. max int. co-occurrence figure
  cooc_max_total <- plot_res(filter(meander_prep$cooc_max, month == "total"), F,
                             "co-occurrence\nproportion", "cooc", coords)
  # cooc_max_total
  ggsave(cooc_max_total, filename = paste0("figures/",region,"_cooc_max_total.pdf"), 
         width = fig_width, height = fig_height)
  
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

