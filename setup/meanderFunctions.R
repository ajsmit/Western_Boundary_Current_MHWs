# This script houses the functions used in "correlate/correlate_meanders.R"


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(RcppRoll)
library(data.table)
library(doMC); doMC::registerDoMC(cores = 50)


# Files -------------------------------------------------------------------

# Bounding boxes
source("setup/regionDefinition.R")

# The AVISO NetCDF files
AVISO_files <- dir(path = "../data/AVISO", pattern = "CMEMS", full.names = T)

# The MHW result files
MHW_files <- dir(path = "../data/MHW", pattern = "MHW.calc.", full.names = T)

# The MHW region files
MHW_region_files <- dir(path = "../data/WBC", pattern = "MHW_sub", full.names = T)

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
  MHW_clim <- MHW_res_sub %>% 
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
masks <- function(AVISO, MHW){
  
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


# Calculate cooccurrence and correlation ----------------------------------

# testers...
# region <- "EAC"
# df <- slice(MHW_sub, 1)

# Wrapper function to run on MHW results
MHW_clim_sub <- function(df, mask_sub){
  res <- df %>%
  filter(paste(lon, lat) %in% paste(mask_sub$lon,
                                    mask_sub$lat)) %>% 
  filter(event == TRUE) %>%
  mutate(intensity = temp-thresh)
  return(res)
}

# Wrapper function for screening out days below the 90th MKE
screen_90 <- function(df){
  res <- df %>% 
    group_by(lon, lat) %>%
    mutate(mke_days = n()) %>%
    filter(mke >= masks$mke_90$mke_90[1]) %>%
    mutate(mke_90_days = n()) %>% 
    filter(event == TRUE) %>%
    mutate(MHW_days = n())
  return(res)
}

# Wrapper function for finding co-occurrence rates
cooc_90 <- function(df){
  res <- df %>% 
    group_by(lon, lat, mke_days, mke_90_days, MHW_days) %>% 
    summarise(cooc_flat = MHW_days[1]/mke_days[1],
              cooc_90 = MHW_days[1]/mke_90_days[1])
  return(res)
}

# Run all the calculations etc.
meander_co_calc <- function(region){
  
  # Load AVISO, MHW, and masks
  AVISO_KE <- fread(paste0("../data/WBC/AVISO_KE_",region,".csv"))
  MHW_clim <- fread(paste0("../data/WBC/MHW_clim_",region,".csv"))
  masks <- readRDS(paste0("masks/masks_",region,".Rds"))

  # First grab MHW as this will be left_join() to AVISO data
  MHW_50 <- MHW_clim_sub(MHW_clim, masks$mke_50)
  MHW_max <- MHW_clim_sub(MHW_clim, masks$max_90)
  rm(MHW_clim); gc()
  
  # Then create the base for the calculations
  AVISO_50 <- AVISO_KE %>% 
    filter(paste(lon, lat) %in% paste(masks$mke_50$lon, 
                                      masks$mke_50$lat)) %>% 
    left_join(MHW_50, by = c("lon", "lat", "t"))
  AVISO_max <- AVISO_KE %>% 
    filter(paste(lon, lat) %in% paste(masks$max_90$lon, 
                                      masks$max_90$lat)) %>% 
    left_join(MHW_max, by = c("lon", "lat", "t"))
  rm(AVISO_KE); gc()
  
  # Prep the data for further calculations
  # Screening out MKE below 90th perc. and days with no MHWs
  calc_50_base <- screen_90(AVISO_50)
  calc_max_base <- screen_90(AVISO_max)
  
  # Calculate the proportion of days when MKE is in the 
  # 90th percentile and MHWs are occurring
  cooc_50 <- cooc_90(calc_50_base)
  cooc_max <- cooc_90(calc_max_base)
  
  # Calculate correlations
  corr_50 <- calc_50_base %>% 
    group_by(lon, lat) %>% 
    summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))
  corr_max <- calc_max_base %>% 
    group_by(lon, lat) %>% 
    summarise(mke_intensity_r = cor(intensity, mke, use = "pairwise.complete.obs"))

  # Merge, save, and clean up
  meander_res <- list(cooc_50 = cooc_50,
                      cooc_max = cooc_max,
                      corr_50 = corr_50,
                      corr_max = corr_max)
  saveRDS(meander_res, file = paste0("correlate/meander_res_",region,".Rds"))
  rm(AVISO_50, AVISO_max, calc_50_base, calc_max_base, cooc_50, cooc_max, 
     corr_50, corr_max,meander_res); gc()
}


# Visualise results -------------------------------------------------------

# testers...
# region <- "GS"

meander_vis <- function(region){
  
  # Determine coordinates and figure dimensions
  coords <- bbox[colnames(bbox) == region][1:4,]
  if(coords[3] > 180) coords[3] <- coords[3]-360
  if(coords[4] > 180) coords[4] <- coords[4]-360
  fig_height <- length(seq(coords[1], coords[2], 1))/3
  fig_width <- length(seq(coords[3], coords[4], 1))/5
  
  # Prep the data
  meander_res <- readRDS(paste0("correlate/meander_res_",region,".Rds")) 
  
  meander_50 <- left_join(meander_res$cooc)%>% 
    ungroup() %>% 
    select(lon, lat, cooc_flat, cooc_90) %>% 
    gather(key = metric, value = val, -lon, -lat) %>% 
    mutate(lon = ifelse(lon > 180, lon-360, lon),
           lat = ifelse(lat > 180, lat-360, lat))
  
  # Create the figure and save
  ggplot(meander_res, aes(x = lon, y = lat)) +
    geom_tile(aes(fill = val), alpha = 0.7) +
    borders(fill = "grey80", colour = "black") +
    coord_equal(xlim = c(coords[3:4]), ylim = c(coords[1:2])) +
    labs(x = NULL, y = NULL) +
    ggtitle(region) +
    scale_fill_viridis_c("co-occurrence\nproportion") +
    facet_wrap(~metric, ncol = 1) +
    theme(legend.position = "bottom")
  ggsave(filename = paste0("figures/",region,"_cooccurrence.pdf"), 
         width = fig_width, height = fig_height)
}

