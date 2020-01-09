# eddies/eddy_vs_MHW
# The purpose of this script is to quantify how often eddies orginiating in
# the 50th perc. MKE area of WBCs are also in pixels that are flagged as MHWs.
# The hypothesis is that eddies orihinating in areas of higher MKE will be warmer
# than the waters surrounding the WBCs, as the eddies escape the central thrust
# of the WBC and enter the adjoining waters these higher temperatures should then
# be detected as MHWs as the eddy will not be representative of the climatology of the pixel

# This is done with the following steps:
  # 1: Load data
    # A: Eddy trajectories
    # B: MHW results by pixel
    # C: MKE masks, 50th perc.
  # 2: Select eddies that originate in the MKE 50th perc.
  # 3: Merge eddy track and MHW pixels by t/lon/lat
  # 4: Find rate of eddy to MHW co-occurrence


# Setup -------------------------------------------------------------------

# NB: The tikoraluk server does not allow users to install packages
# so first one must activate a user friendly library destination
.libPaths(c("~/R-packages", .libPaths()))

source("setup/meanderFunctions.R")


# Load data ---------------------------------------------------------------

# region <- "AC"
eddy_vs_MHW <- function(region){
  
  # Load pre-processed detailed eddy data
  # NB: These data are only on tikoraluk
  # They are created in 'correlate/correlate_meanders.R'
  # They use the radius of the eddy to find which pixels in the OISST data would be affected
  print("Loading eddy data")
  AVISO_eddy <- fread(paste0("~/data/WBC/AVISO_eddy_mask_",region,".csv"), nThread = 30) %>% 
    # unite(lon_mask, lat_mask, col = "coord_index", sep = "", remove = F) %>% 
    dplyr::rename(lon_eddy = lon, lat_eddy = lat,
                  lon = lon_mask, lat = lat_mask,
                  t = time) %>%
    mutate(lon = ifelse(lon < 0, lon+360, lon),
           lon_eddy = ifelse(lon_eddy < 0, lon_eddy+360, lon_eddy)) %>% 
    mutate_if(is.numeric, round, 3)
  
  start_eddy <- AVISO_eddy %>% 
    group_by(track) %>% 
    filter(observation_number == 0) %>% 
    ungroup()

  # Load MKE masks
  print("Finding eddies within 50th. perc MKE")
  mke_masks <- readRDS(paste0("masks/masks_",region,".Rds"))
  load(paste0("masks/", region, "-mask_points.Rdata"))
  load(paste0("masks/", region, "-mask_polys.RData"))
  
  # set same proj4string for points as is present in poly
  poly <- mask.list$mke
  start_pts <- SpatialPoints(unique(start_eddy[, c(2,1)]))
  pixel_pts <- SpatialPoints(unique(AVISO_eddy[, c(11,12)]))
  proj4string(start_pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(pixel_pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  # find eddies that start in the WBCs
  start_pts_in <- as.data.frame(start_pts[!is.na(over(start_pts, poly)), ])
  pixel_pts_in <- as.data.frame(pixel_pts[!is.na(over(pixel_pts, poly)), ])
  pixel_pts_in$MKE_50 <- TRUE
  
  tracks <- dplyr::right_join(start_eddy, start_pts_in) %>%
    # dplyr::filter(observation_number == 0) %>%
    dplyr::select(track) %>% 
    unique()
  
  WBC_eddies <- AVISO_eddy %>%
    dplyr::filter(track %in% tracks$track) %>%
    dplyr::mutate(cyclonic_type = as.factor(cyclonic_type),
                  t = as.Date(t))
  
  # Load MHW dresults
  print("Loading MHW clim data")
  MHW_clim <- fread(paste0("~/data/WBC/MHW_clim_",region,".csv"), nThread = 30) %>%
    filter(event == TRUE) %>%
    mutate(intensity = temp-thresh,
           t = as.Date(t)) %>%
    mutate_if(is.numeric, round, 3)
  
  # Join eddy and MHW data
  print("Joining eddy and MHW data")
  eddy_MHW <- left_join(WBC_eddies, MHW_clim, by = c("lon", "lat", "t")) %>% 
    left_join(pixel_pts_in, by = c("lon", "lat")) %>% 
    mutate(MKE_50 = ifelse(is.na(MKE_50), FALSE, MKE_50))
  
  # Calculate the proportion of pixels experiencing a MHW per day
  eddy_MHW_daily <- eddy_MHW %>% 
    group_by(lon_eddy, lat_eddy, t, amplitude, cyclonic_type, observation_number, 
             observed_flag, speed_average, speed_radius, track) %>% 
    # tally()
    mutate(pixel_daily = n(),
           pixel_daily_MHW = length(which(event_no > 0)),
           pixel_daily_MHW_prop = round(pixel_daily_MHW/pixel_daily, 4))

  
  # Save and clean up
  fwrite(AVISO_MHW_eddy_50, file = paste0("~/data/WBC/eddy_vs_MHW_",region,".csv"), nThread = 10)
  rm(AVISO_MHW_eddy_50, AVISO_MHW_eddy_max); gc()
}