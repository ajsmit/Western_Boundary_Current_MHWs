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

require(rgeos); require(maptools) # maptools must be loaded after rgeos
library(PBSmapping)


# Function ----------------------------------------------------------------

# region <- "AC"
eddy_vs_MHW <- function(region) {

  # Load pre-processed detailed eddy data
  # NB: These data are only on tikoraluk
  # They are created in 'correlate/correlate_meanders.R'
  # They use the radius of the eddy to find which pixels in the OISST data would be affected
  print("Loading eddy data")
  AVISO_eddy <- fread(paste0("~/data/WBC/AVISO_eddy_mask_",region,".csv"), nThread = 30) %>%
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

  # Set same proj4string for points as is present in poly
  poly <- mask.list$mke
  start_pts <- SpatialPoints(unique(start_eddy[, c(2,1)]))
  pixel_pts <- SpatialPoints(unique(AVISO_eddy[, c(11,12)]))
  proj4string(start_pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(pixel_pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # Find eddies that start in the WBCs
  start_pts_in <- as.data.frame(start_pts[!is.na(over(start_pts, poly)), ])
  pixel_pts_in <- as.data.frame(pixel_pts[!is.na(over(pixel_pts, poly)), ])
  pixel_pts_in$MKE_50 <- TRUE

  tracks <- dplyr::right_join(start_eddy, start_pts_in) %>%
    dplyr::select(track) %>%
    unique()

  WBC_eddies <- AVISO_eddy %>%
    dplyr::filter(track %in% tracks$track) %>%
    dplyr::mutate(cyclonic_type = as.factor(cyclonic_type),
                  t = as.Date(t))

  # Load MHW results
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
    group_by(lon_eddy, lat_eddy, t, cyclonic_type, track) %>%
    summarise(pixel_daily = n(),
              pixel_MHW = length(which(event_no > 0)),
              pixel_daily_MHW_prop = round(pixel_MHW/pixel_daily, 4),
              pixel_daily_MHW_prop = replace_na(pixel_daily_MHW_prop, 0),
              # Calculate proportion occuring within MKE 50 perc region
              pixel_MKE = length(which(MKE_50 == TRUE)),
              pixel_MKE_MHW = length(which(MKE_50 == TRUE & event_no > 0)),
              pixel_MKE_MHW_prop = round(pixel_MKE_MHW/pixel_MKE, 4),
              pixel_MKE_MHW_prop = replace_na(pixel_MKE_MHW_prop, 0),
              # Calculate proportion outside of MKE 50th perc.
              pixel_no_MKE = length(which(MKE_50 == FALSE)),
              pixel_no_MKE_MHW = length(which(MKE_50 == FALSE & event_no > 0)),
              pixel_no_MKE_MHW_prop = round(pixel_no_MKE_MHW/pixel_no_MKE, 4),
              pixel_no_MKE_MHW_prop = replace_na(pixel_no_MKE_MHW_prop, 0),)

  # Calculate the overall proportions of MHW co-occurrence for each eddy
  eddy_MHW_total <- eddy_MHW_daily %>%
    group_by(track, cyclonic_type) %>%
    summarise_if(is.numeric, mean) %>%
    mutate_if(is.numeric, round, 2)

  # Save and clean up
  eddy_MHW_res <- list(daily = eddy_MHW_daily,
                       total = eddy_MHW_total)
  save(eddy_MHW_res, file = paste0("eddies/eddy_vs_MHW_",region,".Rdata"))
  rm(AVISO_eddy, MHW_clim); gc()
}


# Calculate ---------------------------------------------------------------

# Run sequentially
for(i in 1:(ncol(bbox)-1)) {

  # Determine region
  region <- colnames(bbox)[i]
  print(paste0("Began run on ",region," at ",Sys.time()))

  # Mask the regions
     # NB: This runs the MKE and Eddy masks
  eddy_vs_MHW(region)
  print(paste0("Finished run on ",region," at ",Sys.time()))

  # Clear up some RAM
  gc()
} # ~ 2.5 minutes each


# Visualise ---------------------------------------------------------------


