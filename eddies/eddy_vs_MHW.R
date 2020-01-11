# eddies/eddy_vs_MHW
# The purpose of this script is to quantify how often eddies orginiating in
# the 90th perc. MKE area of WBCs are also in pixels that are flagged as MHWs.
# The hypothesis is that eddies orihinating in areas of higher MKE will be warmer
# than the waters surrounding the WBCs, as the eddies escape the central thrust
# of the WBC and enter the adjoining waters these higher temperatures should then
# be detected as MHWs as the eddy will not be representative of the climatology of the pixel

# This is done with the following steps:
  # 1: Load data
    # A: Eddy trajectories
    # B: MHW results by pixel
    # C: MKE masks, 90th perc.
  # 2: Select eddies that originate in the MKE 90th perc.
  # 3: Merge eddy track and MHW pixels by t/lon/lat
  # 4: Find rate of eddy to MHW co-occurrence


# Setup -------------------------------------------------------------------

# NB: The tikoraluk server does not allow users to install packages
# so first one must activate a user friendly library destination
.libPaths(c("~/R-packages", .libPaths()))

source("setup/meanderFunctions.R")

require(rgeos); require(maptools) # maptools must be loaded after rgeos
library(PBSmapping)
library(ggpubr)


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
  print("Finding eddies within 90th. perc MKE")
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
  pixel_pts_in$MKE_90 <- TRUE

  tracks <- dplyr::right_join(start_eddy, start_pts_in,
                              by = c("lat_eddy", "lon_eddy")) %>%
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
    mutate(MKE_90 = ifelse(is.na(MKE_90), FALSE, MKE_90))

  # Calculate the proportion of pixels experiencing a MHW per day
  eddy_MHW_daily <- eddy_MHW %>%
    group_by(lon_eddy, lat_eddy, t, cyclonic_type, track) %>%
    summarise(pixel_daily = n(),
              pixel_MHW = length(which(event_no > 0)),
              pixel_daily_MHW_prop = round(pixel_MHW/pixel_daily, 4),
              pixel_daily_MHW_prop = replace_na(pixel_daily_MHW_prop, 0),
              # Calculate proportion occuring within MKE 90 perc. region
              pixel_MKE = length(which(MKE_90 == TRUE)),
              pixel_MKE_MHW = length(which(MKE_90 == TRUE & event_no > 0)),
              pixel_MKE_MHW_prop = round(pixel_MKE_MHW/pixel_MKE, 4),
              pixel_MKE_MHW_prop = replace_na(pixel_MKE_MHW_prop, 0),
              # Calculate proportion outside of MKE 90th perc.
              pixel_no_MKE = length(which(MKE_90 == FALSE)),
              pixel_no_MKE_MHW = length(which(MKE_90 == FALSE & event_no > 0)),
              pixel_no_MKE_MHW_prop = round(pixel_no_MKE_MHW/pixel_no_MKE, 4),
              pixel_no_MKE_MHW_prop = replace_na(pixel_no_MKE_MHW_prop, 0),) %>%
    ungroup()

  # Calculate the overall proportions of MHW co-occurrence for each eddy
  eddy_MHW_total <- eddy_MHW_daily %>%
    group_by(track, cyclonic_type) %>%
    summarise_if(is.numeric, mean) %>%
    ungroup() %>%
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

# Map showing all of the eddy tracks
# region <- "BC"
eddy_track_map <- function(region) {

  # Load data
  load(paste0("masks/", region, "-mask_polys.RData"))
  load(paste0("eddies/eddy_vs_MHW_",region,".Rdata"))

  # Prep
  mke <- fortify(mask.list$mke) %>%
    mutate(long = ifelse(long > 180, long - 360, long))
  eddy_MHW_res_daily <- eddy_MHW_res$daily %>%
    ungroup() %>%
    mutate(lon_eddy = ifelse(lon_eddy > 180, lon_eddy - 360, lon_eddy))
  coords <- bbox %>%
    select({{region}}) %>%
    mutate(coord = row.names(.)) %>%
    spread(key = coord, value = {{region}}) %>%
    mutate(lonmin = ifelse(lonmin > 180, lonmin - 360, lonmin),
           lonmax = ifelse(lonmax > 180, lonmax - 360, lonmax))

  # Plot
  ggplot(data = eddy_MHW_res_daily, aes(x = lon_eddy, y = lat_eddy)) +
    geom_point(aes(colour = pixel_no_MKE_MHW_prop), alpha = 0.2) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 colour = "red3", size = 0.4, fill = NA) +
    borders(fill = "grey80") +
    scale_colour_viridis_c(option = "D") +
    coord_cartesian(xlim = c(coords$lonmin, coords$lonmax),
                    ylim = c(coords$latmin, coords$latmax)) +
    labs(x = NULL, y = NULL, colour = "Eddy pixels\nflagged as MHW\n(proportion)")
}

AC_map <- eddy_track_map("AC")
BC_map <- eddy_track_map("BC")
EAC_map <- eddy_track_map("EAC")
GS_map <- eddy_track_map("GS")
KC_map <- eddy_track_map("KC")

eddy_map <- ggarrange(AC_map, BC_map, EAC_map, GS_map, KC_map,
                      ncol = 1, nrow = 5, labels = list("c", "g", "k", "o", "s"))
ggplot2::ggsave("eddies/eddy_map.jpg", plot = eddy_map,
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

# Histogram showing the proportion of MHW days for eddies
# region <- "AC"
eddy_hist_plot <- function(region){
  load(paste0("eddies/eddy_vs_MHW_",region,".Rdata"))
  eddy_MHW_res$total$cyclonic_type <- factor(eddy_MHW_res$total$cyclonic_type,
                                             labels = c("Cyclonic", "Anticyclonic"))
  ggplot(data = eddy_MHW_res$total, aes(x = pixel_daily_MHW_prop)) +
    geom_histogram(aes(fill = cyclonic_type), position = "dodge", binwidth = 0.1) +
    geom_label(aes(label = paste0(region,"\n",
                                  "Cyclonic = ",length(which(eddy_MHW_res$total$cyclonic_type == "Cyclonic")),
                                  "\nAnticyclonic = ",length(which(eddy_MHW_res$total$cyclonic_type == "Anticyclonic"))),
                   x = 0.6, y = nrow(eddy_MHW_res$total)/6)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(expand = F) +
    labs(x = "Eddy pixels flagged as a MHW (proportion)",
         fill = "Eddy type")
}

AC_hist <- eddy_hist_plot("AC")
BC_hist <- eddy_hist_plot("BC")
EAC_hist <- eddy_hist_plot("EAC")
GS_hist <- eddy_hist_plot("GS")
KC_hist <- eddy_hist_plot("KC")

eddy_hist <- ggarrange(AC_hist, BC_hist, EAC_hist, GS_hist, KC_hist,
                       ncol = 1, nrow = 5, labels = list("c", "g", "k", "o", "s"))
ggplot2::ggsave("eddies/eddy_hist.jpg", plot = eddy_hist,
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

