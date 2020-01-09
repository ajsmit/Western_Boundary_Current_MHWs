#

library(ncdf4)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(data.table)
library(colorspace)
library(OneR) # for binning


# Miscellaneous things to load --------------------------------------------

# Read in the region-specific components and functions
source("setup/functions.R")
source("setup/regionDefinition.R")
source("setup/plot.layers.R")


# The data ----------------------------------------------------------------

# eddies
# for converting julian date to calendar date, see:
# https://stackoverflow.com/questions/20789612/convert-julian-date-to-calendar-dates-within-a-data-frame-in-r

ncFile <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/value-added/eddy-trajectory/eddy_trajectory_2.0exp_19930101_20180118.nc"
nc <- nc_open(ncFile)
eddies <- tibble(lat = round(ncvar_get(nc, varid = "latitude"), 1), # observation longitude
                 lon = round(ncvar_get(nc, varid = "longitude"), 1), # observation latitude
                 # days since 1950-01-01 00:00:00 UTC:
                 time = as.Date(ncvar_get(nc, varid = "time"), origin = "1950-01-01"),
                 # magnitude of the height difference between (Obs) cm the extremum
                 # of SLA within the eddy andthe SLA around the contour defining the
                 # eddy perimeter:
                 amplitude = round(ncvar_get(nc, varid = "amplitude"), 1),
                 # flow orientation of eddies -1 is Cyclonic and 1 is Anticyclonic:
                 cyclonic_type = ncvar_get(nc, varid = "cyclonic_type"),
                 # eddy identification number:
                 track = ncvar_get(nc, varid = "track"),
                 # observation sequence number, days from eddy start:
                 observation_number = ncvar_get(nc, varid = "observation_number"),
                 # flag indicating if the value is interpolated between two observations
                 # or not (0: observed, 1: interpolated):
                 observed_flag = ncvar_get(nc, varid = "observed_flag"),
                 # average speed of the contour defining the radius scale L:
                 speed_average = round(ncvar_get(nc, varid = "speed_average"), 1),
                 # radius of a circle whose area is equal to that enclosed by the
                 # contour of maximum circum-average speed:
                 speed_radius = round(ncvar_get(nc, varid = "speed_radius"), 1))
nc_close(nc); rm(nc)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- vroom::vroom(file = paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv"),
                          delim = ",",
                          col_select = list(lon, lat, event_no, duration, date_start),
                          col_types = list(c))

BC_events <- as_tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as_tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as_tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as_tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))


# Bin MHWs to .1 degree pixels --------------------------------------------

AC_events <- AC_events %>%
  dplyr::mutate(lon = round(lon, 1),
                lat = round(lat, 1))




# The function ------------------------------------------------------------

# region <- "AC"
# plot.parameters <- AC.layers
# ev.date <- AC.ev

eddy_plot <- function(eddies, region, plot.parameters, ev.date) {
  # bathy data
  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))

  # some plots of eddies
  # add column woth duration of each eddy

  eddies <- eddies %>%
    dplyr::filter(lat >= bbox[1, region] & lat <= bbox[2, region]) %>%
    dplyr::filter(lon >= bbox[3, region] & lon <= bbox[4, region]) %>%
    dplyr::group_by(track) %>%
    dplyr::mutate(duration = max(observation_number + 1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin = bin(duration, method = "content"))

  # find eddies that start in the WBCs
  load(paste0("masks/", region, "-mask_points.Rdata"))
  load(paste0("masks/", region, "-mask_polys.Rdata"))

  # polys
  poly <- mask.list$mke90

  # set same proj4string for points as is present in poly
  pts <- SpatialPoints(eddies[, c(2,1)])
  proj4string(pts) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # apply over function to extract points that fall within poly
  pts_in <- as.data.frame(pts[!is.na(over(pts, poly)), ])

  tracks <- dplyr::right_join(eddies, pts_in) %>%
    dplyr::filter(observation_number == 0) %>%
    dplyr::select(track)

  WBC_eddies <- eddies %>%
    dplyr::filter(track %in% c(tracks$track)) %>%
    dplyr::mutate(cyclonic_type = as.factor(cyclonic_type))

  WBC_eddies_ev <- WBC_eddies %>%
    dplyr::filter(time >= as.Date(ev.date[1]) & time <= as.Date(ev.date[2]) |
                    time >= as.Date(ev.date[3]) & time <= as.Date(ev.date[4]) |
                    time >= as.Date(ev.date[5]) & time <= as.Date(ev.date[6])) %>%
    dplyr::select(track) %>%
    unique

  ev_eddies <- dplyr::right_join(WBC_eddies, WBC_eddies_ev)

  load(paste0("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/", region, "-mask_polys.RData"))
  mke <- fortify(mask.list$mke90)

  eddy_plt <- ggplot(data = WBC_eddies, aes(x = lon, y = lat)) +
    geom_path(aes(group = track), col = "#225887", size = 0.3, alpha = 0.5) +
    # scale_color_continuous_sequential(palette = "Emrld", na.value = "#011789", rev = FALSE) +
    geom_path(data = ev_eddies, aes(group = track), size = 0.3, col = "#e2fc69") +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 colour = "red3", size = 0.4, fill = NA) +
    scale_alpha_continuous(range = c(0.5, 1.0)) +
    guides(alpha = "none",
           size = "none",
           colour = "none") +
    theme_map() +
    plot.parameters

  # ggplot2::ggsave("eddies/Fig018_Eddy_trajectories_AC.jpg",
  #                 width = 3.5 * (1/3), height = 2.6 * (1/3), scale = 3.7)

  return(eddy_plt)
}

AC.Fig018 <- eddy_plot(eddies, "AC", AC.layers, AC.ev)
BC.Fig018 <- eddy_plot(eddies, "BC", BC.layers, BC.ev)
EAC.Fig018 <- eddy_plot(eddies, "EAC", EAC.layers, EAC.ev)
GS.Fig018 <- eddy_plot(eddies, "GS", GS.layers, GS.ev)
KC.Fig018 <- eddy_plot(eddies, "KC", KC.layers, KC.ev)

Fig018 <- ggarrange(AC.Fig018,
                    BC.Fig018,
                    EAC.Fig018,
                    GS.Fig018,
                    KC.Fig018,
                    ncol = 1, nrow = 5, labels = list("c", "g", "k", "o", "s"))
ggplot2::ggsave("eddies/Fig018_Eddy_trajectories.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
