#

library(ncdf4)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(data.table)
library(colorspace)
library(OneR) # for binning


# Miscellaneous things to load --------------------------------------------

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# Read in the region-specific components
source("setup/regionDefinition.R")
source("setup/plot.layers.R")


# The data ----------------------------------------------------------------

# eddies
# for converting julian date to calendar date, see:
# https://stackoverflow.com/questions/20789612/convert-julian-date-to-calendar-dates-within-a-data-frame-in-r

ncFile <- "~/lustre/Aviso/global/delayed-time/value-added/eddy-trajectory/eddy_trajectory_2.0exp_19930101_20180118.nc"
nc <- nc_open(ncFile)
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


# The function ------------------------------------------------------------

# define even dates for each region
AC.ev <- c("2003-06-25", "2003-07-12", "2014-08-02", "2014-09-15", "2018-03-01", "2018-03-17")
BC.ev <- c("2015-06-11", "2015-10-03", "2017-02-17", "2017-08-15", "2014-08-30", "2015-01-02")
EAC.ev <- c("2005-10-21", "2005-11-20", "2000-12-04", "2000-12-26", "2017-01-07", "2017-01-18")
GS.ev <- c("2014-03-01", "2014-04-2", "2014-11-23", "2015-02-24", "2018-03-05", "2018-05-02")
KC.ev <- c("2013-03-31", "2013-05-05", "2007-06-04", "2007-06-15", "2013-09-09", "2013-09-15")

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
  pts <- fread(paste0("setup/masks/", region, "-eddy_pts.csv"))

  tracks <- dplyr::right_join(eddies, pts) %>%
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

  eddy_plt <- ggplot(data = WBC_eddies, aes(x = lon, y = lat)) +
    geom_path(aes(group = track, col = duration), size = 0.3, alpha = 0.5) +
    scale_color_continuous_sequential(palette = "Emrld", na.value = "#011789", rev = FALSE) +
    geom_path(data = ev_eddies, aes(group = track), size = 0.4, col = "red") +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z, alpha = ..level..),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    scale_alpha_continuous(range = c(0.5, 1.0)) +
    guides(alpha = "none",
           size = "none",
           colour = guide_colourbar(title = "Duration",
                                    frame.colour = "black",
                                    frame.linewidth = 0.4,
                                    ticks.colour = "black",
                                    barheight = unit(2, units = "mm"),
                                    barwidth = unit(28, units = "mm"),
                                    draw.ulim = F,
                                    title.position = 'left',
                                    title.hjust = 0.5,
                                    label.hjust = 0.5)) +
    theme_map() +
    plot.parameters

  return(eddy_plt)
}

plt1 <- eddy_plot(eddies, "AC", AC.layers, AC.ev)
plt2 <- eddy_plot(eddies, "BC", BC.layers, BC.ev)
plt3 <- eddy_plot(eddies, "EAC", EAC.layers, EAC.ev)
plt4 <- eddy_plot(eddies, "GS", GS.layers, GS.ev)
plt5 <- eddy_plot(eddies, "KC", KC.layers, KC.ev)

fig <- ggarrange(plt1,
                   plt2,
                   plt3,
                   plt4,
                   plt5,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("setup/eddies/WBC_eddies.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
