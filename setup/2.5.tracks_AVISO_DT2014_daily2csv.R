# 2.5.tracks_AVISO_DT2014_daily2csv.R

library(ncdf4)
library(tidyverse)
library(lubridate)
library(data.table)
library(colorspace)


# Miscellaneous things to load --------------------------------------------

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# Read in the region-specific components
source("setup/regionDefinition.R")
source("setup/plot.layers.R")

# The bathy data
bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
bathy_AC <- fread(paste0(bathyDir, "/AC_bathy.csv"))


# Read eddy db ------------------------------------------------------------

# ncFile <- "/Users/ajsmit/Downloads/tracks_AVISO_DT2014_daily_web.nc"
# nc <- nc_open(ncFile)
# dat <- tibble(A = round(ncvar_get(nc, varid = "A"), 3),
#               cyc = ncvar_get(nc, varid = "cyc"),
#               j1 = as.Date(ncvar_get(nc, varid = "j1") - 2448623, origin = "1992-01-01"),
#               L = round(ncvar_get(nc, varid = "L"), 3),
#               lat = round(ncvar_get(nc, varid = "lat"), 4),
#               lon = round(ncvar_get(nc, varid = "lon"), 4),
#               n = ncvar_get(nc, varid = "n"),
#               track = ncvar_get(nc, varid = "track"),
#               U = round(ncvar_get(nc, varid = "U"), 3))
#
# csvDir <- "/Users/ajsmit/Downloads"



# Setup -------------------------------------------------------------------

# for converting julian date to calendar date, see:
# https://stackoverflow.com/questions/20789612/convert-julian-date-to-calendar-dates-within-a-data-frame-in-r

ncFile <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/value-added/eddy-trajectory/eddy_trajectory_2.0exp_19930101_20180118.nc"
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
# fwrite(dat,
#        file = paste0(csvDir, "/tracks_AVISO_DT2014_daily_web", ".csv", sep = ""),
#        append = FALSE, col.names = TRUE)
# rm(dat)

# some plots of eddies
# add column woth duration of each eddy
library(OneR)

eddies2 <- eddies %>%
  dplyr::filter(lat >= bbox$AC[1] & lat <= bbox$AC[2]) %>%
  dplyr::filter(lon >= bbox$AC[3] & lon <= bbox$AC[4]) %>%
  # for the original Chelton et al. 2011 version
  # dplyr::filter(lon >= bbox$AC[3] + 360 & lon <= bbox$AC[4] + 360) %>%
  # dplyr::mutate(lon = lon - 360,
  #               row_ID = seq(1:n())) %>%
  # dplyr::mutate(row_ID = seq(1:n())) %>%
  dplyr::group_by(track) %>%
  dplyr::mutate(duration = max(observation_number + 1)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(bin = bin(duration, method = "content"))


# Find eddies starting the the AC retroflection region --------------------

# find only eddies that originate in the Agulhas retroflection region:
# 15 to 25°E and -45 to -35°S
start_eddy_ID <- eddies2 %>%
  dplyr::filter(observation_number == 0) %>%
  dplyr::filter(lon >= 10 & lon <= 25) %>%
  dplyr::filter(lat >= -45 & lat <= -35) %>%
  dplyr::select(track)

start_eddies <- dplyr::right_join(eddies2, start_eddy_ID)

start_eddies %>%
  # dplyr::filter(duration <= quantile(dur, 0.25)) %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_path(aes(colour = cyclonic_type, group = track), size = 0.2) +
  stat_contour(data = bathy_AC, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "magenta", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none") +
  facet_wrap(~ bin, nrow = 3) +
  theme_map() +
  AC.layers


# Find eddies ending in the AC retroflection region -----------------------

# find eddies that end in Agulhas retroflection region
end_eddy_ID <- eddies2 %>%
  dplyr::group_by(track) %>%
  dplyr::filter(observation_number == max(observation_number)) %>%
  dplyr::filter(lon >= 10 & lon <= 25) %>%
  dplyr::filter(lat >= -45 & lat <= -35) %>%
  dplyr::ungroup() %>%
  dplyr::select(track)

end_eddies <- dplyr::right_join(eddies2, end_eddy_ID)

fig_end <- end_eddies %>%
  # dplyr::filter(duration <= quantile(dur, 0.25)) %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_path(aes(colour = cyclonic_type, group = track), size = 0.2) +
  stat_contour(data = bathy_AC, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "magenta", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none") +
  facet_wrap(~ bin, nrow = 3) +
  theme_map() +
  AC.layers

fig_end


# Eddies coinciding with MHWs ---------------------------------------------

end_eddies_ev <- end_eddies %>%
  dplyr::filter(time >= as.Date("2003-06-25") & time <= as.Date("2003-07-12") |
                  time >= as.Date("2014-08-02") & time <= as.Date("2014-09-15") |
                  time >= as.Date("2018-03-01") & time <= as.Date("2018-03-17")) %>%
  dplyr::select(track) %>%
  unique

ev_eddies <- dplyr::right_join(eddies2, end_eddies_ev)

fig_end2 <- ggplot(data = end_eddies, aes(x = lon, y = lat)) +
  geom_path(aes(group = track, col = duration), size = 0.3, alpha = 0.5) +
  scale_color_continuous_sequential(palette = "Emrld", na.value = "#011789", rev = FALSE) +
  geom_path(data = ev_eddies, aes(group = track), size = 0.4, col = "red") +
  stat_contour(data = bathy_AC, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none") +
  # facet_wrap(~ bin, nrow = 3) +
  theme_map() +
  AC.layers

fig_end2


# Find eddies that start in the AC ----------------------------------------

# find only eddies that originate in the Agulhas retroflection region:
# 15 to 25°E and -45 to -35°S
AC_pts <- fread("setup/masks/AC-eddy_pts.csv")

AC_eddy_tracks <- dplyr::right_join(eddies2, AC_pts) %>%
  dplyr::filter(observation_number == 0) %>%
  dplyr::select(track)
summary(AC_eddy_tracks)
c(AC_eddy_tracks)

AC_eddies <- eddies2 %>%
  dplyr::filter(track %in% c(AC_eddy_tracks$track)) %>%
  dplyr::mutate(cyclonic_type = as.factor(cyclonic_type))

AC_eddies %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_path(aes(colour = cyclonic_type, group = track), size = 0.2) +
  stat_contour(data = bathy_AC, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "magenta", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none") +
  facet_wrap(~ bin, nrow = 5) +
  theme_map() +
  AC.layers

AC_eddies_ev <- AC_eddies %>%
  dplyr::filter(time >= as.Date("2003-06-25") & time <= as.Date("2003-07-12") |
                  time >= as.Date("2014-08-02") & time <= as.Date("2014-09-15") |
                  time >= as.Date("2018-03-01") & time <= as.Date("2018-03-17")) %>%
  dplyr::select(track) %>%
  unique

AC_ev_eddies <- dplyr::right_join(AC_eddies, AC_eddies_ev)


AC_eddy_plt <- ggplot(data = AC_eddies, aes(x = lon, y = lat)) +
  geom_path(aes(group = track, col = duration), size = 0.3, alpha = 0.5) +
  scale_color_continuous_sequential(palette = "Emrld", na.value = "#011789", rev = FALSE) +
  geom_path(data = AC_ev_eddies, aes(group = track), size = 0.4, col = "red") +
  stat_contour(data = bathy_AC, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none",
         size = "none",
         colour = guide_colourbar(title = "Duration",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(35, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
  theme_map() +
  AC.layers

AC_eddy_plt
