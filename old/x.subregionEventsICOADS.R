# 6.subregionEvents.R

require(mgcv)
library(readr)
library(dplyr)
library(plyr)
library(lubridate)
library(tidyverse)
library(data.table)
library(timeSeries)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)
source("setup/theme.R")
theme_set(grey_update) # Set theme
source("setup/functions.R")

samples <- data.frame(region = rep("AC", 6),
                      type = c("current", "current", "current", "control", "control", "control"),
                      ID = c("sub1", "sub2", "sub3", "sub4", "sub5", "sub6"),
                      lonmin = c(23.5, 22, 20.5, 25, 26.5, 28.5),
                      lonmax = c(25, 23.5, 22, 26.5, 28, 30),
                      latmin = c(-36.5, -36.5, -38, -40, -38, -40),
                      latmax = c(-35, -35, -36.5, -38.5, -36.5, -38.5))

ICOADS.dir <- "/Volumes/Benguela/OceanData/ICOADS/ICOADS_1_Degree_Enhanced"

nc <- nc_open(paste0(ICOADS.dir, "/", "sst.mean.nc"))
LatIdx <- which(nc$dim$lat$vals %in% c(-34.5, -35.5, -37.5))
LonIdx <- which(nc$dim$lon$vals %in% c(24.5, 22.5, 21.5))
time <- seq(1:length(nc$dim$time$vals))
sst.nc <- ncvar_get(nc,
                 varid = "sst",
                 start = c(LonIdx[1], LatIdx[1], time[1]),
                 count = c(length(LonIdx), length(LatIdx), length(time))) %>%
  round(4)
dimnames(sst.nc) <- list(lon = nc$dim$lon$vals[LonIdx],
                      lat = nc$dim$lat$vals[LatIdx],
                      time = nc$dim$time$vals[time])

sst <-
  as.data.frame(melt(sst.nc, value.name = "temp"), row.names = NULL) %>%
  mutate(time = as.Date(time, origin = "1800-1-1 00:00:0.0")) %>%
  group_by(time) %>%
  summarise(temp = mean(temp, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(anom = temp - mean(temp, na.rm = TRUE))
summary(sst)

nc_close(nc)

anomaly.gam <- gam(anom ~ s(c(1:length(anom))), data = sst)
anomaly.pred <- data.frame(date = sst$time[1:length(anomaly.gam$y)],
                           anom = sst$anom,
                           fitted = anomaly.gam$fitted.values)

# Fit a gam and detect MHWs in the deseasoned and detrended data
# (Probably remove the seasons before detrending...)
anomaly.detr <- data.frame(doy = anomaly$doy[1:length(anomaly.gam$y)],
                           date = anomaly$date[1:length(anomaly.gam$y)],
                           temp = residuals(anomaly.gam))

# Plot the deseasoned data
anomaly.pl <- ggplot(anomaly.pred, aes(x = date, y = anom)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  geom_line(col = "grey25", size = 0.2) +
  geom_line(aes(y = fitted, x = date), col = "white", alpha = 1.0, size = 0.9) +
  xlab(NULL) +
  ylab("Temperature anomaly") + ggtitle("a. 'De-seasoned'")

# Plot the deseasoned and detrended data
anomaly.detr.pl <- ggplot(anomaly.detr, aes(x = date, y = temp)) +
  geom_line(col = "grey25", size = 0.2) +
  geom_line(data = anomaly.detr.ev$clim, aes(x = date, y = seas_clim_year),
            col = "blue3", alpha = 1, size = 0.5) +
  geom_line(data = anomaly.detr.ev$clim, aes(x = date, y = thresh_clim_year), col = "red", alpha = 0.7) +
  geom_rug(data = anomaly.detr.ev$event, aes(x = date_start, y = int_max_abs), col = "red", sides = "b") +
  xlab(NULL) +
  ylab("Temperature anomaly") + ggtitle("b. 'De-seasoned' and detrended")

ggarrange(anomaly.pl + theme(legend.position = 'none'),
          anomaly.detr.pl + theme(legend.position = 'none'),
          ncol = 1, nrow = 2)
ggplot2::ggsave("figures/Fig13_AC_trendPlots_current.pdf", width = 7, height = 3.5, scale = 1.4)

