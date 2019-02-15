# 6.subregionEvents.R

library(doMC) # multi-core
require(mgcv)
library(readr)
library(data.table)
library(dplyr)
library(plyr)
library(lubridate)
library(tidyverse)
library(data.table)
library(timeSeries)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)
source("setup/theme.R")
source("setup/functions.R")
theme_set(grey_update) # Set theme

# sst <- read_csv("/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily/AC-avhrr-only-v2.19810901-20171019.csv")
sst <- read_csv("/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily/BC-avhrr-only-v2.19810901-20171019.csv")
# sst <- read_csv("/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily/EAC-avhrr-only-v2.19810901-20171019.csv")
# sst <- read_csv("/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily/KC-avhrr-only-v2.19810901-20171019.csv")
# sst <- read_csv("/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily/GS-avhrr-only-v2.19810901-20171019.csv")

samples <- data.frame(region = rep("AC", 6),
                      type = c("current", "current", "current", "control", "control", "control"),
                      ID = c("sub1", "sub2", "sub3", "sub4", "sub5", "sub6"),
                      lonmin = c(23.5, 22, 20.5, 25, 26.5, 28.5),
                      lonmax = c(25, 23.5, 22, 26.5, 28, 30),
                      latmin = c(-36.5, -36.5, -38, -40, -38, -40),
                      latmax = c(-35, -35, -36.5, -38.5, -36.5, -38.5))
fwrite(samples, "setup/subRegions.csv")

makeSubSample <- function(data, sampleFrame, region, ID, type = NULL) {
  bbox <- sampleFrame[sampleFrame$region == region & sampleFrame$ID == ID, c(4:7)]
  x <- data %>%
    dplyr::filter(lon > bbox$lonmin & lon < bbox$lonmax | lat > bbox$latmin & lat < bbox$latmax) %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
    dplyr::mutate(ID = ID) %>%
    dplyr::mutate(type = type) %>%
    dplyr::mutate(region = region) %>%
    dplyr::ungroup()
return(x)
}

sub1 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub1", type = "current")
sub2 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub2", type = "current")
sub3 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub3", type = "current")
sub4 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub4", type = "control")
sub5 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub5", type = "control")
sub6 <- makeSubSample(data = sst, sampleFrame = samples, region = "AC", ID = "sub6", type = "control")
sub.c <- rbind(sub1, sub2, sub3, sub4, sub5, sub6)
remove(sub1, sub2, sub3, sub4, sub5, sub6)

# Calculate mean temp across the replicate samples within the 'type' group
# (i.e. either 'control' or 'current')
raw <- sub.c %>%
  dplyr::filter(type == "current") %>%
  dplyr::group_by(type, t) %>%
  dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
  select(t, temp) %>%
  dplyr::ungroup() %>%
  make_whole()

# Detect MHWs in the raw time series
raw.ev <- RmarineHeatWaves::detect(raw, climatology_start = 1983,
                                   climatology_end = 2012,
                                   min_duration = 5)

# Fit a gam and detect MHWs in the deseasoned data
anomaly <- raw.ev$clim %>%
  dplyr::ungroup() %>%
  dplyr::mutate(anom = temp - seas_clim_year) %>%
  dplyr::select(doy, date, anom) %>%
  dplyr::rename(temp = anom)
anomaly.gam <- gam(temp ~ s(c(1:length(temp))), data = anomaly)
anomaly.pred <- data.frame(doy = anomaly$doy[1:length(anomaly.gam$y)],
                           date = anomaly$date[1:length(anomaly.gam$y)],
                           fitted = anomaly.gam$fitted.values)
anomaly.ev <- RmarineHeatWaves::detect(anomaly,
                                       climatology_start = 1983,
                                       climatology_end = 2012,
                                       min_duration = 5)

# Fit a gam and detect MHWs in the deseasoned and detrended data
# (Probably remove the seasons before detrending...)
anomaly.detr <- data.frame(doy = anomaly$doy[1:length(anomaly.gam$y)],
                           date = anomaly$date[1:length(anomaly.gam$y)],
                           temp = residuals(anomaly.gam))
anomaly.detr.ev <- RmarineHeatWaves::detect(anomaly.detr,
                                            climatology_start = 1983,
                                            climatology_end = 2012,
                                            min_duration = 5)

# Plot the deseasoned data
anomaly.pl <- ggplot(anomaly, aes(x = date, y = temp)) +
  geom_line(col = "grey25", size = 0.2) +
  geom_line(data = anomaly.ev$clim, aes(x = date, y = seas_clim_year),
            col = "blue3", alpha = 1, size = 0.5) +
  geom_line(data = anomaly.ev$clim, aes(x = date, y = thresh_clim_year), col = "red", alpha = 0.7) +
  geom_line(data = anomaly.pred, aes(y = fitted, x = date), col = "white", alpha = 1.0, size = 0.9) +
  geom_rug(data = anomaly.ev$event, aes(x = date_start, y = int_max_abs), col = "red", sides = "b") +
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

