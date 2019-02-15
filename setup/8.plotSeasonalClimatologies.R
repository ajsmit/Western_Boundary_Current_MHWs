# 8.plotSeasonalClimatologies.R

library(doMC) # multi-core
require(mgcv)
library(readr)
library(dplyr)
library(plyr)
library(lubridate)
library(tidyverse)
library(viridis)


# Source the regression function
source("setup/functions.R")
source("setup/theme.R")

# Make world coastline
library(PBSmapping)
gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"
lats <- c(-90, 90)
lons <- c(0, 360)
shore <- importGSHHS(paste0(gshhsDir, "/gshhs_i.b"), xlim = lons, ylim = lats, maxLevel = 1, useWest = FALSE)

# Load the data
csvDir <- "/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily"
AC <- read_csv(paste0(csvDir, "/AC-avhrr-only-v2.19810901-20171019.csv"))
BC <- read_csv(paste0(csvDir, "/BC-avhrr-only-v2.19810901-20171019.csv"))
EAC <- read_csv(paste0(csvDir, "/EAC-avhrr-only-v2.19810901-20171019.csv"))
KC <- read_csv(paste0(csvDir, "/KC-avhrr-only-v2.19810901-20171019.csv"))
GS <- read_csv(paste0(csvDir, "/GS-avhrr-only-v2.19810901-20171019.csv"))

# a function to create monthly data from daily data
daily2DJF_JJA <- function(dailyData) {
  seasonData <- dailyData %>%
    dplyr::mutate(t = month(t)) %>%
    dplyr::filter(t %in% c(12, 1, 2, 6, 7, 8)) %>%
    dplyr::mutate(t = dplyr::recode_factor(t, `12` = "DJF", `1` = "DJF", `2` = "DJF",
                                           `6` = "JJA", `7` = "JJA", `8` = "JJA")) %>%
    dplyr::group_by(lon, lat, t) %>%
    dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
    dplyr::ungroup()
  return(seasonData)
}

AC_seas <- daily2DJF_JJA(AC)

AC_fig1 <- ggplot(AC_seas, aes(x = lon, y = lat)) +
  facet_wrap(~ t) +
  geom_raster(aes(fill = temp), interpolate = FALSE) +
  scale_fill_viridis(name = "°C") +
  geom_vline(xintercept = seq(10, 45, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-45, -20, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(-1.875, 53.125), ylim = c(-52.5, -12.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = ""),) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = ""),) +
  labs(title = "A", x = NULL, y = NULL)
