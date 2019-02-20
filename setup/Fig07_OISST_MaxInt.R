library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggpubr)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

# read in the region-specific components
source("setup/plot.layers.R")

# Calculate the total maximum intensity

MaxInt <- function(events) {
  totMax <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(intensity_max, na.rm = TRUE)) %>%
    as.tibble()
}

AC_MaxInt <- MaxInt(AC_events)
BC_MaxInt <- MaxInt(BC_events)
EAC_MaxInt <- MaxInt(EAC_events)
KC_MaxInt <- MaxInt(KC_events)
GS_MaxInt <- MaxInt(GS_events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuRd", na.value = "#011789", rev = TRUE) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.fig07 <- gv.plot(AC_MaxInt, AC.layers, AC_bathy)
BC.fig07 <- gv.plot(BC_MaxInt, BC.layers, BC_bathy)
EAC.fig07 <- gv.plot(EAC_MaxInt, EAC.layers, EAC_bathy)
KC.fig07 <- gv.plot(KC_MaxInt, KC.layers, KC_bathy)
GS.fig07 <- gv.plot(GS_MaxInt, GS.layers, GS_bathy)

fig07 <- ggarrange(AC.fig07,
                   BC.fig07,
                   EAC.fig07,
                   GS.fig07,
                   KC.fig07,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig07_OISST_MaxInt.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
