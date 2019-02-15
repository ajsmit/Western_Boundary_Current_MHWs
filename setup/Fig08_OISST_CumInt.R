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

# read in the region-specific components
source("setup/plot.layers.R")

# Calculate the total maximum intensity

CumInt <- function(events) {
  cum <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(intensity_cumulative, na.rm = TRUE)) %>%
    as.tibble()
}

AC_CumInt <- CumInt(AC_events)
BC_CumInt <- CumInt(BC_events)
EAC_CumInt <- CumInt(EAC_events)
KC_CumInt <- CumInt(KC_events)
GS_CumInt <- CumInt(GS_events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuRd", na.value = "#011789", rev = TRUE) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C.days]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(26, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig08 <- gv.plot(AC_CumInt, AC.layers)
BC.fig08 <- gv.plot(BC_CumInt, BC.layers)
EAC.fig08 <- gv.plot(EAC_CumInt, EAC.layers)
KC.fig08 <- gv.plot(KC_CumInt, KC.layers)
GS.fig08 <- gv.plot(GS_CumInt, GS.layers)

fig08 <- ggarrange(AC.fig08,
                   BC.fig08,
                   EAC.fig08,
                   GS.fig08,
                   KC.fig08,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig08_OISST_CumInt.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
