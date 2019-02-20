library(data.table)
library(tidyverse)
library(ggplot2)
library(colorspace)
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

MeanInt <- function(events) {
  totMax <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(intensity_mean, na.rm = TRUE)) %>%
    as.tibble()
}

AC_MeanInt <- MeanInt(AC_events)
BC_MeanInt <- MeanInt(BC_events)
EAC_MeanInt <- MeanInt(EAC_events)
KC_MeanInt <- MeanInt(KC_events)
GS_MeanInt <- MeanInt(GS_events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuRd", na.value = "#011789", rev = TRUE) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(35, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig06 <- gv.plot(AC_MeanInt, AC.layers)
BC.fig06 <- gv.plot(BC_MeanInt, BC.layers)
EAC.fig06 <- gv.plot(EAC_MeanInt, EAC.layers)
KC.fig06 <- gv.plot(KC_MeanInt, KC.layers)
GS.fig06 <- gv.plot(GS_MeanInt, GS.layers)

# l <- get_legend(AC.fig06)
fig06 <- ggarrange(AC.fig06,
                   BC.fig06,
                   EAC.fig06,
                   GS.fig06,
                   KC.fig06,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig06_OISST_MeanInt.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
