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

Duration <- function(events) {
  dur <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(duration, na.rm = TRUE)) %>%
    as.tibble()
}

AC_Duration <- Duration(AC_events)
BC_Duration <- Duration(BC_events)
EAC_Duration <- Duration(EAC_events)
KC_Duration <- Duration(KC_events)
GS_Duration <- Duration(GS_events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuBuGn", na.value = "#011789", rev = TRUE) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[days/yr]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(27.5, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig09 <- gv.plot(AC_Duration, AC.layers)
BC.fig09 <- gv.plot(BC_Duration, BC.layers)
EAC.fig09 <- gv.plot(EAC_Duration, EAC.layers)
KC.fig09 <- gv.plot(KC_Duration, KC.layers)
GS.fig09 <- gv.plot(GS_Duration, GS.layers)

fig09 <- ggarrange(AC.fig09,
                   BC.fig09,
                   EAC.fig09,
                   GS.fig09,
                   KC.fig09,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig09_OISST_Duration.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
