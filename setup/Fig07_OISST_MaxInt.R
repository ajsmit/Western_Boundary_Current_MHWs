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

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_viridis(limit = c(0, 5.5)) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "Mean maximum intensity\n[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(70, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'top',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig07 <- gv.plot(AC_MaxInt, AC.layers)
BC.fig07 <- gv.plot(BC_MaxInt, BC.layers)
EAC.fig07 <- gv.plot(EAC_MaxInt, EAC.layers)
KC.fig07 <- gv.plot(KC_MaxInt, KC.layers)
GS.fig07 <- gv.plot(GS_MaxInt, GS.layers)

l <- get_legend(AC.fig07)
fig07 <- ggarrange(AC.fig07 + theme(legend.position = 'none'),
                   BC.fig07 + theme(legend.position = 'none'),
                   EAC.fig07 + theme(legend.position = 'none'),
                   l,
                   KC.fig07 + theme(legend.position = 'none'),
                   GS.fig07 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig07,
                top = text_grob("OISST MHW maximum intensity (mean: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig07_OISST_MaxInt.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
