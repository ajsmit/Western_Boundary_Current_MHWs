library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))

# read in the region-specific components
source("setup/plot.layers.R")

# Calculte the number of events per annum

eventCumIntFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(intensity_cumulative)) %>%
    dplyr::ungroup() %>%
    as.tibble()
}

# Do the annual CumInts

AC_annualCumInt <- eventCumIntFun(AC_events)
BC_annualCumInt <- eventCumIntFun(BC_events)
EAC_annualCumInt <- eventCumIntFun(EAC_events)
KC_annualCumInt <- eventCumIntFun(KC_events)
GS_annualCumInt <- eventCumIntFun(GS_events)

# Source the regression function
source("setup/functions.R")

# Calculate the trend for CumInts per year
AC_annualCumIntTrend <- ddply(AC_annualCumInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualCumIntTrend <- ddply(BC_annualCumInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualCumIntTrend <- ddply(EAC_annualCumInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualCumIntTrend <- ddply(KC_annualCumInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualCumIntTrend <- ddply(GS_annualCumInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1, limits = c(-60, 60)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C.days/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(24, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig012 <- gv.plot(AC_annualCumIntTrend, AC.layers)
BC.fig012 <- gv.plot(BC_annualCumIntTrend, BC.layers)
EAC.fig012 <- gv.plot(EAC_annualCumIntTrend, EAC.layers)
KC.fig012 <- gv.plot(KC_annualCumIntTrend, KC.layers)
GS.fig012 <- gv.plot(GS_annualCumIntTrend, GS.layers)

l <- get_legend(AC.fig012)
fig012 <- ggarrange(AC.fig12,
                    BC.fig12,
                    EAC.fig12,
                    GS.fig12,
                    KC.fig12,
                    ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig012_OISST_CumIntTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
