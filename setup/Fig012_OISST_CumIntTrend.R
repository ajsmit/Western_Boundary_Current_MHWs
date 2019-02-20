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


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

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

gv.plot <- function(data, plot.parameters, bathy) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey50", size = 0.2) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C.days/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(20, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.fig012 <- gv.plot(AC_annualCumIntTrend, AC.layers, AC_bathy)
BC.fig012 <- gv.plot(BC_annualCumIntTrend, BC.layers, BC_bathy)
EAC.fig012 <- gv.plot(EAC_annualCumIntTrend, EAC.layers, EAC_bathy)
KC.fig012 <- gv.plot(KC_annualCumIntTrend, KC.layers, KC_bathy)
GS.fig012 <- gv.plot(GS_annualCumIntTrend, GS.layers, GS_bathy)

fig012 <- ggarrange(AC.fig012,
                    BC.fig012,
                    EAC.fig012,
                    GS.fig012,
                    KC.fig012,
                    ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig012_OISST_CumIntTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
