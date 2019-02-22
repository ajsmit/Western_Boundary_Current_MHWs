library(data.table)
library(tidyverse)
library(ggplot2)
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


# Source the regression function
source("setup/functions.R")

# read in the region-specific components
source("setup/plot.layers.R")

# Calculte the number of events per annum

eventMeanIntFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(intensity_mean)) %>%
    dplyr::ungroup() %>%
    as.tibble()
}

# Do the annual MeanInts
AC_annualMeanInt <- eventMeanIntFun(AC_events)
BC_annualMeanInt <- eventMeanIntFun(BC_events)
EAC_annualMeanInt <- eventMeanIntFun(EAC_events)
KC_annualMeanInt <- eventMeanIntFun(KC_events)
GS_annualMeanInt <- eventMeanIntFun(GS_events)

# Calculate the trend for MeanInts per year
AC_annualMeanIntTrend <- ddply(AC_annualMeanInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualMeanIntTrend <- ddply(BC_annualMeanInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualMeanIntTrend <- ddply(EAC_annualMeanInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualMeanIntTrend <- ddply(KC_annualMeanInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualMeanIntTrend <- ddply(GS_annualMeanInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)

# Make a function to create the basic plot components ---------------------

# read in the region-specific components
source("setup/plot.layers.R")

gv.plot <- function(data, plot.parameters, bathy) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1,
                         values = scales::rescale(c(min(data$slope), 0, max(data$slope)))) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "white", size = 0.2) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    guides(fill = guide_colourbar(title = "[°C/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(27, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.fig010 <- gv.plot(AC_annualMeanIntTrend, AC.layers, AC_bathy)
BC.fig010 <- gv.plot(BC_annualMeanIntTrend, BC.layers, BC_bathy)
EAC.fig010 <- gv.plot(EAC_annualMeanIntTrend, EAC.layers, EAC_bathy)
KC.fig010 <- gv.plot(KC_annualMeanIntTrend, KC.layers, KC_bathy)
GS.fig010 <- gv.plot(GS_annualMeanIntTrend, GS.layers, GS_bathy)

fig010 <- ggarrange(AC.fig010,
                   BC.fig010,
                   EAC.fig010,
                   GS.fig010,
                   KC.fig010,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig010_OISST_MeanIntTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

