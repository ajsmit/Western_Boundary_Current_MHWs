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

eventOnsetFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(rate_onset)) %>%
    dplyr::ungroup() %>%
    as.tibble()
}

# Do the annual Onsets

AC_annualOnset <- eventOnsetFun(AC_events)
BC_annualOnset <- eventOnsetFun(BC_events)
EAC_annualOnset <- eventOnsetFun(EAC_events)
KC_annualOnset <- eventOnsetFun(KC_events)
GS_annualOnset <- eventOnsetFun(GS_events)

# Source the regression function
source("setup/functions.R")

# Calculate the trend for Onsets per year
AC_annualOnsetTrend <- ddply(AC_annualOnset, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualOnsetTrend <- ddply(BC_annualOnset, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualOnsetTrend <- ddply(EAC_annualOnset, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualOnsetTrend <- ddply(KC_annualOnset, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualOnsetTrend <- ddply(GS_annualOnset, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


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
           fill = guide_colourbar(title = "[°C/d/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(21, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.fig014 <- gv.plot(AC_annualOnsetTrend, AC.layers, AC_bathy)
BC.fig014 <- gv.plot(BC_annualOnsetTrend, BC.layers, BC_bathy)
EAC.fig014 <- gv.plot(EAC_annualOnsetTrend, EAC.layers, EAC_bathy)
KC.fig014 <- gv.plot(KC_annualOnsetTrend, KC.layers, KC_bathy)
GS.fig014 <- gv.plot(GS_annualOnsetTrend, GS.layers, GS_bathy)

fig014 <- ggarrange(AC.fig014,
                    BC.fig014,
                    EAC.fig014,
                    GS.fig014,
                    KC.fig014,
                    ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig014_OISST_OnsetTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
