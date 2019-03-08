library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
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

# Calculte the number of events per annum

eventMaxIntFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(intensity_max)) %>%
    dplyr::ungroup() %>%
    as.tibble()
}

# Do the annual MaxInts

AC_annualMaxInt <- eventMaxIntFun(AC_events)
BC_annualMaxInt <- eventMaxIntFun(BC_events)
EAC_annualMaxInt <- eventMaxIntFun(EAC_events)
KC_annualMaxInt <- eventMaxIntFun(KC_events)
GS_annualMaxInt <- eventMaxIntFun(GS_events)

# Calculate the trend for MaxInts per year
AC_annualMaxIntTrend <- ddply(AC_annualMaxInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualMaxIntTrend <- ddply(BC_annualMaxInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualMaxIntTrend <- ddply(EAC_annualMaxInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualMaxIntTrend <- ddply(KC_annualMaxInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualMaxIntTrend <- ddply(GS_annualMaxInt, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

# read in the region-specific components
source("setup/plot.layers.R")

gv.plot <- function(data, plot.parameters, bathy) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1_51,
                         values = scales::rescale(c(min(data$slope), 0, max(data$slope)))) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "white", size = 0.2) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C/dec]",
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

AC.fig011 <- gv.plot(AC_annualMaxIntTrend, AC.layers, AC_bathy)
BC.fig011 <- gv.plot(BC_annualMaxIntTrend, BC.layers, BC_bathy)
EAC.fig011 <- gv.plot(EAC_annualMaxIntTrend, EAC.layers, EAC_bathy)
KC.fig011 <- gv.plot(KC_annualMaxIntTrend, KC.layers, KC_bathy)
GS.fig011 <- gv.plot(GS_annualMaxIntTrend, GS.layers, GS_bathy)

fig011 <- ggarrange(AC.fig011,
                    BC.fig011,
                    EAC.fig011,
                    GS.fig011,
                    KC.fig011,
                    ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig011_OISST_MaxIntTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

