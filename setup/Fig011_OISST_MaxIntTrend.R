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

gv.plot <- function(data, plot.parameters) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1, limits = c(-1.0, 1.0)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "Trend: maximum intensity\n[°C/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.2,
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

AC.fig011 <- gv.plot(AC_annualMaxIntTrend, AC.layers)
BC.fig011 <- gv.plot(BC_annualMaxIntTrend, BC.layers)
EAC.fig011 <- gv.plot(EAC_annualMaxIntTrend, EAC.layers)
KC.fig011 <- gv.plot(KC_annualMaxIntTrend, KC.layers)
GS.fig011 <- gv.plot(GS_annualMaxIntTrend, GS.layers)

l <- get_legend(AC.fig011)
fig011 <- ggarrange(AC.fig011 + theme(legend.position = 'none'),
                    BC.fig011 + theme(legend.position = 'none'),
                    EAC.fig011 + theme(legend.position = 'none'),
                    l,
                    KC.fig011 + theme(legend.position = 'none'),
                    GS.fig011 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig011,
                top = text_grob("OISST MHW maximum intensity (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig011_OISST_MaxIntTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)

