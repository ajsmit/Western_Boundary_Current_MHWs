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
           fill = guide_colourbar(title = "Trend: Mean intensity\n[°C/dec]",
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

AC.fig010 <- gv.plot(AC_annualMeanIntTrend, AC.layers)
BC.fig010 <- gv.plot(BC_annualMeanIntTrend, BC.layers)
EAC.fig010 <- gv.plot(EAC_annualMeanIntTrend, EAC.layers)
KC.fig010 <- gv.plot(KC_annualMeanIntTrend, KC.layers)
GS.fig010 <- gv.plot(GS_annualMeanIntTrend, GS.layers)

l <- get_legend(AC.fig010)
fig010 <- ggarrange(AC.fig010 + theme(legend.position = 'none'),
                    BC.fig010 + theme(legend.position = 'none'),
                    EAC.fig010 + theme(legend.position = 'none'),
                    l,
                    KC.fig010 + theme(legend.position = 'none'),
                    GS.fig010 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig010,
                top = text_grob("OISST MHW Mean intensity (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig010_OISST_MeanIntTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)

