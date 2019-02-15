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

gv.plot <- function(data, plot.parameters) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1, limits = c(-0.8, 0.8)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "Trend: rate of onset\n[°C/d/dec]",
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

AC.fig014 <- gv.plot(AC_annualOnsetTrend, AC.layers)
BC.fig014 <- gv.plot(BC_annualOnsetTrend, BC.layers)
EAC.fig014 <- gv.plot(EAC_annualOnsetTrend, EAC.layers)
KC.fig014 <- gv.plot(KC_annualOnsetTrend, KC.layers)
GS.fig014 <- gv.plot(GS_annualOnsetTrend, GS.layers)

l <- get_legend(AC.fig014)
fig014 <- ggarrange(AC.fig014 + theme(legend.position = 'none'),
                    BC.fig014 + theme(legend.position = 'none'),
                    EAC.fig014 + theme(legend.position = 'none'),
                    l,
                    KC.fig014 + theme(legend.position = 'none'),
                    GS.fig014 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig014,
                top = text_grob("OISST MHW rate of onset (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig014_OISST_OnsetTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
