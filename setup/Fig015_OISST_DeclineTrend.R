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


# Calculte the number of events per annum

eventDeclineFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(rate_decline, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    as.tibble()
}

# Do the annual Declines

AC_annualDecline <- eventDeclineFun(AC_events)
BC_annualDecline <- eventDeclineFun(BC_events)
EAC_annualDecline <- eventDeclineFun(EAC_events)
KC_annualDecline <- eventDeclineFun(KC_events)
GS_annualDecline <- eventDeclineFun(GS_events)

# Source the regression function
source("setup/functions.R")

# Calculate the trend for Declines per year
AC_annualDeclineTrend <- ddply(AC_annualDecline, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualDeclineTrend <- ddply(BC_annualDecline, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualDeclineTrend <- ddply(EAC_annualDecline, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualDeclineTrend <- ddply(KC_annualDecline, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualDeclineTrend <- ddply(GS_annualDecline, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

# read in the region-specific components
source("setup/plot.layers.R")

gv.plot <- function(data, plot.parameters) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1, limits = c(-0.36, 0.36)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "Trend: rate of decline\n[°C/d/dec]",
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

AC.fig015 <- gv.plot(AC_annualDeclineTrend, AC.layers)
BC.fig015 <- gv.plot(BC_annualDeclineTrend, BC.layers)
EAC.fig015 <- gv.plot(EAC_annualDeclineTrend, EAC.layers)
KC.fig015 <- gv.plot(KC_annualDeclineTrend, KC.layers)
GS.fig015 <- gv.plot(GS_annualDeclineTrend, GS.layers)

l <- get_legend(AC.fig015)
fig015 <- ggarrange(AC.fig015 + theme(legend.position = 'none'),
                    BC.fig015 + theme(legend.position = 'none'),
                    EAC.fig015 + theme(legend.position = 'none'),
                    l,
                    KC.fig015 + theme(legend.position = 'none'),
                    GS.fig015 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig015,
                top = text_grob("OISST rate of decline (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig015_OISST_DeclineTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
