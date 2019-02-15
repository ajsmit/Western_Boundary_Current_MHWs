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

# calculte the number of events per annum
eventDurFun <- function(events) {
  dur <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(duration))
}

# Do the annual counts
AC_annualCount <- eventDurFun(AC_events)
BC_annualCount <- eventDurFun(BC_events)
EAC_annualCount <- eventDurFun(EAC_events)
KC_annualCount <- eventDurFun(KC_events)
GS_annualCount <- eventDurFun(GS_events)

# Source the regression function
source("setup/functions.R")

# Calculate the trend for counts per year
AC_annualDurationTrend <- ddply(AC_annualCount, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualDurationTrend <- ddply(BC_annualCount, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualDurationTrend <- ddply(EAC_annualCount, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualDurationTrend <- ddply(KC_annualCount, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualDurationTrend <- ddply(GS_annualCount, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    scale_fill_gradientn(colours = col1, limits = c(-13, 13)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "Trend: duration\n[days/dec]",
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

AC.fig013 <- gv.plot(AC_annualDurationTrend, AC.layers)
BC.fig013 <- gv.plot(BC_annualDurationTrend, BC.layers)
EAC.fig013 <- gv.plot(EAC_annualDurationTrend, EAC.layers)
KC.fig013 <- gv.plot(KC_annualDurationTrend, KC.layers)
GS.fig013 <- gv.plot(GS_annualDurationTrend, GS.layers)

l <- get_legend(AC.fig013)
fig013 <- ggarrange(AC.fig013 + theme(legend.position = 'none'),
                    BC.fig013 + theme(legend.position = 'none'),
                    EAC.fig013 + theme(legend.position = 'none'),
                    l,
                    KC.fig013 + theme(legend.position = 'none'),
                    GS.fig013 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig013,
                top = text_grob("OISST MHW duration (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig013_OISST_DurationTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
