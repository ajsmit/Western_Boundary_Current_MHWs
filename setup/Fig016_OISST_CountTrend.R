library(data.table)
library(tidyverse)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# read in the region-specific components
source("setup/plot.layers.R")


# Define functions --------------------------------------------------------

# calculte the number of events per annum
# (events per year per decade)
eventCntFun <- function(events) {
  freq <- events %>%
    dplyr::filter(event == TRUE) %>%
    dplyr::mutate(t = fastDate(t)) %>%
    dplyr::mutate(year = year(t)) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = n())
}


# Read in OISST trend data ------------------------------------------------

# ... and do the annual counts
csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_proto_events.csv")))
AC_annualCount <- eventCntFun(AC_events); rm(AC_events)
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_proto_events.csv")))
BC_annualCount <- eventCntFun(BC_events); rm(BC_events)
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_proto_events.csv")))
EAC_annualCount <- eventCntFun(EAC_events); rm(EAC_events)
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_proto_events.csv")))
KC_annualCount <- eventCntFun(KC_events); rm(KC_events)
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_proto_events.csv")))
GS_annualCount <- eventCntFun(GS_events); rm(GS_events)


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

# Source the regression function
source("setup/functions.R")

# Calculate the trend for counts per year
AC_annualCountTrend <- ddply(AC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
BC_annualCountTrend <- ddply(BC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
EAC_annualCountTrend <- ddply(EAC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
KC_annualCountTrend <- ddply(KC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
GS_annualCountTrend <- ddply(GS_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

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
           fill = guide_colourbar(title = "[events/dec]",
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

AC.fig016 <- gv.plot(AC_annualCountTrend, AC.layers, AC_bathy)
BC.fig016 <- gv.plot(BC_annualCountTrend, BC.layers, BC_bathy)
EAC.fig016 <- gv.plot(EAC_annualCountTrend, EAC.layers, EAC_bathy)
KC.fig016 <- gv.plot(KC_annualCountTrend, KC.layers, KC_bathy)
GS.fig016 <- gv.plot(GS_annualCountTrend, GS.layers, GS_bathy)

fig016 <- ggarrange(AC.fig016,
                    BC.fig016,
                    EAC.fig016,
                    GS.fig016,
                    KC.fig016,
                    ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig016_OISST_CountTrend_v2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
