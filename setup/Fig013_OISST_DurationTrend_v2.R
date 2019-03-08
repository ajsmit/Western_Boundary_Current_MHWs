# Fig013_OISST_DurationTrend_v2.R

library(data.table)
library(tidyverse)
library(colorspace)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# read in the region-specific components
source("setup/plot.layers.R")

# Source the regression function
source("setup/functions.R")


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))


# Define functions --------------------------------------------------------

eventDurFun <- function(events) {
  dur <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_end, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = mean(duration))
}

# Do the annual counts
AC_annualDur <- eventDurFun(AC_events)
BC_annualDur <- eventDurFun(BC_events)
EAC_annualDur <- eventDurFun(EAC_events)
KC_annualDur <- eventDurFun(KC_events)
GS_annualDur <- eventDurFun(GS_events)

# Calculate the trend for counts per year
AC_annualDurationTrend <- ddply(AC_annualDur, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
BC_annualDurationTrend <- ddply(BC_annualDur, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
EAC_annualDurationTrend <- ddply(EAC_annualDur, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
KC_annualDurationTrend <- ddply(KC_annualDur, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)
GS_annualDurationTrend <- ddply(GS_annualDur, .(lon, lat), linFun, poissan = FALSE, .parallel = TRUE)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy, region) {

  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))

  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke <- fortify(mask.list$mke90)
  eke <- fortify(mask.list$eke90)
  int <- fortify(mask.list$int90)

  data$slope <- round(data$slope, 2)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope * 10)) +
    # scale_fill_gradientn(colours = col1,
    #                      values = scales::rescale(c(min(data$slope), 0, max(data$slope)))) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE, n_interp = 11) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_contour(aes(x = lon, y = lat, z = pval),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "white", size = 0.2) +
    # geom_polygon(data = int, aes(long, lat, group = group),
    #              fill = NA, colour = "purple", size = 0.3) +
    # geom_polygon(data = mke, aes(long, lat, group = group),
    #              fill = NA, colour = "red3", size = 0.3) +
    # geom_polygon(data = eke, aes(long, lat, group = group),
    #              fill = NA, colour = "navy", size = 0.3) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "(days/dec)",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(25, units = "mm"),
                                  draw.ulim = FALSE,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.Fig013_v2 <- gv.plot(AC_annualDurationTrend, AC.layers, AC_bathy, "AC")
BC.Fig013_v2 <- gv.plot(BC_annualDurationTrend, BC.layers, BC_bathy, "BC")
EAC.Fig013_v2 <- gv.plot(EAC_annualDurationTrend, EAC.layers, EAC_bathy, "EAC")
GS.Fig013_v2 <- gv.plot(GS_annualDurationTrend, GS.layers, GS_bathy, "GS")
KC.Fig013_v2 <- gv.plot(KC_annualDurationTrend, KC.layers, KC_bathy, "KC")

Fig013_v2 <- ggarrange(AC.Fig013_v2,
                       BC.Fig013_v2,
                       EAC.Fig013_v2,
                       GS.Fig013_v2,
                       KC.Fig013_v2,
                       ncol = 1, nrow = 5, labels = "auto")
ggplot2::ggsave("figures/Fig013_OISST_DurationTrend_v2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
