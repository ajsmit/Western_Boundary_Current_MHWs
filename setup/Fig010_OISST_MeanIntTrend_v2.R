library(data.table)
library(tidyverse)
library(ggpubr)
library(colorspace)
library(doMC); doMC::registerDoMC(cores = 4)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as_tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as_tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as_tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as_tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as_tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))

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

gv.plot <- function(data, plot.parameters, bathy, region) {

  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))

  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke <- fortify(mask.list$mke90)
  eke <- fortify(mask.list$eke90)
  int <- fortify(mask.list$int90)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = slope)) +
    # scale_fill_gradientn(colours = col1,
    #                      values = scales::rescale(c(min(data$slope), 0, max(data$slope)))) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE, n_interp = 11) +
    # geom_contour(aes(x = lon, y = lat, z = pval),
    #              binwidth = 0.05, breaks = c(0.05),
    #              colour = "white", size = 0.2) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = NA, colour = "purple", size = 0.3) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 fill = NA, colour = "red3", size = 0.3) +
    geom_polygon(data = eke, aes(long, lat, group = group),
                 fill = NA, colour = "navy", size = 0.3) +
    guides(fill = guide_colourbar(title = "(°C/dec)",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28, units = "mm"),
                                  draw.ulim = FALSE,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.Fig010 <- gv.plot(AC_annualMeanIntTrend, AC.layers, AC_bathy, "AC")
BC.Fig010 <- gv.plot(BC_annualMeanIntTrend, BC.layers, BC_bathy, "BC")
EAC.Fig010 <- gv.plot(EAC_annualMeanIntTrend, EAC.layers, EAC_bathy, "EAC")
GS.Fig010 <- gv.plot(GS_annualMeanIntTrend, GS.layers, GS_bathy, "GS")
KC.Fig010 <- gv.plot(KC_annualMeanIntTrend, KC.layers, KC_bathy, "KC")

Fig010 <- ggarrange(AC.Fig010,
                    BC.Fig010,
                    EAC.Fig010,
                    GS.Fig010,
                    KC.Fig010,
                    ncol = 1, nrow = 5, labels = "auto")
ggplot2::ggsave("figures/Fig010_OISST_MeanIntTrend_v2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

