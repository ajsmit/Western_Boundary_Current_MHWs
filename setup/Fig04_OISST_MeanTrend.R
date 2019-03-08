library(data.table)
library(tidyverse)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
library(ggpubr)


# Read in OISST trend data ------------------------------------------------

inDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/monthly-rate"
AC_trend <- as.tibble(fread(paste0(inDir, "/AC-rate-19810901-20180901.csv")))
BC_trend <- as.tibble(fread(paste0(inDir, "/BC-rate-19810901-20180901.csv")))
EAC_trend <- as.tibble(fread(paste0(inDir, "/EAC-rate-19810901-20180901.csv")))
KC_trend <- as.tibble(fread(paste0(inDir, "/KC-rate-19810901-20180901.csv")))
GS_trend <- as.tibble(fread(paste0(inDir, "/GS-rate-19810901-20180901.csv")))

# read in the region-specific components
source("setup/plot.layers.R")


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy, region) {

  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))

  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke <- fortify(mask.list$mke90)
  eke <- fortify(mask.list$eke90)
  int <- fortify(mask.list$int90)

  data$p_bracket <- cut(data$p_trend, breaks = c(0, 0.001, 0.01, 0.05, 1))
  data$DT_interval <- round(data$DT_model, 1)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE, n_interp = 11) +
    # stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
    #              col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_contour(aes(x = lon, y = lat, z = p_trend),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "white", size = 0.2) +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = NA, colour = "purple", size = 0.3) +
    # geom_polygon(data = mke, aes(long, lat, group = group),
    #              fill = NA, colour = "red3", size = 0.3) +
    # geom_polygon(data = eke, aes(long, lat, group = group),
    #              fill = NA, colour = "navy", size = 0.3) +
    guides(alpha = "none",
           size = "none",
           fill = guide_colourbar(title = "(°C/dec)",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28.75, units = "mm"),
                                  draw.ulim = FALSE,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() +
    plot.parameters
  return(fig)
}

AC.Fig04 <- gv.plot(AC_trend, AC.layers, AC_bathy, "AC")
BC.Fig04 <- gv.plot(BC_trend, BC.layers, BC_bathy, "BC")
EAC.Fig04 <- gv.plot(EAC_trend, EAC.layers, EAC_bathy, "EAC")
GS.Fig04 <- gv.plot(GS_trend, GS.layers, GS_bathy, "GS")
KC.Fig04 <- gv.plot(KC_trend, KC.layers, KC_bathy, "KC")

Fig04 <- ggarrange(AC.Fig04,
                   BC.Fig04,
                   EAC.Fig04,
                   GS.Fig04,
                   KC.Fig04,
                   ncol = 1, nrow = 5, labels = "auto")
ggplot2::ggsave("figures/Fig04_OISST_MeanTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
