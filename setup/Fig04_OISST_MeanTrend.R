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


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy) {
  data$p_bracket <- cut(data$p_trend, breaks = c(0, 0.001, 0.01, 0.05, 1))
  data$DT_interval <- round(data$DT_model, 1)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_gradientn(colours = col1, limits = c(-1.1, 1.1), breaks = c(-1.0, 0, 1.0)) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_contour(aes(x = lon, y = lat, z = p_trend),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "white", size = 0.2) +
    guides(alpha = "none",
           size = "none",
           fill = guide_colourbar(title = "[°C/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() +
    plot.parameters
  return(fig)
}

AC.fig04 <- gv.plot(AC_trend, AC.layers, AC_bathy)
BC.fig04 <- gv.plot(BC_trend, BC.layers, BC_bathy)
EAC.fig04 <- gv.plot(EAC_trend, EAC.layers, EAC_bathy)
GS.fig04 <- gv.plot(GS_trend, GS.layers, GS_bathy)
KC.fig04 <- gv.plot(KC_trend, KC.layers, KC_bathy)

fig04 <- ggarrange(AC.fig04,
                   BC.fig04,
                   EAC.fig04,
                   GS.fig04,
                   KC.fig04,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig04_OISST_MeanTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
