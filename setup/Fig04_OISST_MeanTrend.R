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

gv.plot <- function(data, plot.parameters) {
  data$p_bracket <- cut(data$p_trend, breaks = c(0, 0.001, 0.01, 0.05, 1))
  data$DT_interval <- round(data$DT_model, 1)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_gradientn(colours = col1, limits = c(-1.1, 1.1), breaks = c(-1.0, 0, 1.0)) +
    geom_contour(aes(x = lon, y = lat, z = p_trend),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C/dec]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    theme_map() +
    plot.parameters
  return(fig)
}

AC.fig04 <- gv.plot(AC_trend, AC.layers)
BC.fig04 <- gv.plot(BC_trend, BC.layers)
EAC.fig04 <- gv.plot(EAC_trend, EAC.layers)
KC.fig04 <- gv.plot(KC_trend, KC.layers)
GS.fig04 <- gv.plot(GS_trend, GS.layers)

# l <- get_legend(AC.fig04)
fig04 <- ggarrange(AC.fig04,
                   BC.fig04,
                   EAC.fig04,
                   GS.fig04,
                   KC.fig04,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig04_OISST_MeanTrend.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
