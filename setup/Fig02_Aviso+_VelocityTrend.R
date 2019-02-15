# Fig02_Aviso+_VelocityTrend.R

library(data.table)
library(tidyverse)
# The devtool version of ggplot2 has some nifty vector plotting options, that can be used to
# create very professional-looking vectors (see especially the) linejoin = "mitre" argument
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(ggpubr)


# Read in Aviso trend data ------------------------------------------------

inDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/monthly-rate"
AC_trend <- as.tibble(fread(paste0(inDir, "/AC-velocity-rate-19930101-20181001.csv")))
BC_trend <- as.tibble(fread(paste0(inDir, "/BC-velocity-rate-19930101-20181001.csv")))
EAC_trend <- as.tibble(fread(paste0(inDir, "/EAC-velocity-rate-19930101-20181001.csv")))
EAC_trend <- EAC_trend %>%
  filter(lat <= -15)
KC_trend <- as.tibble(fread(paste0(inDir, "/KC-velocity-rate-19930101-20181001.csv")))
GS_trend <- as.tibble(fread(paste0(inDir, "/GS-velocity-rate-19930101-20181001.csv")))

# read in the region-specific components
source("setup/plot.layers.R")


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_gradientn(colours = col1,
                         limit = c(-0.2, 0.2)) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(Velocity~trend~(m~s^{-1}~dec^{-1})),
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

AC.fig02 <- gv.plot(AC_trend, AC.layers)
BC.fig02 <- gv.plot(BC_trend, BC.layers)
EAC.fig02 <- gv.plot(EAC_trend, EAC.layers)
KC.fig02 <- gv.plot(KC_trend, KC.layers)
GS.fig02 <- gv.plot(GS_trend, GS.layers)

l <- get_legend(AC.fig02)
fig02 <- ggarrange(AC.fig02 + theme(legend.position = 'none'),
                   BC.fig02 + theme(legend.position = 'none'),
                   EAC.fig02 + theme(legend.position = 'none'),
                   l,
                   KC.fig02 + theme(legend.position = 'none'),
                   GS.fig02 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig02,
                top = text_grob("Aviso+ surface current velocity (trend: 1993-01-01 to 2018-10-23)"))
ggplot2::ggsave("figures/Fig02_Aviso_VelocityTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
