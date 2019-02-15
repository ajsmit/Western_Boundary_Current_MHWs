# Fig03_Aviso+_EKETrend.R

library(data.table)
library(tidyverse)
# The devtool version of ggplot2 has some nifty vector plotting options, that can be used to
# create very professional-looking vectors (see especially the) linejoin = "mitre" argument
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(ggpubr)


# Read in Aviso trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/monthly-rate"
AC_trend <- as.tibble(fread(paste0(csvDir, "/AC-eke-rate-19930101-20181001.csv")))
BC_trend <- as.tibble(fread(paste0(csvDir, "/BC-eke-rate-19930101-20181001.csv")))
EAC_trend <- as.tibble(fread(paste0(csvDir, "/EAC-eke-rate-19930101-20181001.csv")))
EAC_trend <- EAC_trend %>%
  filter(lat <= -15)
KC_trend <- as.tibble(fread(paste0(csvDir, "/KC-eke-rate-19930101-20181001.csv")))
GS_trend <- as.tibble(fread(paste0(csvDir, "/GS-eke-rate-19930101-20181001.csv")))

# read in the region-specific components
source("setup/plot.layers.R")


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_gradientn(colours = col1,
                         limit = c(-0.15, 0.15)) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(EKE~trend~(m^2~s^{-2}~dec^{-1})),
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

AC.fig03 <- gv.plot(AC_trend, AC.layers)
BC.fig03 <- gv.plot(BC_trend, BC.layers)
EAC.fig03 <- gv.plot(EAC_trend, EAC.layers)
KC.fig03 <- gv.plot(KC_trend, KC.layers)
GS.fig03 <- gv.plot(GS_trend, GS.layers)

l <- get_legend(AC.fig03)
fig03 <- ggarrange(AC.fig03 + theme(legend.position = 'none'),
                   BC.fig03 + theme(legend.position = 'none'),
                   EAC.fig03 + theme(legend.position = 'none'),
                   l,
                   KC.fig03 + theme(legend.position = 'none'),
                   GS.fig03 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig03,
                top = text_grob("Aviso+ eddy kinetic energy (trend: 1993-01-01 to 2018-10-23)"))
ggplot2::ggsave("figures/Fig03_Aviso_EKETrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
