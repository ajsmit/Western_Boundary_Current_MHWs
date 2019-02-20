# Fig01_Aviso+_VelocityMean.R

# Create the figure for the Aviso+ mean velocity plots.

library(data.table)
library(tidyverse)
library(ggpubr)
library(colorspace)

# Source the regression function
source("setup/functions.R")

# read in the region-specific components
source("setup/plot.layers.R")


# Read in monthly Aviso+ --------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/mean"
mAC.sl <- fread(paste0(csvDir, "/mean_AC_sla-19930101-20181023.csv"))
mBC.sl <- fread(paste0(csvDir, "/mean_BC_sla-19930101-20181023.csv"))
mEAC.sl <- fread(paste0(csvDir, "/mean_EAC_sla-19930101-20181023.csv"))
mKC.sl <- fread(paste0(csvDir, "/mean_KC_sla-19930101-20181023.csv"))
mGS.sl <- fread(paste0(csvDir, "/mean_GS_sla-19930101-20181023.csv"))


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))


# Make a function to create the basic plot components ---------------------


gv.plot <- function(data, plot.parameters, bathy) {
  current_uv_scalar <- 2.25
  vec <- oceFun2(data)
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = velocity)) +
    scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(300, 3000, 6000)) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    # geom_segment(data = vec,
    #              aes(xend = lon + ugos * current_uv_scalar,
    #                  yend = lat + vgos * current_uv_scalar, alpha = velocity),
    #              colour = "black",
    #              arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "open"),
    #              linejoin = "mitre", size = 0.15, show.legend = FALSE) +
    scale_alpha_continuous(range = c(0, 1.0)) + theme_map() +
    guides(alpha = "none",
           size = "none",
           fill = guide_colourbar(title = expression(Velocity~(m~s^{-1})),
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(28, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.fig01 <- gv.plot(mAC.sl, AC.layers, AC_bathy)
BC.fig01 <- gv.plot(mBC.sl, BC.layers, BC_bathy)
EAC.fig01 <- gv.plot(mEAC.sl, EAC.layers, EAC_bathy)
KC.fig01 <- gv.plot(mKC.sl, KC.layers, KC_bathy)
GS.fig01 <- gv.plot(mGS.sl, GS.layers, GS_bathy)

fig01 <- ggarrange(AC.fig01,
                   BC.fig01,
                   EAC.fig01,
                   GS.fig01,
                   KC.fig01,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("figures/Fig01_Aviso+_VelocityMean.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
