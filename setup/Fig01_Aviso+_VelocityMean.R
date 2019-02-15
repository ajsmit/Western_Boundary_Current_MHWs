# Fig01_Aviso+_VelocityMean.R

# Create the figure for the Aviso+ mean velocity plots.

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Source the regression function
source("setup/functions.R")


# Read in monthly Aviso+ --------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/mean"
mAC.sl <- fread(paste0(csvDir, "/mean_AC_sla-19930101-20181023.csv"))
mBC.sl <- fread(paste0(csvDir, "/mean_BC_sla-19930101-20181023.csv"))
mEAC.sl <- fread(paste0(csvDir, "/mean_EAC_sla-19930101-20181023.csv"))
mKC.sl <- fread(paste0(csvDir, "/mean_KC_sla-19930101-20181023.csv"))
mGS.sl <- fread(paste0(csvDir, "/mean_GS_sla-19930101-20181023.csv"))


# Make a function to create the basic plot components ---------------------

# read in the region-specific components
source("setup/plot.layers.R")

gv.plot <- function(data, plot.parameters) {
  current_uv_scalar <- 2.25
  vec <- oceFun2(data)
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = velocity)) +
    scale_fill_gradientn(colours = rev(rainbow(7, end = 4/6)),
                         space = "Lab", limits = c(0, 1.45)) +
    geom_segment(data = vec,
                 aes(xend = lon + u * current_uv_scalar,
                     yend = lat + v * current_uv_scalar, alpha = velocity),
                 colour = "black",
                 arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "open"),
                 linejoin = "mitre", size = 0.35, show.legend = FALSE) +
    scale_alpha_continuous(range = c(0, 1.0)) + theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(Velocity~(m~s^{-1})),
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

AC.fig01 <- gv.plot(mAC.sl, AC.layers) + theme(legend.position = "bottom")
BC.fig01 <- gv.plot(mBC.sl, BC.layers)
EAC.fig01 <- gv.plot(mEAC.sl, EAC.layers)
KC.fig01 <- gv.plot(mKC.sl, KC.layers)
GS.fig01 <- gv.plot(mGS.sl, GS.layers)

l <- get_legend(AC.fig01)
fig01 <- ggarrange(AC.fig01 + theme(legend.position = 'none'),
                   BC.fig01 + theme(legend.position = 'none'),
                   EAC.fig01 + theme(legend.position = 'none'),
                   l,
                   KC.fig01 + theme(legend.position = 'none'),
                   GS.fig01 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig01,
                top = text_grob("Aviso+ surface currents (mean: 1993-01-01 to 2018-10-23)"))
ggplot2::ggsave("figures/Fig01_Aviso_VelocityMean.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
