# Fig00_Aviso+_EKEMean.R

# Create the figure for the Aviso+ mean velocity plots.

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)


# Functions: velocities, EKE, etc. ----------------------------------------

source("setup/functions.R")
source("setup/plot.layers.R")


# Read in monthly Aviso+ --------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/mean"
mAC.sl <- fread(paste0(csvDir, "/mean_AC_sla-19930101-20181023.csv"))
mBC.sl <- fread(paste0(csvDir, "/mean_BC_sla-19930101-20181023.csv"))
mEAC.sl <- fread(paste0(csvDir, "/mean_EAC_sla-19930101-20181023.csv"))
mKC.sl <- fread(paste0(csvDir, "/mean_KC_sla-19930101-20181023.csv"))
mGS.sl <- fread(paste0(csvDir, "/mean_GS_sla-19930101-20181023.csv"))


# Make a function to create the basic plot components ---------------------


gv.plot <- function(data, plot.parameters) {
  current_uv_scalar <- 2.25
  vec <- oceFun2(data)
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = eke)) +
    scale_fill_gradientn(colours = rev(rainbow(7, end = 4/6)),
                         space = "Lab", limits = c(0, 1)) +
    geom_segment(data = vec,
                 aes(xend = lon + ugos * current_uv_scalar,
                     yend = lat + vgos * current_uv_scalar, alpha = velocity),
                 colour = "black",
                 arrow = arrow(angle = 20, length = unit(0.1, "cm"), type = "open"),
                 linejoin = "mitre", size = 0.35, show.legend = FALSE) +
    scale_alpha_continuous(range = c(0, 1.0)) + theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(EKE~(m^2~s^{-2})),
                                  frame.colour = "black",
                                  frame.linewidth = 0.2,
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

AC.fig00 <- gv.plot(mAC.sl, AC.layers) + theme(legend.position = "bottom")
# ggplot2::ggsave("Changing_nature_of_ocean_heat/figures/Fig00_Aviso_EKEMean_AC.png", width = 8.25, height = 5.5, scale = 0.7)
BC.fig00 <- gv.plot(mBC.sl, BC.layers)
EAC.fig00 <- gv.plot(mEAC.sl, EAC.layers)
KC.fig00 <- gv.plot(mKC.sl, KC.layers)
GS.fig00 <- gv.plot(mGS.sl, GS.layers)

l <- get_legend(AC.fig00)
fig00 <- ggarrange(AC.fig00 + theme(legend.position = 'none'),
                   BC.fig00 + theme(legend.position = 'none'),
                   EAC.fig00 + theme(legend.position = 'none'),
                   l,
                   KC.fig00 + theme(legend.position = 'none'),
                   GS.fig00 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig00,
                top = text_grob("Aviso+ EKE (mean: 1993-01-01 to 2018-10-23)"))
ggplot2::ggsave("figures/Fig00_Aviso_EKEMean .pdf", width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)


