# Fig017_MODIS_Aqua_ChlMean.R

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# source the regression function
source("setup/functions.R")


# read in monthly MODIS chl-a ---------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/regridded"
AC.chl <- fread(paste0(csvDir, "/AC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"))
BC.chl <- fread(paste0(csvDir, "/BC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"))
EAC.chl <- fread(paste0(csvDir, "/EAC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"))
KC.chl <- fread(paste0(csvDir, "/KC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"))
GS.chl <- fread(paste0(csvDir, "/GS.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"))


# function to calculate mean MODIS chl-a ----------------------------------

chlFun <- function(data) {
  out <- data %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(chlor_a = round(mean(chlor_a, na.rm = TRUE), 4)) %>%
    dplyr::ungroup()
  return(out)
}


# calc mean chl-a ---------------------------------------------------------

mAC.chl <- chlFun(AC.chl)
mBC.chl <- chlFun(BC.chl)
mEAC.chl <- chlFun(EAC.chl)
mGS.chl <- chlFun(GS.chl)
mKC.chl <- chlFun(KC.chl)


# read in the region-specific components ----------------------------------

source("setup/plot.layers.R")


# make a function to create the basic plot components ---------------------

library(sp) # for the colours...
cols <- get_col_regions()

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = chlor_a)) +
    scale_fill_gradientn(colours = cols,
                         space = "Lab", limits = c(0.025, 20),
                         breaks = c(0.02, 0.1, 1.0, 5.0, 20.0),
                         trans = 'log10') +
    guides(fill = guide_colourbar(title = "Mean chlorophyll-a\n[mg/m3]",
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
    theme_map() +
    plot.parameters
  return(fig)
}

AC.fig017 <- gv.plot(mAC.chl, AC.layers) + theme(legend.position = "bottom")
BC.fig017 <- gv.plot(mBC.chl, BC.layers)
EAC.fig017 <- gv.plot(mEAC.chl, EAC.layers)
KC.fig017 <- gv.plot(mKC.chl, KC.layers)
GS.fig017 <- gv.plot(mGS.chl, GS.layers)

l <- get_legend(AC.fig017)
fig017 <- ggarrange(AC.fig017 + theme(legend.position = 'none'),
                    BC.fig017 + theme(legend.position = 'none'),
                    EAC.fig017 + theme(legend.position = 'none'),
                    l,
                    KC.fig017 + theme(legend.position = 'none'),
                    GS.fig017 + theme(legend.position = 'none'),
                    ncol = 2, nrow = 3)
annotate_figure(fig017,
                top = text_grob("MODIS-Aqua chl-a (mean: 2002-07-05 to 2018-11-18)"))
ggplot2::ggsave("figures/Fig017_MODIS_Aqua_ChlMean.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
