# Fig018_MODIS_Aqua_MeanTrend.R

library(data.table)
library(tidyverse)
# The devtool version of ggplot2 has some nifty vector plotting options, that can be used to
# create very professional-looking vectors (see especially the) linejoin = "mitre" argument
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(ggpubr)


# read in MODIS Aqua data -------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/regridded"
AC.chl <- fread(paste0(csvDir, "/AC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
BC.chl <- fread(paste0(csvDir, "/BC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
EAC.chl <- fread(paste0(csvDir, "/EAC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
KC.chl <- fread(paste0(csvDir, "/KC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
GS.chl <- fread(paste0(csvDir, "/GS.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))


# output Dir for chlor_a trend --------------------------------------------

outDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/monthly-rate"


# calc time series of monthly means ---------------------------------------

daily2monthly <- function(dailyData) {
  monthlyData <- dailyData %>%
    dplyr::mutate(t = fastDate(time)) %>%
    dplyr::mutate(month = floor_date(time, unit = "month")) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(chlor_a = mean(chlor_a, na.rm = TRUE)) %>%
    dplyr::mutate(year = year(month)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(num = seq(1:length(temp))) %>%
    dplyr::ungroup()
  return(monthlyData)
}


# function for calculating chl-a trend ------------------------------------

chlLinFun <- function(annualEvent) {
  mod <- lm(y ~ year, data = annualEvent)
  trend <- data.frame(slope = summary(mod)$coefficients[2,1] * 10, pval = summary(mod)$coefficients[2,4])
  # trend <- summary(mod)$coefficients[2,c(1,4)]
  trend$pbracket <- cut(trend$pval, breaks = c(0, 0.001, 0.01, 0.05, 1))
  return(trend)
}

# read in the region-specific components ----------------------------------

source("setup/plot.layers.R")


# make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters) {
  data$p_bracket <- cut(data$p_trend, breaks = c(0, 0.001, 0.01, 0.05, 1))
  data$DT_interval <- round(data$DT_model, 1)

  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = DT_model)) +
    scale_fill_gradientn(colours = col1, limits = c(-1.0, 1.0)) +
    geom_contour(aes(x = lon, y = lat, z = p_trend),
                 binwidth = 0.05, breaks = c(0.05),
                 colour = "grey20", size = 0.2) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "SST trend\n[°C/dec]",
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

AC.fig04 <- gv.plot(AC_trend, AC.layers)
BC.fig04 <- gv.plot(BC_trend, BC.layers)
EAC.fig04 <- gv.plot(EAC_trend, EAC.layers)
KC.fig04 <- gv.plot(KC_trend, KC.layers)
GS.fig04 <- gv.plot(GS_trend, GS.layers)

l <- get_legend(AC.fig04)
fig04 <- ggarrange(AC.fig04 + theme(legend.position = 'none'),
                   BC.fig04 + theme(legend.position = 'none'),
                   EAC.fig04 + theme(legend.position = 'none'),
                   l,
                   KC.fig04 + theme(legend.position = 'none'),
                   GS.fig04 + theme(legend.position = 'none'),
                   ncol = 2, nrow = 3)
annotate_figure(fig04,
                top = text_grob("OISST (trend: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig04_OISST_MeanTrend.jpg",
                width = 7.0 * (1/3), height = 8.3 * (1/3), scale = 3.5)
