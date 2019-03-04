library(data.table)
library(tidyverse)
library(colorspace)
library(ggpubr)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))


# Misc setup --------------------------------------------------------------

# read in the region-specific components
source("setup/plot.layers.R")

# Calculate the total count

totalCntFun <- function(events) {
  freq <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = n())
}

AC_totalCount <- totalCntFun(AC_events)
BC_totalCount <- totalCntFun(BC_events)
EAC_totalCount <- totalCntFun(EAC_events)
KC_totalCount <- totalCntFun(KC_events)
GS_totalCount <- totalCntFun(GS_events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy, region) {
  load(paste0("/Users/ajsmit/Dropbox/R/WBCs/masks/", region, "-masks.RData"))
  eke <- fortify(mask.list$eke90)
  mke <- fortify(mask.list$mke90)
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuBuGn", na.value = "#011789", rev = TRUE) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = eke, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "navy", size = 0.3) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[n]",
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

AC.fig05 <- gv.plot(AC_totalCount, AC.layers, AC_bathy, "AC")
BC.fig05 <- gv.plot(BC_totalCount, BC.layers, BC_bathy, "BC")
EAC.fig05 <- gv.plot(EAC_totalCount, EAC.layers, EAC_bathy, "EAC")
GS.fig05 <- gv.plot(GS_totalCount, GS.layers, GS_bathy, "GS")
KC.fig05 <- gv.plot(KC_totalCount, KC.layers, KC_bathy, "KC")

# l <- get_legend(AC.fig05)
fig05 <- ggarrange(AC.fig05,
                   BC.fig05,
                   EAC.fig05,
                   GS.fig05,
                   KC.fig05,
                   ncol = 1, nrow = 5, labels = list("a", "c", "e", "g", "i"))
# annotate_figure(fig05,
# top = text_grob("OISST MHW count (total: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig05_OISST_TotalCount.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
