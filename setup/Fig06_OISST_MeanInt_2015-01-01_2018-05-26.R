library(data.table)
library(tidyverse)
library(colorspace)
library(ggpubr)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC.events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC.events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC.events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
GS.events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))
KC.events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC.bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC.bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC.bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS.bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC.bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

# read in the region-specific components
source("setup/plot.layers.R")

# Calculate the total maximum intensity

MeanInt <- function(events) {
  totMax <- events %>%
    dplyr::filter(date_start >= "2015-07-01" & date_end <= "2018-05-26") %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(intensity_mean, na.rm = TRUE)) %>%
    as.tibble()
}

AC.MeanInt <- MeanInt(AC.events)
BC.MeanInt <- MeanInt(BC.events)
EAC.MeanInt <- MeanInt(EAC.events)
KC.MeanInt <- MeanInt(KC.events)
GS.MeanInt <- MeanInt(GS.events)


# Make a function to create the basic plot components ---------------------

gv.plot <- function(data, plot.parameters, bathy, region) {
  load(paste0("/Users/ajsmit/Dropbox/R/WBCs/masks/", region, "-mask_polys.RData"))
  int <- fortify(mask.list$int90)
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuBuGn", na.value = "#011789", rev = TRUE) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = NA, colour = "purple", size = 0.4) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(32, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters
  return(fig)
}

AC.Fig06 <- gv.plot(AC.MeanInt, AC.layers, AC.bathy, "AC")
BC.Fig06 <- gv.plot(BC.MeanInt, BC.layers, BC.bathy, "BC")
EAC.Fig06 <- gv.plot(EAC.MeanInt, EAC.layers, EAC.bathy, "EAC")
GS.Fig06 <- gv.plot(GS.MeanInt, GS.layers, GS.bathy, "GS")
KC.Fig06 <- gv.plot(KC.MeanInt, KC.layers, KC.bathy, "KC")

Fig06_sub <- ggarrange(AC.Fig06,
                   BC.Fig06,
                   EAC.Fig06,
                   GS.Fig06,
                   KC.Fig06,
                   ncol = 1, nrow = 5, labels = list("D", "H", "L", "P", "T"))
ggplot2::ggsave("figures/Fig06_OISST_MeanInt_2015-07-01_2018-05-26.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
