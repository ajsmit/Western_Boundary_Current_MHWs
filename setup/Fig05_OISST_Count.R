library(data.table)
library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggpubr)


# Read in OISST trend data ------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
AC_events <- as.tibble(fread(paste0(csvDir, "/", "AC-avhrr-only-v2.19810901-20180930_events.csv")))
BC_events <- as.tibble(fread(paste0(csvDir, "/", "BC-avhrr-only-v2.19810901-20180930_events.csv")))
EAC_events <- as.tibble(fread(paste0(csvDir, "/", "EAC-avhrr-only-v2.19810901-20180930_events.csv")))
KC_events <- as.tibble(fread(paste0(csvDir, "/", "KC-avhrr-only-v2.19810901-20180930_events.csv")))
GS_events <- as.tibble(fread(paste0(csvDir, "/", "GS-avhrr-only-v2.19810901-20180930_events.csv")))

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

gv.plot <- function(data, plot.parameters) {
  fig <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = y)) +
    scale_fill_continuous_sequential(palette = "PuBuGn", na.value = "#011789", rev = TRUE) +
    theme_map() +
    guides(alpha = "none",
           fill = guide_colourbar(title = "[n]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(37, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    plot.parameters
  return(fig)
}

AC.fig05 <- gv.plot(AC_totalCount, AC.layers)
BC.fig05 <- gv.plot(BC_totalCount, BC.layers)
EAC.fig05 <- gv.plot(EAC_totalCount, EAC.layers)
KC.fig05 <- gv.plot(KC_totalCount, KC.layers)
GS.fig05 <- gv.plot(GS_totalCount, GS.layers)

# l <- get_legend(AC.fig05)
fig05 <- ggarrange(AC.fig05,
                   BC.fig05,
                   EAC.fig05,
                   GS.fig05,
                   KC.fig05,
                   ncol = 1, nrow = 5)
# annotate_figure(fig05,
# top = text_grob("OISST MHW count (total: 1981-09-01 to 2018-09-01)"))
ggplot2::ggsave("figures/Fig05_OISST_TotalCount.jpg",
                width = 3.6 * (1/3), height = 13 * (1/3), scale = 3.7)
