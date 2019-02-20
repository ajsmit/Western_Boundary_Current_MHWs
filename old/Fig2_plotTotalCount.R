# Fig2_countTrends.R

library(readr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
require(rgeos); require(maptools) # maptools must be loaded after rgeos

theme_set(theme_grey())
grey_update <- theme_grey() +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)
  )
theme_set(grey_update) # Set theme

# Make world coastline
library(PBSmapping)
gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"
lats <- c(-90, 90)
lons <- c(0, 360)
shore <- importGSHHS(paste0(gshhsDir, "/gshhs_i.b"), xlim = lons, ylim = lats, maxLevel = 1, useWest = FALSE)

AC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/AC-avhrr-only-v2.19810901-20171019_events.csv")
BC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/BC-avhrr-only-v2.19810901-20171019_events.csv")
EAC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/EAC-avhrr-only-v2.19810901-20171019_events.csv")
KC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/KC-avhrr-only-v2.19810901-20171019_events.csv")
GS_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/GS-avhrr-only-v2.19810901-20171019_events.csv")

# calculate the total count
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

# Agulhas Current
AC_fig2 <- ggplot(AC_totalCount, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y), interpolate = FALSE) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = seq(10, 45, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-45, -20, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "floralwhite", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(-1.875, 53.125), ylim = c(-52.5, -12.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  # labs(title = "A", x = NULL, y = NULL)
  labs(title = "Total event count (1983-2012)", x = NULL, y = NULL)

ggplot2::ggsave("figures/Fig2_totalCount_AC.png", width = 8.25, height = 5.5, scale = 0.7)

# Brazil Current
BC_fig2 <- ggplot(BC_totalCount, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y), interpolate = FALSE) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = seq(300, 335, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-40, -10, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(290, 345), ylim = c(-45, -5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  labs(title = "B", x = NULL, y = NULL)

# East Australian Current
EAC_fig2 <- ggplot(EAC_totalCount, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y), interpolate = FALSE) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = c(145, 150, 155, 160), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-40, -15, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(125, 180), ylim = c(-48.75, -8.75), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  labs(title = "C", x = NULL, y = NULL)

# Kuroshio Current
KC_fig2 <- ggplot(KC_totalCount, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y), interpolate = FALSE) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = seq(125, 170, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(20, 45, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(120, 175), ylim = c(12.5, 52.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°N", sep = "")) +
  labs(title = "D", x = NULL, y = NULL)

# Gulf Stream Current
GS_fig2 <- ggplot(GS_totalCount, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y), interpolate = FALSE) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = seq(270, 320, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(20, 50, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(267.5, 322.5), ylim = c(15, 55), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°N", sep = "")) +
  labs(title = "E", x = NULL, y = NULL)

library(ggpubr)
l <- get_legend(AC_fig2)
ggarrange(AC_fig2 + theme(legend.position = 'none'),
          BC_fig2 + theme(legend.position = 'none'),
          EAC_fig2 + theme(legend.position = 'none'),
          KC_fig2 + theme(legend.position = 'none'),
          GS_fig2 + theme(legend.position = 'none'),
          l,
          ncol = 2, nrow = 3)
ggplot2::ggsave("figures/Fig2_totalCount.pdf", width = 8.25, height = 9, scale = 1.1)
