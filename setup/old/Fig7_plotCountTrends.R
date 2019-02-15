# Fig7_plotCountTrends.R

library(readr)
library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
require(rgeos); require(maptools) # maptools must be loaded after rgeos
library(doMC); doMC::registerDoMC(cores = 4)

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

# Load the data
AC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/AC-avhrr-only-v2.19810901-20171019_events.csv")
BC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/BC-avhrr-only-v2.19810901-20171019_events.csv")
EAC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/EAC-avhrr-only-v2.19810901-20171019_events.csv")
KC_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/KC-avhrr-only-v2.19810901-20171019_events.csv")
GS_events <- read_csv("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs/GS-avhrr-only-v2.19810901-20171019_events.csv")

# calculte the number of events per annum
eventCntFun <- function(events) {
  freq <- events %>%
    dplyr::mutate(year = year(date_start)) %>%
    dplyr::select(-date_start, -date_stop, -date_peak) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(y = n())
}

# Do the annual counts
AC_annualCount <- eventCntFun(AC_events)
BC_annualCount <- eventCntFun(BC_events)
EAC_annualCount <- eventCntFun(EAC_events)
KC_annualCount <- eventCntFun(KC_events)
GS_annualCount <- eventCntFun(GS_events)

# Source the regression function
source("setup/functions.R")

# Calculate the trend for counts per year
AC_annualCountTrend <- ddply(AC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
BC_annualCountTrend <- ddply(BC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
EAC_annualCountTrend <- ddply(EAC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
KC_annualCountTrend <- ddply(KC_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)
GS_annualCountTrend <- ddply(GS_annualCount, .(lon, lat), linFun, poissan = TRUE, .parallel = TRUE)

# AC_annualCountTrend$slope <- round(AC_annualCountTrend$slope, 1)
AC_fig7 <- ggplot(AC_annualCountTrend, aes(x = lon, y = lat)) +
  coord_fixed(ratio = 1, xlim = c(-1.875, 53.125), ylim = c(-52.5, -12.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  # labs(title = "A", x = NULL, y = NULL) +
  labs(title = "Trend: count/dec (1983-2012)", x = NULL, y = NULL) +
  geom_raster(aes(fill = slope), interpolate = FALSE) +
  scale_fill_gradient2(name = "Trend:\ncount\n(/dec)", low = "#243c80", mid = "gray98",
                       high = "#b6272f", midpoint = 0) +
  geom_contour(aes(x = lon, y = lat, z = pval), binwidth = 0.05, breaks = c(0.05), colour = "grey20", size = 0.2) +
  geom_vline(xintercept = seq(10, 45, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-45, -20, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "floralwhite", size = 0.2)

ggplot2::ggsave("figures/Fig7_trendCount_AC.png", width = 8.25, height = 5.5, scale = 0.7)

BC_annualCountTrend$slope <- round(BC_annualCountTrend$slope, 1)
BC_fig7 <- ggplot(BC_annualCountTrend, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = slope), interpolate = FALSE) +
  scale_fill_gradient2(name = "Trend:\ncount\n(/dec)", low = "navy", mid = "gray98",
                       high = "red4", midpoint = 0, limits = c(-1, 3)) +
  geom_contour(aes(x = lon, y = lat, z = pval), binwidth = 0.05, breaks = c(0.05), colour = "red4", size = 0.2) +
  # geom_point(aes(x = lon, y = lat, colour = pbracket), shape = "/", size = 0.25) +
  # scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
  #                     values = c("black", "black", "black", NA),
  #                     name = "p-value", guide = FALSE) +
  geom_vline(xintercept = seq(300, 335, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-40, -10, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(290, 345), ylim = c(-45, -5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  labs(title = "B", x = NULL, y = NULL)

EAC_annualCountTrend$slope <- round(annualCountTrend$slope, 1)
EAC_fig7 <- ggplot(EAC_annualCountTrend, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = slope), interpolate = FALSE) +
  scale_fill_gradient2(name = "Trend:\ncount\n(/dec)", low = "navy", mid = "gray98",
                       high = "red4", midpoint = 0, limits = c(-1, 3)) +
  geom_contour(aes(x = lon, y = lat, z = pval), binwidth = 0.05, breaks = c(0.05), colour = "red4", size = 0.2) +
  # geom_point(aes(x = lon, y = lat, colour = pbracket), shape = "/", size = 0.25) +
  # scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
  #                     values = c("black", "black", "black", NA),
  #                     name = "p-value", guide = FALSE) +
  geom_vline(xintercept = c(145, 150, 155, 160), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-40, -15, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(125, 180), ylim = c(-48.75, -8.75), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°S", sep = "")) +
  labs(title = "C", x = NULL, y = NULL)

KC_annualCountTrend$slope <- round(annualCountTrend$slope, 1)
KC_fig7 <- ggplot(KC_annualCountTrend, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = slope), interpolate = FALSE) +
  scale_fill_gradient2(name = "Trend:\ncount\n(/dec)", low = "navy", mid = "gray98",
                       high = "red4", midpoint = 0, limits = c(-1, 3)) +
  geom_contour(aes(x = lon, y = lat, z = pval), binwidth = 0.05, breaks = c(0.05), colour = "red4", size = 0.2) +
  # geom_point(aes(x = lon, y = lat, colour = pbracket), shape = "/", size = 0.25) +
  # scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
  #                     values = c("black", "black", "black", NA),
  #                     name = "p-value", guide = FALSE) +
  geom_vline(xintercept = seq(125, 170, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(20, 45, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(120, 175), ylim = c(12.5, 52.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°N", sep = "")) +
  labs(title = "D", x = NULL, y = NULL)

GS_annualCountTrend$slope <- round(annualCountTrend$slope, 1)
GS_fig7 <- ggplot(GS_annualCountTrend, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = slope), interpolate = FALSE) +
  scale_fill_gradient2(name = "Trend:\ncount\n(/dec)", low = "navy", mid = "gray98",
                       high = "red4", midpoint = 0, limits = c(-1, 3)) +
  geom_contour(aes(x = lon, y = lat, z = pval), binwidth = 0.05, breaks = c(0.05), colour = "red4", size = 0.2) +
  # geom_point(aes(x = lon, y = lat, colour = pbracket), shape = "/", size = 0.25) +
  # scale_colour_manual(breaks = c("(0,0.001]", "(0.001,0.01]", "(0.01,0.05]", "(0.05,1]"),
  #                     values = c("black", "black", "black", NA),
  #                     name = "p-value", guide = FALSE) +
  geom_vline(xintercept = seq(270, 320, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(20, 50, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(267.5, 322.5), ylim = c(15, 55), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°N", sep = "")) +
  labs(title = "E", x = NULL, y = NULL)

# ggplot(data = shore, aes(x = X, y = Y, group = PID)) +
#   geom_polygon() + coord_fixed(ratio = 1) +
#   scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

library(ggpubr)
l <- get_legend(AC_fig7)
ggarrange(AC_fig7 + theme(legend.position = 'none'),
          BC_fig7 + theme(legend.position = 'none'),
          EAC_fig7 + theme(legend.position = 'none'),
          KC_fig7 + theme(legend.position = 'none'),
          GS_fig7 + theme(legend.position = 'none'),
          l,
          ncol = 2, nrow = 3)
ggplot2::ggsave("figures/Fig7_trendCount.pdf", width = 8.25, height = 9, scale = 1.1)
