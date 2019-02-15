# Make world coastline ----------------------------------------------------

require(rgeos); require(maptools) # maptools must be loaded after rgeos
library(PBSmapping)
library(ggplot2)
library(mapproj) # use the Lambert azimuthal equal-area projection

gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"
lats <- c(-90, 90)
lons <- c(0, 360)
shore <- importGSHHS(paste0(gshhsDir, "/gshhs_l.b"), xlim = lons, ylim = lats, maxLevel = 1, useWest = FALSE)


# create a theme ----------------------------------------------------------

# library(extrafont)
theme_map <- function(...) {
  theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.text = element_text(size = 13, colour = "black"),
      axis.ticks = element_line(size = 0.5),
      plot.title = element_text(size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = "#f5f5f2", color = NA),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "white", size = 0.3),
      # legend.position = "right",
      legend.position = c(0.27, 0.9),
      legend.title = element_text(size = 10),
      panel.border = element_blank(),
      ...
    )
}


# Topo map ----------------------------------------------------------------

# library(raster)
# library(rgdal)
# library(sp)
# library(dplyr)
#
# file_in <- "/Users/ajsmit/spatial/Natural_Earth_Data/SR_HR/SR_HR.tif"
# dem <- raster(file_in)
#
# # Prepare for ggplot2 (long format) and tidy up a bit
# dem_df <- as_tibble(as.data.frame(dem, xy = TRUE))
# names(dem_df) <- c("lon", "lat", "z")
# rm(dem)
xlim_AC <- c(-1.875, 53.125)
ylim_AC <- c(-52.5, -12.5)
# dem_AC <- dem_df %>%
#   filter(lon >= xlim_AC[1] & lon <=xlim_AC[2]) %>%
#   filter(lat >= ylim_AC[1] & lat <=ylim_AC[2])


# Find the subregion box extents ------------------------------------------

bx <- read.csv2("/Volumes/GoogleDrive/My Drive/WBCs/setup/subRegions.csv")
AC.bx <- c(bx[1:9, c(4:7)])
BC.bx <- c(bx[10:18, c(4:7)])
EAC.bx <- c(bx[19:27, c(4:7)])
KC.bx <- c(bx[28:36, c(4:7)])
GS.bx <- c(bx[37:45, c(4:7)])


# Set the default plot options --------------------------------------------

AC.layers = list(
  geom_vline(xintercept = seq(10, 45, 5), linetype = "dashed", size = 0.2),
  geom_hline(yintercept = seq(-45, -20, 5), linetype = "dashed", size = 0.2),
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "#272822", fill = "#c8c8c8", size = 0.4),
  coord_fixed(ratio = 1, xlim = c(-1.875, 53.125), ylim = c(-52.5, -12.5), expand = TRUE),
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "")),
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = ""))
)

BC.layers = list(
  geom_vline(xintercept = seq(300, 335, 5), linetype = "dashed", size = 0.2),
  geom_hline(yintercept = seq(-40, -10, 5), linetype = "dashed", size = 0.2),
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "#272822", fill = "#c8c8c8", size = 0.4),
  coord_fixed(ratio = 1, xlim = c(290, 345), ylim = c(-45, -5), expand = TRUE),
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "")),
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = ""))
)

EAC.layers <- list(
  geom_vline(xintercept = c(145, 150, 155, 160), linetype = "dashed", size = 0.2),
  geom_hline(yintercept = seq(-40, -15, 5), linetype = "dashed", size = 0.2),
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "#272822", fill = "#c8c8c8", size = 0.4),
  coord_fixed(ratio = 1, xlim = c(125, 180), ylim = c(-48.75, -8.75), expand = TRUE),
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "")),
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = ""))
)

KC.layers <- list(
  geom_vline(xintercept = seq(125, 170, 5), linetype = "dashed", size = 0.2),
  geom_hline(yintercept = seq(20, 45, 5), linetype = "dashed", size = 0.2),
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "#272822", fill = "#c8c8c8", size = 0.4),
  coord_fixed(ratio = 1, xlim = c(120, 175), ylim = c(12.5, 52.5), expand = TRUE),
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "")),
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°N", sep = ""))
)

GS.layers <- list(
  geom_vline(xintercept = seq(270, 320, 5), linetype = "dashed", size = 0.2),
  geom_hline(yintercept = seq(20, 50, 5), linetype = "dashed", size = 0.2),
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "#272822", fill = "#c8c8c8", size = 0.4),
  coord_fixed(ratio = 1, xlim = c(267.5, 322.5), ylim = c(15, 55), expand = TRUE),
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "")),
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°N", sep = ""))
)
