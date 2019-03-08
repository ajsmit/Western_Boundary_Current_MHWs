# Calculate the mean (long-term) MKE, EKE and mean intensity for masking

library(data.table)
library(raster)
library(sf)
library(smoothr)
library(tidyverse)
library(ncdf4)
library(RcppRoll)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

mask.fun <- function(region) {
  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(4)

  # Prepare the data.table
  dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                         lat = nc$dim$lat$vals,
                         t = nc$dim$time$vals)

  ke <- as.data.table(melt(ugos, value.name = "ugos"))
  ke[, c("vgos", "ugosa", "vgosa") := .( # slow
    # as.Date(t, origin = "1950-01-01 00:00:00"), # not necessary to convert as no date calcs needed
    as.vector(ncvar_get(nc, varid = "vgos")) * 100,
    as.vector(ncvar_get(nc, varid = "ugosa")) * 100,
    as.vector(ncvar_get(nc, varid = "vgosa")) * 100
  )
  ]
  nc_close(nc)

  geo <- ke %>%
    na.omit() %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))) %>%
    dplyr::filter(eke > -999) %>%
    dplyr::summarise(mke90 = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke90 = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mke75 = mke90,
                  eke75 = eke90) %>%
    data.table()

  # Read in OISST data
  csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
  events <- fread(paste0(csvDir, "/", region, "-avhrr-only-v2.19810901-20180930_events.csv"))
  int <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(int90 = mean(intensity_mean, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(int75 = int90) %>%
    data.table()

  # Join them
  out <- na.omit(int[geo, on=c("lon==lon", "lat==lat")])

  # return(out)
  # define the regions to mask based on MKE and EKE thresholds
  out$int75[which(out$int75 >= quantile(out$int75, probs = 0.75))] <- -999
  out$mke75[which(out$mke75 >= quantile(out$mke75, probs = 0.75))] <- -999
  out$eke75[which(out$eke75 >= quantile(out$eke75, probs = 0.75))] <- -999
  out$int90[which(out$int90 >= quantile(out$int90, probs = 0.9))] <- -999
  out$mke90[which(out$mke90 >= quantile(out$mke90, probs = 0.9))] <- -999
  out$eke90[which(out$eke90 >= quantile(out$eke90, probs = 0.9))] <- -999

  coordinates(out) <- ~ lon + lat
  gridded(out) <- TRUE
  proj4string(out) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # out.st <- stack(out)
  int75 <- raster(out, layer = "int75")
  mke75 <- raster(out, layer = "mke75")
  eke75 <- raster(out, layer = "eke75")
  int90 <- raster(out, layer = "int90")
  mke90 <- raster(out, layer = "mke90")
  eke90 <- raster(out, layer = "eke90")

  int75.mask <- smooth(rasterToPolygons(int75, fun = function(int75) { int75 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  mke75.mask <- smooth(rasterToPolygons(mke75, fun = function(mke75) { mke75 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  eke75.mask <- smooth(rasterToPolygons(eke75, fun = function(eke75) { eke75 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  int90.mask <- smooth(rasterToPolygons(int90, fun = function(int90) { int90 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  mke90.mask <- smooth(rasterToPolygons(mke90, fun = function(mke90) { mke90 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  eke90.mask <- smooth(rasterToPolygons(eke90, fun = function(eke90) { eke90 == -999 }, dissolve = TRUE),
                       method = "ksmooth", smoothness = 3)
  mask.list <- list(int75 = int75.mask,
                    mke75 = mke75.mask,
                    eke75 = eke75.mask,
                    int90 = int90.mask,
                    mke90 = mke90.mask,
                    eke90 = eke90.mask)
  save(mask.list, file = paste0("/Users/ajsmit/Dropbox/R/WBCs/masks/", region, "-mask_polys.RData"))
  return(mask.list)
}

AC.mask <- mask.fun("AC")
BC.mask <- mask.fun("BC")
EAC.mask <- mask.fun("EAC")
GS.mask <- mask.fun("GS")
KC.mask <- mask.fun("KC")


# Make some plots ---------------------------------------------------------

# read in the region-specific components
source("setup/plot.layers.R")

# plot function
# plot.parameters <- EAC.layers
# bathy <- AC_bathy
# data <- mask.list
# str(data$int75)
# region <- "EAC"

poly.plot <- function(region, plot.parameters) {
  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))
  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke90 <- fortify(mask.list$mke90)
  eke90 <- fortify(mask.list$eke90)
  int90 <- fortify(mask.list$int90)
  mke75 <- fortify(mask.list$mke75)
  eke75 <- fortify(mask.list$eke75)
  int75 <- fortify(mask.list$int75)
  plt90 <- ggplot() +
    geom_polygon(data = mke90, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    geom_polygon(data = int90, aes(long, lat, group = group),
                 fill = alpha("#c83c7e", 1.0), colour = NA) +
    geom_polygon(data = eke75, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "navy", size = 0.3,
                 linetype = "dashed") +
    geom_polygon(data = eke90, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "navy", size = 0.3) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    guides(alpha = "none",
           fill = "none",
           size = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters

  # plt75 <- ggplot() +
  #   geom_polygon(data = int75, aes(long, lat, group = group),
  #                fill = alpha("#c83c7e", 1.0), colour = NA) +
  #   geom_polygon(data = eke75, aes(long, lat, group = group),
  #                fill = alpha("gray50", 0.5), colour = "navy", size = 0.3) +
  #   geom_polygon(data = mke75, aes(long, lat, group = group),
  #                fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
  #   stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
  #                col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  #   guides(alpha = "none",
  #          fill = "none",
  #          size = "none") +
  #   theme_map() + labs(x = NULL, y = NULL) +
  #   plot.parameters

  # plt.list <- list(plt90, plt75)
  return(plt90)
}

AC.poly.plt <- poly.plot("AC", AC.layers)
BC.poly.plt <- poly.plot("BC", BC.layers)
EAC.poly.plt <- poly.plot("EAC", EAC.layers)
GS.poly.plt <- poly.plot("GS", GS.layers)
KC.poly.plt <- poly.plot("KC", KC.layers)

fig_masks <- ggarrange(AC.poly.plt,
                       BC.poly.plt,
                       EAC.poly.plt,
                       GS.poly.plt,
                       KC.poly.plt,
                       ncol = 1, nrow = 5, labels = list("a", "b", "c", "d", "e"))
ggplot2::ggsave("publ_plots/Figure_2_75eke.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
