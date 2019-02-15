library(tidyverse)
library(ggplot2)
library(lubridate)
library(raster)
library(colorspace)
library(scales)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# Read in the region-specific components
source("setup/regionDefinition.R")
source("setup/regional_analyses/regional.plot.layers.R")
AvisoDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"
eventDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
sstRateDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/monthly-rate"

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# Miscellaneous functions -------------------------------------------------

# Mean Kinetic Energy (MKE) function
mke.fun <- function(dat) {
  mke <- dat %>%
    dplyr::select(-ugosa, -vgosa) %>%
    dplyr::mutate(ugos = ugos * 100,
                  vgos = vgos * 100,
                  year = lubridate::year(fastDate(time)),
                  month = lubridate::month(fastDate(time))) %>%
    dplyr::group_by(lat, lon) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2)) %>%
    dplyr::ungroup()
  return(mke)
}

# Eddy Kinetic Energy (EKE) function
eke.fun <- function(dat) {
  eke <- dat %>%
    dplyr::mutate(ugos = ugos * 100,
                  vgos = vgos * 100,
                  ugosa = ugosa * 100,
                  vgosa = vgosa * 100,
                  year = lubridate::year(fastDate(time)),
                  month = lubridate::month(fastDate(time))) %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup()
  return(eke)
}

# Mean mean intensity
MeanInt <- function(events) {
  totMax <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = mean(intensity_mean, na.rm = TRUE)) %>%
    as_tibble()
}

# Total count
totalCntFun <- function(events) {
  freq <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(y = n())
  return(freq)
}

# aviso.dat <- "AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
# event.dat <- "AC-avhrr-only-v2.19810901-20180930_events.csv"
#
# ev <- as_tibble(fread(paste0(eventDir, "/", event.dat)))

# An über function to rule them all
matchEvent2KEgrid <- function(aviso.dat, event.dat, sstRate.dat) {
  require(data.table)
  # kinetic energy
  ke <- as_tibble(fread(paste0(AvisoDir, "/", aviso.dat), verbose = FALSE))
  colnames(ke) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")
  mke <- mke.fun(ke)
  eke <- eke.fun(ke)
  rm(ke)

  # mean of mean heat wave intensity
  ev <- as_tibble(fread(paste0(eventDir, "/", event.dat)))
  int <- MeanInt(ev)
  cnt <- totalCntFun(ev)
  rm(ev)

  # rate of change in SST
  trnd <- as_tibble(fread(paste0(sstRateDir, "/", sstRate.dat)))

  # make event intensity into a raster (interpolation grid)
  coordinates(int) <- ~ lon + lat
  gridded(int) <- TRUE
  int.grd <- raster(int)
  proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  crs(int.grd) <- proj

  # make count into a raster (interpolation grid)
  coordinates(cnt) <- ~ lon + lat
  gridded(cnt) <- TRUE
  cnt.grd <- raster(cnt)
  crs(cnt.grd) <- proj

  # make rate of change in SST into a raster
  coordinates(trnd) <- ~ lon + lat
  gridded(trnd) <- TRUE
  trnd.grd <- raster(trnd)
  crs(trnd.grd) <- proj

  # make eddy kinetic energy data into a raster
  coordinates(eke) <- ~ lon + lat
  gridded(eke) <- TRUE
  eke.grd <- raster(eke)
  crs(eke.grd) <- proj

  # make mean kinetic energy data into a raster
  coordinates(mke) <- ~ lon + lat
  gridded(mke) <- TRUE
  mke.grd <- raster(mke)
  crs(mke.grd) <- proj

  # regrid mke and eke to mean event intensity grid
  eke.grd <- resample(eke.grd, int.grd)
  mke.grd <- resample(mke.grd, int.grd)

  out <- tibble(y.int = as.vector(int.grd@data@values),
                       y.cnt = as.vector(cnt.grd@data@values),
                       y.trnd = as.vector(trnd.grd@data@values),
                       x.mke = as.vector(mke.grd@data@values),
                       x.eke = as.vector(eke.grd@data@values))
  out <- out[complete.cases(out), ]

  return(out)
}


# Prepare data ------------------------------------------------------------

AC.dat <- matchEvent2KEgrid(aviso.dat = "AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv",
                            event.dat = "AC-avhrr-only-v2.19810901-20180930_events.csv",
                            sstRate.dat = "AC-rate-19810901-20180901.csv")

BC.dat <- matchEvent2KEgrid(aviso.dat = "BC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv",
                            event.dat = "BC-avhrr-only-v2.19810901-20180930_events.csv",
                            sstRate.dat = "BC-rate-19810901-20180901.csv")

EAC.dat <- matchEvent2KEgrid(aviso.dat = "EAC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv",
                             event.dat = "EAC-avhrr-only-v2.19810901-20180930_events.csv",
                             sstRate.dat = "EAC-rate-19810901-20180901.csv")

GS.dat <- matchEvent2KEgrid(aviso.dat = "GS-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv",
                            event.dat = "GS-avhrr-only-v2.19810901-20180930_events.csv",
                            sstRate.dat = "GS-rate-19810901-20180901.csv")

KC.dat <- matchEvent2KEgrid(aviso.dat = "KC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv",
                            event.dat = "KC-avhrr-only-v2.19810901-20180930_events.csv",
                            sstRate.dat = "KC-rate-19810901-20180901.csv")


# 'Pattern' correlations --------------------------------------------------

cor(AC.dat$x.mke, AC.dat$y.int)
cor(AC.dat$x.eke, AC.dat$y.int)
cor(AC.dat$x.mke, AC.dat$y.cnt)
cor(AC.dat$x.eke, AC.dat$y.cnt)
cor(AC.dat$x.mke, AC.dat$y.trnd)
cor(AC.dat$x.eke, AC.dat$y.trnd)

cor(BC.dat$x.mke, BC.dat$y.int)
cor(BC.dat$x.eke, BC.dat$y.int)
cor(BC.dat$x.mke, BC.dat$y.cnt)
cor(BC.dat$x.eke, BC.dat$y.cnt)
cor(BC.dat$x.mke, BC.dat$y.trnd)
cor(BC.dat$x.eke, BC.dat$y.trnd)

cor(EAC.dat$x.mke, EAC.dat$y.int)
cor(EAC.dat$x.eke, EAC.dat$y.int)
cor(EAC.dat$x.mke, EAC.dat$y.cnt)
cor(EAC.dat$x.eke, EAC.dat$y.cnt)
cor(EAC.dat$x.mke, EAC.dat$y.trnd)
cor(EAC.dat$x.eke, EAC.dat$y.trnd)

cor(GS.dat$x.mke, GS.dat$y.int)
cor(GS.dat$x.eke, GS.dat$y.int)
cor(GS.dat$x.mke, GS.dat$y.cnt)
cor(GS.dat$x.eke, GS.dat$y.cnt)
cor(GS.dat$x.mke, GS.dat$y.trnd)
cor(GS.dat$x.eke, GS.dat$y.trnd)

cor(KC.dat$x.mke, KC.dat$y.int)
cor(KC.dat$x.eke, KC.dat$y.int)
cor(KC.dat$x.mke, KC.dat$y.cnt)
cor(KC.dat$x.eke, KC.dat$y.cnt)
cor(KC.dat$x.mke, KC.dat$y.trnd)
cor(KC.dat$x.eke, KC.dat$y.trnd)


# import shoreline
library(maptools)
gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"
AC.shore <- Rgshhs(paste0(gshhsDir, "/gshhs_l.b"),
                   xlim = c(bbox$AC[3], bbox$AC[4]),
                   ylim = c(bbox$AC[1], bbox$AC[2]), level = 1)$SP
BC.shore <- Rgshhs(paste0(gshhsDir, "/gshhs_l.b"),
                   xlim = c(bbox$BC[3], bbox$BC[4]),
                   ylim = c(bbox$BC[1], bbox$BC[2]), level = 1)$SP
EAC.shore <- Rgshhs(paste0(gshhsDir, "/gshhs_l.b"),
                   xlim = c(bbox$EAC[3], bbox$EAC[4]),
                   ylim = c(bbox$EAC[1], bbox$EAC[2]), level = 1)$SP
GS.shore <- Rgshhs(paste0(gshhsDir, "/gshhs_l.b"),
                   xlim = c(bbox$GS[3], bbox$GS[4]),
                   ylim = c(bbox$GS[1], bbox$GS[2]), level = 4)$SP
KC.shore <- Rgshhs(paste0(gshhsDir, "/gshhs_l.b"),
                   xlim = c(bbox$KC[3], bbox$KC[4]),
                   ylim = c(bbox$KC[1], bbox$KC[2]), level = 1)$SP

plot(AC.shore, axes = TRUE)
plot(BC.shore, axes = TRUE)
plot(EAC.shore, axes = TRUE)
plot(KC.shore, axes = TRUE)

cols1 <- sequential_hcl(99, "Plasma")
cols2 <- sequential_hcl(99, "YlGnBu")

ext <- extent(bbox$AC[3], bbox$AC[4], bbox$AC[1], bbox$AC[2])
ext <- as.vector(ext)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
sp::plot(AC.dat$int.grd, col = cols, ext = ext, interpolate = FALSE, )
plot(AC.shore, add = TRUE)
sp::plot(AC.dat$eke.grd, col = cols, ext = ext, interpolate = FALSE)
sp::plot(AC.dat$mke.grd, col = cols, ext = ext, interpolate = FALSE)
par(op)

library(rasterVis)
levelplot(AC.dat$int.grd, col.regions = cols2, colorkey = list(height = 1), contour = TRUE, margin = FALSE)
levelplot(AC.dat$eke.grd, col.regions = cols2, colorkey = list(height = 1), contour = TRUE, margin = FALSE)
levelplot(AC.dat$mke.grd, col.regions = cols2, colorkey = list(height = 1), contour = TRUE, margin = FALSE)

# plots
boundaries <- map('worldHires',
                  xlim = ext[1:2], ylim = ext[3:4],
                  plot = FALSE)
boundaries <- map2SpatialLines(boundaries, proj4string = CRS(projection(int.grd)))


library(sf)
levelplot(AC.dat$eke.grd, contour = TRUE)
levelplot(AC.dat$eke.grd, col.regions = cols, colorkey = list(height = 1), contour = TRUE, margin = FALSE) +
  layer(sp.lines(boundaries, lwd = 0.5))

levelplot(AC.dat$int.grd, layers = 1)

op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2))
plot(AC.dat$int.grd, col = cols, ext = ext, interpolate = FALSE,
     ylim = c(bbox$AC[1], bbox$AC[2]),
     xlim = c(bbox$AC[3], bbox$AC[4]))
plot(AC.shore, add = TRUE, col = "grey80")

par(op)
