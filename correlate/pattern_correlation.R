library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(colorspace)
library(scales)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# Read in the region-specific components
source("setup/regionDefinition.R")
source("regional_analyses/regional.plot.layers.R")
AvisoDir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/daily"
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

# region <- "EAC"

# An über function to rule them all
matchEvent2KEgrid <- function(region) {
  # define mask unions
  load(paste0("/Users/ajsmit/Dropbox/R/WBCs/masks/", region, "-mask_points.Rdata"))
  mke.ev.pts <- unique(rbind(as.data.frame(mask.pts$mke.pt), as.data.frame(mask.pts$ev.pt)))
  eke.ev.pts <- unique(rbind(as.data.frame(mask.pts$eke.pt), as.data.frame(mask.pts$ev.pt)))

  # load Aviso+, event, and SST trend data
  ke <- fread(paste0(AvisoDir, "/", region, "-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"))
  colnames(ke) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")
  ev <- fread(paste0(eventDir, "/", region, "-avhrr-only-v2.19810901-20180930_events.csv"))
  trnd <- fread(paste0(sstRateDir, "/", region, "-rate-19810901-20180901.csv"))
  load(paste0("masks/", region, "-mask_points.Rdata"))

  # calculate kinetic energy
  mke <- mke.fun(ke)
  eke <- eke.fun(ke)
  rm(ke)

  # calculate mean of mean heat wave intensity
  int <- MeanInt(ev)
  cnt <- totalCntFun(ev)
  rm(ev)

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

  int.df <- as.data.frame(int.grd, xy = TRUE, na.rm = FALSE, long = FALSE)
  cnt.df <- as.data.frame(cnt.grd, xy = TRUE, na.rm = FALSE, long = FALSE)
  trnd.df <- as.data.frame(trnd.grd, xy = TRUE, na.rm = FALSE, long = FALSE)
  mke.df <- as.data.frame(mke.grd, xy = TRUE, na.rm = FALSE, long = FALSE)
  eke.df <- as.data.frame(eke.grd, xy = TRUE, na.rm = FALSE, long = FALSE)

  dat.df <- data.frame(lon = mke.df$x,
                       lat = mke.df$y,
                       mke = mke.df$mke,
                       eke = eke.df$eke,
                       int = int.df[,3],
                       cnt = cnt.df[,3],
                       trnd = trnd.df$DT_model)

  mke.dat.df <- dat.df %>%
    dplyr::right_join(mke.ev.pts) %>%
    na.omit()
  eke.dat.df <- dat.df %>%
    dplyr::right_join(eke.ev.pts) %>%
    na.omit()

  out <- list(mke = mke.dat.df,
              eke = eke.dat.df)

  return(out)
}


# Prepare data ------------------------------------------------------------

AC.dat <- matchEvent2KEgrid("AC")
BC.dat <- matchEvent2KEgrid("BC")
EAC.dat <- matchEvent2KEgrid("EAC")
GS.dat <- matchEvent2KEgrid("GS")
KC.dat <- matchEvent2KEgrid("KC")


# 'Pattern' correlations --------------------------------------------------

cor(AC.dat$mke$mke, AC.dat$mke$int)
cor(AC.dat$mke$mke, AC.dat$mke$cnt)
cor(AC.dat$mke$mke, AC.dat$mke$trnd)
cor(AC.dat$eke$eke, AC.dat$eke$int)
cor(AC.dat$eke$eke, AC.dat$eke$cnt)
cor(AC.dat$eke$eke, AC.dat$eke$trnd)

cor(BC.dat$mke$mke, BC.dat$mke$int)
cor(BC.dat$mke$mke, BC.dat$mke$cnt)
cor(BC.dat$mke$mke, BC.dat$mke$trnd)
cor(BC.dat$eke$eke, BC.dat$eke$int)
cor(BC.dat$eke$eke, BC.dat$eke$cnt)
cor(BC.dat$eke$eke, BC.dat$eke$trnd)

cor(EAC.dat$mke$mke, EAC.dat$mke$int)
cor(EAC.dat$mke$mke, EAC.dat$mke$cnt)
cor(EAC.dat$mke$mke, EAC.dat$mke$trnd)
cor(EAC.dat$eke$eke, EAC.dat$eke$int)
cor(EAC.dat$eke$eke, EAC.dat$eke$cnt)
cor(EAC.dat$eke$eke, EAC.dat$eke$trnd)

cor(GS.dat$mke$mke, GS.dat$mke$int)
cor(GS.dat$mke$mke, GS.dat$mke$cnt)
cor(GS.dat$mke$mke, GS.dat$mke$trnd)
cor(GS.dat$eke$eke, GS.dat$eke$int)
cor(GS.dat$eke$eke, GS.dat$eke$cnt)
cor(GS.dat$eke$eke, GS.dat$eke$trnd)

cor(KC.dat$mke$mke, KC.dat$mke$int)
cor(KC.dat$mke$mke, KC.dat$mke$cnt)
cor(KC.dat$mke$mke, KC.dat$mke$trnd)
cor(KC.dat$eke$eke, KC.dat$eke$int)
cor(KC.dat$eke$eke, KC.dat$eke$cnt)
cor(KC.dat$eke$eke, KC.dat$eke$trnd)

