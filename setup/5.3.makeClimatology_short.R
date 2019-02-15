# 5.3.makeClimatology_short.R

library(heatwaveR)
library(fasttime)
library(plyr)
library(data.table)
library(doMC); doMC::registerDoMC(cores = 4)
source("setup/functions.R")

setClass('myDate')
setAs("character", "myDate", function(from) as.Date(fastPOSIXct(from, tz = NULL)))


# Daily OISST data --------------------------------------------------------

OISSTDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
clDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"

ts2clm.grid <- function(data) {
  out <- ts2clm(data, x = t, y = temp,
                climatologyPeriod = c("1983-01-01", "2012-12-31"),
                robust = FALSE, maxPadLength = 3, windowHalfWidth = 5,
                pctile = 90, smoothPercentile = TRUE,
                smoothPercentileWidth = 31, clmOnly = TRUE)
  return(out)
}

# Calculate climatology and threshold etc.

# Agulhas Current
OISST_AC <- fread(paste0(OISSTDir, "/AC-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
AC_clim <- ddply(OISST_AC, .(lon, lat), ts2clm.grid,
                 .parallel = TRUE, .progress = "text")

fwrite(AC_clim,
       file = paste0(clDir, "/AC-avhrr-only-v2.19810901-20180930_climatology_short.csv"),
       append = FALSE)
remove(AC_clim); remove(OISST_AC)


# Brazil Current
OISST_BC <- fread(paste0(OISSTDir, "/BC-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
BC_clim <- ddply(OISST_BC, .(lon, lat), ts2clm.grid,
                 .parallel = TRUE, .progress = "text")

fwrite(BC_clim,
       file = paste0(clDir, "/BC-avhrr-only-v2.19810901-20180930_climatology_short.csv"),
       append = FALSE)
remove(BC_clim); remove(OISST_BC)


# East Australian Current
OISST_EAC <- fread(paste0(OISSTDir, "/EAC-avhrr-only-v2.19810901-20180930.csv"),
                   colClasses = list("numeric" = 1:3, "myDate" = 4),
                   col.names = c("lon", "lat", "temp", "t"))
EAC_clim <- ddply(OISST_EAC, .(lon, lat), ts2clm.grid,
                  .parallel = TRUE, .progress = "text")

fwrite(EAC_clim,
       file = paste0(clDir, "/EAC-avhrr-only-v2.19810901-20180930_climatology_short.csv"),
       append = FALSE)
remove(EAC_clim); remove(OISST_EAC)


# Gulf Stream
OISST_GS <- fread(paste0(OISSTDir, "/GS-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
GS_clim <- ddply(OISST_GS, .(lon, lat), ts2clm.grid,
                 .parallel = TRUE, .progress = "text")

fwrite(GS_clim,
       file = paste0(clDir, "/GS-avhrr-only-v2.19810901-20180930_climatology_short.csv"),
       append = FALSE)
remove(GS_clim); remove(OISST_GS)


# Kuroshio Current
OISST_KC <- fread(paste0(OISSTDir, "/KC-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
KC_clim <- ddply(OISST_KC, .(lon, lat), ts2clm.grid,
                 .parallel = TRUE, .progress = "text")

fwrite(KC_clim,
       file = paste0(clDir, "/KC-avhrr-only-v2.19810901-20180930_climatology_short.csv"),
       append = FALSE)
remove(KC_clim); remove(OISST_KC)
