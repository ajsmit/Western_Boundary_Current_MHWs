# 5.1.makeClimatology_0.1degr.R

library(heatwaveR)
library(tictoc)
library(data.table)
library(doMC); doMC::registerDoMC(cores = 3)
source("setup/functions.R")


# Daily OISST data --------------------------------------------------------

OISSTDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
clDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"

ts2clm.grid <- function(data) {
  out <- ts2clm(data, x = t, y = temp,
                climatologyPeriod = c("1983-01-01", "2012-12-31"),
                robust = FALSE, maxPadLength = 3, windowHalfWidth = 5,
                pctile = 90, smoothPercentile = TRUE,
                smoothPercentileWidth = 31, clmOnly = FALSE)
  return(out)
}

# Calculate climatology and threshold etc.

# Agulhas Current
tic()
OISST_AC <- vroom::vroom(file = paste0(OISSTDir, "/AC-avhrr-only-v2.19810901-20180930.csv"),
                         delim = ",",
                         col_names = c("lon", "lat", "temp", "t"),
                         col_types  = c(lon = "d", lat = "d", temp = "d", t = "D"),
                         na = c("", "NA", "NULL"))
toc()

# Aggregate to 0.1 degree pixels and save again
tic()
OISST_AC <-  plyr::ddply(OISST_AC,
                         .(round(OISST_AC$lon, 1), round(OISST_AC$lat, 1), t),
                         summarize,
                         temp = mean(temp),
                         .progress = "text")
colnames(OISST_AC) <- c("lon", "lat", "temp")
toc()

vroom::vroom_write(OISST_AC,
                   path = paste0(OISSTDir, "/AC-avhrr-only-v2.19810901-20180930_0.1degr.csv"),
                   delim = ",")

AC_clim <- ddply(OISST_AC,
                 .(lon, lat),
                 ts2clm.grid,
                 .progress = "text")

vroom::vroom_write(AC_clim,
                   path = paste0(clDir, "/AC-avhrr-only-v2.19810901-20180930_climatology_0.1degr.csv"),
                   delim = ",")
remove(AC_clim); remove(OISST_AC)


# Brazil Current
# Works with parallel = TRUE
OISST_BC <- fread(paste0(OISSTDir, "/BC-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))

# Aggregate to 0.1 degree pixels and save again
OISST_BC$lon <- round(OISST_BC$lon, 1)
OISST_BC$lat <- round(OISST_BC$lat, 1)
OISST_BC <-  ddply(OISST_BC, .(lon, lat), summarize,
                   temp = mean(temp))

fwrite(OISST_BC,
       file = paste0(OISSTDir, "/BC-avhrr-only-v2.19810901-20180930_0.1degr.csv"),
       append = FALSE)

BC_clim <- ddply(OISST_BC, .(lon, lat), ts2clm.grid,
                 .parallel = FALSE, .progress = "text")

fwrite(BC_clim,
       file = paste0(clDir, "/BC-avhrr-only-v2.19810901-20180930_climatology_0.1degr.csv"),
       append = FALSE)
remove(BC_clim); remove(OISST_BC)


# East Australian Current
tic()
OISST_EAC <- fread(paste0(OISSTDir, "/KC-avhrr-only-v2.19810901-20180930.csv"),
                   colClasses = list("numeric" = 1:3, "myDate" = 4),
                   col.names = c("lon", "lat", "temp", "t"))
toc()

tic()
OISST_EAC[, c("lon", "lat") := .(round(lon, 1), round(lat, 1))]
toc()

fwrite(OISST_EAC, file = paste0(OISSTDir, "/EAC-avhrr-only-v2.19810901-20180930_0.1degr.csv"))

EAC_clim <- OISST_EAC[, .()]

EAC_clim <- ddply(OISST_EAC,
                  .(lon, lat),
                  ts2clm.grid,
                  .progress = "text")

fwrite(EAC_clim, file = paste0(clDir, "/EAC-avhrr-only-v2.19810901-20180930_climatology_0.1degr.csv"))
remove(EAC_clim); remove(OISST_EAC)


# Gulf Stream
# Will only run if parallel = FALSE
OISST_GS <- fread(paste0(OISSTDir, "/GS-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
GS_clim <- ddply(OISST_GS, .(lon, lat), ts2clm.grid,
                 .parallel = FALSE, .progress = "text")

fwrite(GS_clim,
       file = paste0(clDir, "/GS-avhrr-only-v2.19810901-20180930_climatology.csv"),
       append = FALSE)
remove(GS_clim); remove(OISST_GS)


# Kuroshio Current
# Will only run if parallel = FALSE
OISST_KC <- fread(paste0(OISSTDir, "/KC-avhrr-only-v2.19810901-20180930.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4),
                  col.names = c("lon", "lat", "temp", "t"))
KC_clim <- ddply(OISST_KC, .(lon, lat), ts2clm.grid,
                 .parallel = FALSE, .progress = "text")

fwrite(KC_clim,
       file = paste0(clDir, "/KC-avhrr-only-v2.19810901-20180930_climatology.csv"),
       append = FALSE)
remove(KC_clim); remove(OISST_KC)
