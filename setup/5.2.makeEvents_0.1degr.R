# 5.2.makeEvents.R

library(heatwaveR)
library(fasttime)
library(plyr)
library(data.table)
library(doMC); doMC::registerDoMC(cores = 4)

# inDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/daily"
# evDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/MHWs"
# clDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/climatology"

clDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
evDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"

# DETECT FUNCTION using heatwaveR
OISST_detect2 <- function(dat) {
  require(heatwaveR)
  # dat$t <- fastDate(dat$t)
  out <- heatwaveR::detect_event(dat)
  event <- out$event
  return(event)
}

setClass('myDate')
setAs(from = "character", to = "myDate", def = function(from) as.Date(fastPOSIXct(from, tz = NULL)))

# Unlike ts2clm(), detect_event() works fine on large data sets and using multicore

# detect AC events
# parallel = FALSE gives more helpful error messages during testing...
AC.cl <- fread(paste0(clDir, "/AC-avhrr-only-v2.19810901-20180930_climatology.csv"),
                colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

AC.ev <- ddply(AC.cl, .(lon, lat), OISST_detect2, .parallel = TRUE, .progress = "text")
fwrite(AC.ev, file = paste0(evDir, "/AC-avhrr-only-v2.19810901-20180930_events.csv"))
remove(AC.cl); remove(AC.ev)


# detect BC events
BC.cl <- fread(paste0(clDir, "/BC-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

BC.ev <- ddply(BC.cl, .(lon, lat), OISST_detect2, .parallel = TRUE, .progress = "text")
fwrite(BC.ev, file = paste0(evDir, "/BC-avhrr-only-v2.19810901-20180930_events.csv"))
remove(BC.cl); remove(BC.ev)


# detect EAC events
EAC.cl <- fread(paste0(clDir, "/EAC-avhrr-only-v2.19810901-20180930_climatology.csv"),
                colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

EAC.ev <- ddply(EAC.cl, .(lon, lat), OISST_detect2, .parallel = TRUE, .progress = "text")
fwrite(EAC.ev, file = paste0(evDir, "/EAC-avhrr-only-v2.19810901-20180930_events.csv"))
remove(EAC.cl); remove(EAC.ev)


# detect GS events
GS.cl <- fread(paste0(clDir, "/GS-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

GS.ev <- ddply(GS.cl, .(lon, lat), OISST_detect2, .parallel = TRUE, .progress = "text")
fwrite(GS.ev, file = paste0(evDir, "/GS-avhrr-only-v2.19810901-20180930_events.csv"))
remove(GS.cl); remove(GS.ev)


# detect KC events
KC.cl <- fread(paste0(clDir, "/KC-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

KC.ev <- ddply(KC.cl, .(lon, lat), OISST_detect2, .parallel = TRUE, .progress = "text")
fwrite(KC.ev, file = paste0(evDir, "/KC-avhrr-only-v2.19810901-201809300_events.csv"))
remove(KC.cl); remove(KC.ev)
