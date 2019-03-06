# 5.3.makeProtoEvents.R

library(tidyverse)
library(heatwaveR)
library(fasttime)
library(data.table)
library(doMC); doMC::registerDoMC(cores = 3)
# library(furrr) # for future_map()
# plan(multiprocess)


# inDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/daily"
# evDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/MHWs"
# clDir <- "/Users/ajsmit/spatial/processed/OISSTv2/WBC/climatology"

clDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
evDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"

# DETECT FUNCTION using heatwaveR
# OISST_proto_detect2 <- function(x) { heatwaveR::detect_event(x, protoEvents = TRUE) }

setClass('myDate')
setAs(from = "character", to = "myDate", def = function(from) as.Date(fastPOSIXct(from, tz = NULL)))

# Unlike ts2clm(), detect_event() works fine on large data sets and using multicore

OISST_proto_detect <- function(x) { heatwaveR::detect_event(x, protoEvents = TRUE) }

# detect AC events
# parallel = FALSE gives more helpful error messages during testing...
AC.cl <- fread(paste0(clDir, "/AC-avhrr-only-v2.19810901-20180930_climatology.csv"),
                colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

AC.sub <- AC.cl[lon <= 15 & lat <= -43, ]
#
# # try with future_map() in 'furrr' package
# AC.ev <- AC.sub %>%
#   dplyr::group_by(lon, lat) %>%
#   nest() %>%
#   dplyr::mutate(ev = future_map(data, heatwaveR::detect_event, protoEvents = TRUE,
#                                 .progress
#  = TRUE)) %>%
#   dplyr::select(-data) %>%
#   unnest(ev)

AC.ev <- AC.cl %>%
  dplyr::group_by(lon, lat) %>%
  nest() %>%
  dplyr::mutate(ev = map(data, heatwaveR::detect_event, protoEvents = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(ev)

fwrite(AC.ev, file = paste0(evDir, "/AC-avhrr-only-v2.19810901-20180930_proto_events.csv"))
remove(AC.cl); remove(AC.ev)


# detect BC events
BC.cl <- fread(paste0(clDir, "/BC-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

BC.ev <- BC.cl %>%
  dplyr::group_by(lon, lat) %>%
  nest() %>%
  dplyr::mutate(ev = map(data, heatwaveR::detect_event, protoEvents = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(ev)

fwrite(BC.ev, file = paste0(evDir, "/BC-avhrr-only-v2.19810901-20180930_proto_events.csv"))
remove(BC.cl); remove(BC.ev)


# detect EAC events
EAC.cl <- fread(paste0(clDir, "/EAC-avhrr-only-v2.19810901-20180930_climatology.csv"),
                colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

EAC.ev <- EAC.cl %>%
  dplyr::group_by(lon, lat) %>%
  nest() %>%
  dplyr::mutate(ev = map(data, heatwaveR::detect_event, protoEvents = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(ev)

fwrite(EAC.ev, file = paste0(evDir, "/EAC-avhrr-only-v2.19810901-20180930_proto_events.csv"))
remove(EAC.cl); remove(EAC.ev)


# detect GS events
GS.cl <- fread(paste0(clDir, "/GS-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

GS.ev <- GS.cl %>%
  dplyr::group_by(lon, lat) %>%
  nest() %>%
  dplyr::mutate(ev = map(data, heatwaveR::detect_event, protoEvents = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(ev)

fwrite(GS.ev, file = paste0(evDir, "/GS-avhrr-only-v2.19810901-20180930_proto_events.csv"))
remove(GS.cl); remove(GS.ev)


# detect KC events
KC.cl <- fread(paste0(clDir, "/KC-avhrr-only-v2.19810901-20180930_climatology.csv"),
               colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

KC.ev <- KC.cl %>%
  dplyr::group_by(lon, lat) %>%
  nest() %>%
  dplyr::mutate(ev = map(data, heatwaveR::detect_event, protoEvents = TRUE)) %>%
  dplyr::select(-data) %>%
  unnest(ev)

fwrite(KC.ev, file = paste0(evDir, "/KC-avhrr-only-v2.19810901-20180930_proto_events.csv"))
remove(KC.cl); remove(KC.ev)
