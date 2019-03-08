# 4.1.monthlyMeanTrends.R

library(doMC); doMC::registerDoMC(cores = 4) # multi-core
require(mgcv)
library(data.table)
library(readr)
library(dplyr)
library(plyr)
library(lubridate)
library(tidyverse)


# Some functions ----------------------------------------------------------

# Source the regression function
source("setup/functions.R")

# a function to create monthly data from daily data
daily2monthly <- function(dailyData) {
  monthlyData <- dailyData %>%
    dplyr::mutate(t = fastDate(t)) %>%
    dplyr::mutate(month = floor_date(t, unit = "month")) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
    dplyr::mutate(year = year(month)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(num = seq(1:length(temp))) %>%
    dplyr::ungroup()
  return(monthlyData)
}


# Make monthly OISST ------------------------------------------------------

# setAs("character", "myDate", function(from) as.Date(fastPOSIXct(from, tz = NULL)))
# setClass('myDate')

inDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
outDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/monthly-rate"

# 1. calculate monthly means from dailys;
# 2. apply the gls;
# 3. convert the resultant list to a data frame of model params;
# 4. save the data as a csv file
# (sorry, a bit repetitive...)

# the Agulhas Current
AC <- fread(paste0(inDir, "/AC-avhrr-only-v2.19810901-20180930.csv"),
            colClasses = list("numeric" = 1:3, "Date" = 4),
            col.names = c("lon", "lat", "temp", "t"))
# AC$t <- fastDate(AC$t) # faster here, not in the daily2monthly function...
AC_monthly <- daily2monthly(AC)
system.time(AC_gls <- dlply(AC_monthly, .(lon, lat), .progress = "text",
                            .parallel = TRUE, gls_fun))
AC_df <- ldply(AC_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(AC_df,
       paste0(outDir, "/AC-rate-",
              format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(AC_monthly)], ".csv"))
rm(list = c("AC", "AC_monthly", "AC_gls", "AC_df"))

# the Brazil Current
BC <- fread(paste0(inDir, "/BC-avhrr-only-v2.19810901-20180930.csv"),
            colClasses = list("numeric" = 1:3, "Date" = 4),
            col.names = c("lon", "lat", "temp", "t"))
BC_monthly <- daily2monthly(BC)
system.time(BC_gls <- dlply(BC_monthly, .(lon, lat), .progress = "text",
                            .parallel = TRUE, gls_fun))
BC_df <- ldply(BC_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(BC_df,
       paste0(outDir, "/BC-rate-",
              format(as.Date(BC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(BC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(BC_monthly)], ".csv"))
rm(list = c("BC", "BC_monthly", "BC_gls", "BC_df"))

# the East Australian Current
EAC <- fread(paste0(inDir, "/EAC-avhrr-only-v2.19810901-20180930.csv"),
             colClasses = list("numeric" = 1:3, "Date" = 4),
             col.names = c("lon", "lat", "temp", "t"))
EAC_monthly <- daily2monthly(EAC)
system.time(EAC_gls <- dlply(EAC_monthly, .(lon, lat), .progress = "text",
                             .parallel = TRUE, gls_fun))
EAC_df <- ldply(EAC_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(EAC_df,
       paste0(outDir, "/EAC-rate-",
              format(as.Date(EAC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(EAC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(EAC_monthly)], ".csv"))
rm(list = c("EAC", "EAC_monthly", "EAC_gls", "EAC_df"))

# the Gulf Stream
GS <- fread(paste0(inDir, "/GS-avhrr-only-v2.19810901-20180930.csv"),
            colClasses = list("numeric" = 1:3, "Date" = 4),
            col.names = c("lon", "lat", "temp", "t"))
GS_monthly <- daily2monthly(GS)
system.time(GS_gls <- dlply(GS_monthly, .(lon, lat), .progress = "text",
                            .parallel = TRUE, gls_fun))
GS_df <- ldply(GS_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(GS_df,
       paste0(outDir, "/GS-rate-",
              format(as.Date(GS_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(GS_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(GS_monthly)], ".csv"))
rm(list = c("GS", "GS_monthly", "GS_gls", "GS_df"))

# the Kuroshio Current
KC <- fread(paste0(inDir, "/KC-avhrr-only-v2.19810901-20180930.csv"),
            colClasses = list("numeric" = 1:3, "Date" = 4),
            col.names = c("lon", "lat", "temp", "t"))
KC_monthly <- daily2monthly(KC)
system.time(KC_gls <- dlply(KC_monthly, .(lon, lat), .progress = "text",
                            .parallel = TRUE, gls_fun))
KC_df <- ldply(KC_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(KC_df,
       paste0(outDir, "/KC-rate-",
              format(as.Date(KC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(KC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(KC_monthly)], ".csv"))
rm(list = c("KC", "KC_monthly", "KC_gls", "KC_df"))


# Stack data frames and plot as facet_wrap...
