# 4.2.Aviso+_monthlyVelocityTrends.R

library(doMC); doMC::registerDoMC(cores = 4) # multi-core
require(mgcv)
library(readr)
library(plyr)
library(dplyr)
library(lubridate)
library(tidyverse)
library(data.table) # because it reads in csv data faster...


# Some functions ----------------------------------------------------------

# Source the regression function
source("setup/functions.R")

# a function to create monthly data from daily data
daily2monthly <- function(dailyData) {
  monthlyData <- dailyData %>%
    dplyr::mutate(time = fastDate(time)) %>%
    dplyr::mutate(month = floor_date(time, unit = "month")) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(u = mean(u, na.rm = TRUE),
                     v = mean(v, na.rm = TRUE),
                     sla = mean(sla, na.rm = TRUE),
                     adt = mean(adt, na.rm = TRUE)) %>%
    dplyr::mutate(year = year(month)) %>%
    dplyr::mutate(velocity = sqrt((v^2) + (u^2)),
                  dir = atan2(u, v) * 180/pi,
                  eke = 0.5 * (u^2 + v^2),
                  arrow_size = ((abs(u * v) / max(abs(u * v))) + 0.3) / 6) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(num = seq(1:length(u))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(month = month)
  return(monthlyData)
}


# the GLS-AR2 -------------------------------------------------------------

gls_fun <- function(df) {
  out <- tryCatch( {
    model <- gls(
      velocity ~ num,
      correlation = corARMA(form = ~ 1 | year, p = 2),
      method = "REML",
      data = df, na.action = na.exclude
    )
    stats <-
      data.frame(
        DT_model = round(as.numeric(coef(model)[2]) * 120, 3),
        se_trend = round(summary(model)[["tTable"]][[2, 2]], 10),
        sd_initial = round(sd(df$velocity, na.rm = TRUE), 2),
        sd_residual = round(sd(model$residuals), 2),
        r = NA,
        R2 = NA,
        p_trend = summary(model)[["tTable"]][[2, 4]],
        p_seas = NA,
        length = length(df$velocity),
        model = "gls"
      )
    return(stats)
  },
  error = function(cond) {
    stats <- data.frame(
      DT_model = NA,
      se_trend = NA,
      sd_initial = NA,
      sd_residual = NA,
      r = NA,
      R2 = NA,
      p_trend = NA,
      p_seas = NA,
      length = length(df$velocity),
      model = "gls"
    )
    return(stats)
  })
  rm(model)
  return(out)
}


# Calculate velocity trends -----------------------------------------------

# 1. calculate monthly means from dailys;
# 2. apply the gls;
# 3. convert the resultant list to a data frame of model params;
# 4. save the data as a csv file
# (sorry, a bit repetitive...)

inDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"
outDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/monthly-rate"
AC <- as_tibble(fread(paste0(inDir, "/AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                      colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                      col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
AC_monthly <- daily2monthly(AC)
system.time(AC_gls <- dlply(AC_monthly, .(lon, lat), .progress = "text", .parallel = TRUE, gls_fun))
AC_df <- ldply(AC_gls, data.frame, .progress = "text")
write_csv(AC_df,
          paste0(outDir, "/AC-velocity-rate-",
                 format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
                 format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(AC_monthly)], ".csv"))
rm(list = c("AC", "AC_df", "AC_gls", "AC_monthly"))

BC <- as_tibble(fread(paste0(inDir, "/BC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                      colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                      col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
BC_monthly <- daily2monthly(BC)
system.time(BC_gls <- dlply(BC_monthly, .(lon, lat), .progress = "text", .parallel = TRUE, gls_fun))
BC_df <- ldply(BC_gls, data.frame, .progress = "text")
write_csv(BC_df,
          paste0(outDir, "/BC-velocity-rate-",
                 format(as.Date(BC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
                 format(as.Date(BC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(BC_monthly)], ".csv"))
rm(list = c("BC", "BC_df", "BC_gls", "BC_monthly"))

EAC <- as_tibble(fread(paste0(inDir, "/EAC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                       colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                       col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
system.time(EAC_monthly <- daily2monthly(EAC))
system.time(EAC_gls <- dlply(EAC_monthly, .(lon, lat), .progress = "text",
                             .parallel = TRUE, gls_fun))
EAC_df <- ldply(EAC_gls, data.frame, .progress = "text")
write_csv(EAC_df,
          paste0(outDir, "/EAC-velocity-rate-",
                 format(as.Date(EAC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
                 format(as.Date(EAC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(EAC_monthly)], ".csv"))
rm(list = c("EAC", "EAC_df", "EAC_gls", "EAC_monthly"))

GS <- as_tibble(fread(paste0(inDir, "/GS-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                      colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                      col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
GS_monthly <- daily2monthly(GS)
system.time(GS_gls <- dlply(GS_monthly, .(lon, lat), .progress = "text", .parallel = TRUE, gls_fun))
GS_df <- ldply(GS_gls, data.frame, .progress = "text")
write_csv(GS_df,
          paste0(outDir, "/GS-velocity-rate-",
                 format(as.Date(GS_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
                 format(as.Date(GS_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(GS_monthly)], ".csv"))
rm(list = c("GS", "GS_df", "GS_gls", "GS_monthly"))

KC <- as_tibble(fread(paste0(inDir, "/KC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                      colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                      col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
KC_monthly <- daily2monthly(KC)
system.time(KC_gls <- dlply(KC_monthly, .(lon, lat), .progress = "text", .parallel = TRUE, gls_fun))
KC_df <- ldply(KC_gls, data.frame, .progress = "text")
write_csv(KC_df,
          paste0(outDir, "/KC-velocity-rate-",
                 format(as.Date(KC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
                 format(as.Date(KC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(KC_monthly)], ".csv"))
rm(list = c("KC", "KC_df", "KC_gls", "KC_monthly"))
