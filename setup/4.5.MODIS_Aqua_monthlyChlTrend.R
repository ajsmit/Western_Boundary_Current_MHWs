# 4.5.MODIS_Aqua_monthlyChlTrendTrend.R

library(data.table)
library(tidyverse)
# The devtool version of ggplot2 has some nifty vector plotting options, that can be used to
# create very professional-looking vectors (see especially the) linejoin = "mitre" argument
# devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(ggpubr)


# some functions ----------------------------------------------------------

source("setup/functions.R")

# read in MODIS Aqua data -------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/regridded"
AC.chl <- fread(paste0(csvDir, "/AC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
BC.chl <- fread(paste0(csvDir, "/BC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
EAC.chl <- fread(paste0(csvDir, "/EAC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
KC.chl <- fread(paste0(csvDir, "/KC.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))
GS.chl <- fread(paste0(csvDir, "/GS.MODIS.Aqua.L3m_8D_CHL_chlor_a_9km_2002-07-05-2018-11-18.csv"))


# output Dir for chlor_a trend --------------------------------------------

outDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/monthly-rate"


# calc time series of monthly means ---------------------------------------

# dailyData <- AC.chl
daily2monthly <- function(dailyData) {
  monthlyData <- dailyData %>%
    dplyr::mutate(time = fastDate(time)) %>%
    dplyr::mutate(month = floor_date(time, unit = "month")) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(chlor_a = mean(chlor_a, na.rm = TRUE)) %>%
    dplyr::mutate(year = year(month)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(num = seq(1:length(chlor_a))) %>%
    dplyr::ungroup()
  return(monthlyData)
}

daily2yearly <- function(dailyData) {
  yearlyData <- dailyData %>%
    dplyr::mutate(time = fastDate(time)) %>%
    dplyr::mutate(year = floor_date(time, unit = "year")) %>%
    dplyr::group_by(lon, lat, year) %>%
    dplyr::summarise(chlor_a = mean(chlor_a, na.rm = TRUE)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(num = seq(1:length(chlor_a))) %>%
    dplyr::ungroup()
  return(yearlyData)
}

# function for calculating chl-a trend ------------------------------------

linFun <- function(data) {
  mod <- lm(chlor_a ~ num, data = data)
  trend <- data.frame(slope = summary(mod)$coefficients[2,1] * 10, pval = summary(mod)$coefficients[2,4])
  # trend <- summary(mod)$coefficients[2,c(1,4)]
  trend$pbracket <- cut(trend$pval, breaks = c(0, 0.001, 0.01, 0.05, 1))
  return(trend)
}

gls_fun <- function(df) {
  out <- tryCatch( {
    model <- gls(
      chlor_a ~ num,
      correlation = corARMA(form = ~ 1 | year, p = 2),
      method = "REML",
      data = df, na.action = na.exclude
    )
    stats <-
      data.frame(
        DT_model = round(as.numeric(coef(model)[2]) * 120, 3),
        se_trend = round(summary(model)[["tTable"]][[2, 2]], 10),
        sd_initial = round(sd(df$chlor_a, na.rm = TRUE), 2),
        sd_residual = round(sd(model$residuals), 2),
        r = NA,
        R2 = NA,
        p_trend = summary(model)[["tTable"]][[2, 4]],
        p_seas = NA,
        length = length(df$chlor_a),
        model = "gls"
      )
    return(stats)
  },
  error = function(cond) {
    stats <- data.frame(
      DT_model = NA,
      se_trend = NA,
      sd_initial = round(sd(df$chlor_a, na.rm = TRUE), 2),
      sd_residual = NA,
      r = NA,
      R2 = NA,
      p_trend = NA,
      p_seas = NA,
      length = length(df$chlor_a),
      model = "gls"
    )
    return(stats)
  })
  rm(model)
  return(out)
}

# calculate and save ------------------------------------------------------

AC_monthly <- daily2monthly(AC.chl)
AC_yearly <- daily2yearly(AC.chl)
AC_yearly <- AC_yearly[complete.cases(AC_yearly), ]
summary(AC_yearly)
system.time(AC_lm <- dlply(AC_yearly, .(lon, lat), .progress = "text",
                            .parallel = TRUE, gls_fun))
AC_df <- ldply(AC_gls, data.frame, .parallel = TRUE, .progress = "text")
fwrite(AC_df,
       paste0(outDir, "/AC-rate-",
              format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[1], "-",
              format(as.Date(AC_monthly$month, format = "%Y%m%d"), "%Y%m%d")[nrow(AC_monthly)], ".csv"))
rm(list = c("AC", "AC_monthly", "AC_gls", "AC_df"))
