# Create a time series of KE and MHW intensity for a small region within each
# of the areas of maximal MWH activity

library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(astrochron)
library(WaveletComp)
library(heatwaveR)
library(grid)
library(doMC); doMC::registerDoMC(cores = 4)


# Setup -------------------------------------------------------------------

source("/Users/ajsmit/Dropbox/R/papers/diatom_age/custom_theme.R")
source("/Users/ajsmit/Dropbox/R/WBCs/setup/functions.R")

region <- "EAC"

# Bounding boxes
bbox <- data.frame("AC" = c(-42.5, -40, 12.5, 15),
                   "BC" = c(-40, -37.5, 305, 307.5),
                   "EAC" = c( -37.5, -35, 152.5, 155),
                   "GS" = c(40, 42.5, 295, 297.5),
                   "KC" = c(40, 42.5, 145, 147.5),
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))


# Function to make time series of bbox region -----------------------------

make_ts <- function(region) {
  # The Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged" # /AC, /BA, etc.
  aviso.dat <- "outfile_regridded.nc"
  coords <- bbox[, region]
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  ugos <- ncvar_get(nc, varid = "ugos",
                    start = c(LonIdx[1], LatIdx[1], 1),
                    count = c(length(LonIdx), length(LatIdx), length(nc$dim$time$vals))) * 100 %>%
    round(4)

  dimnames(ugos) <- list(lon = nc$dim$lon$vals[LonIdx],
                         lat = nc$dim$lat$vals[LatIdx],
                         date = nc$dim$time$vals)
  ke <- as.data.table(melt(ugos, value.name = "ugos"))
  ke[, c("date", "vgos", "ugosa", "vgosa") := .( # slow
    as.Date(date, origin = "1950-01-01 00:00:00"),
    as.vector(ncvar_get(nc, varid = "vgos",
                        start = c(LonIdx[1], LatIdx[1], 1),
                        count = c(length(LonIdx), length(LatIdx), length(nc$dim$time$vals)))) * 100,
    as.vector(ncvar_get(nc, varid = "ugosa",
                        start = c(LonIdx[1], LatIdx[1], 1),
                        count = c(length(LonIdx), length(LatIdx), length(nc$dim$time$vals)))) * 100,
    as.vector(ncvar_get(nc, varid = "vgosa",
                        start = c(LonIdx[1], LatIdx[1], 1),
                        count = c(length(LonIdx), length(LatIdx), length(nc$dim$time$vals)))) * 100
  )
  ]
  nc_close(nc)

  ke <- ke[, lapply(.SD, mean), by = c("date")]
  ke[, c("lon", "lat") := NULL]

  ke[, eke := (0.5 * ((vgosa)^2 + (ugosa)^2))]
  ke[, c("mke_s", "eke_s") := .(
    0.5 * (frollmean(vgos, n = 30, align = "center", na.rm = TRUE)^2 + # frollmean is a bloody fast data.table function!
             frollmean(ugos, n = 30, align = "center", na.rm = TRUE)^2),
    frollmean(eke, n = 30, align = "center", na.rm = TRUE)
  )]
  # ke <- na.omit(ke)
  # ke[, c("ugos", "vgos", "ugosa", "vgosa") := NULL]

  # The OISST data
  clim.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
  clim.dat <- "-avhrr-only-v2.19810901-20180930_climatology.csv"
  clim <- fread(paste0(clim.dir, "/", region, clim.dat))
  colnames(clim)[colnames(clim) == "t"] <- "date"
  clim[, date := fastDate(date)]
  clim <- clim[lon >= bbox["lonmin", region] & lon <= bbox["lonmax", region], ]
  clim <- clim[lat >= bbox["latmin", region] & lat <= bbox["latmax", region], ]
  clim <- clim[, lapply(.SD, mean), by = c("date")]
  clim[, c("ex_s", "ex") := .(frollmean((temp - thresh), n = 30, align = "center", na.rm = TRUE),
                             temp - thresh)]
  clim[, c("lon", "lat") := NULL]

  setkey(ke, date)
  setkey(clim, date)
  dat <- clim[ke, nomatch = 0]
  dat <- na.omit(dat)

  ts.dir <- "/Volumes/Benguela/spatial/processed/WBC/misc_results"
  fwrite(dat, file = paste0(ts.dir, "/", region, "-wavelet-ts.csv"), showProgress = TRUE, nThread = 6)

  return(dat)
}


# Apply function to create and save ts ------------------------------------

AC.ts <- make_ts("AC")
BC.ts <- make_ts("BC")
EAC.ts <- make_ts("EAC")
GS.ts <- make_ts("GS")
KC.ts <- make_ts("KC")


# Load data and do wavelet analysis ---------------------------------------

ts.dir <- "/Volumes/Benguela/spatial/processed/WBC/misc_results"
AC.ts <- fread(file = paste0(ts.dir, "/", "AC", "-wavelet-ts.csv"))
BC.ts <- fread(file = paste0(ts.dir, "/", "BC", "-wavelet-ts.csv"))
EAC.ts <- fread(file = paste0(ts.dir, "/", "EAC", "-wavelet-ts.csv"))
GS.ts <- fread(file = paste0(ts.dir, "/", "GS", "-wavelet-ts.csv"))
KC.ts <- fread(file = paste0(ts.dir, "/", "KC", "-wavelet-ts.csv"))

# check for autocorreltion using 'auto.arima()' in the 'forecast' package...
# library(forecast)
#
# auto.arima(AC.ts$eke_roll2, max.p = 3, max.q = 3, stationary = FALSE,
#            seasonal = FALSE)

# apply pre-whitening to the data; this effectively removes the above
# autocorrelation structure and the residuals are then used for the remainder
# of the analyses; this allows us to easily identify the embedded spectral
# frequencies


# ts <- EAC.ts
# region <- "EAC"

wl_plts <- function(ts, region) {
  ts <- as.data.frame(ts)
  ts$n <- c(1:nrow(ts))

  ev <- detect_event(data = ts[, c("date", "temp", "seas", "thresh")], x = date)
  ev.n <- ev$event %>%
    dplyr::arrange(desc(intensity_max_relThresh)) %>%
    dplyr::select(date_start) %>%
    head(10) %>%
    c()
  tlabs <- c(ev.n$date_start)

  # first the exceedence data
  ts.ex <- prewhiteAR(ts[,c("n", "ex_s")], order = 3, method = "ols", aic = TRUE,
                      genplot = FALSE, verbose = FALSE)
  colnames(ts.ex) <- c("seq", "ex_s")
  row.names(ts.ex) <- ts$date[4:nrow(ts)]

  # then the eke data
  ts.eke <- prewhiteAR(ts[,c("n", "eke_s")], order = 3, method = "ols", aic = TRUE,
                       genplot = FALSE, verbose = FALSE)
  colnames(ts.eke) <- c("seq", "eke_s")
  row.names(ts.eke) <- ts$date[4:nrow(ts)]

  ts.w <- ts.ex
  ts.w$eke_s <- ts.eke[, 2]

  # coherence - whitened time series
  wl.w <- analyze.coherency(ts.w, c("ex_s", "eke_s"),
                             loess.span = 0,
                             dt = 1, dj = 1/50,
                             # window.size.t = 1, window.size.s = 1/2,
                             lowerPeriod = 6,
                             make.pval = TRUE, n.sim = 30,
                             date.format = "%Y-%m-%d")


  # wl.ex <- analyze.wavelet(ts.ex, my.series = "ex_s", loess.span = 0, dt = 1,
  #                           dj = 1/50, lowerPeriod = 6, make.pval = TRUE, n.sim = 30,
  #                           method = "Fourier.rand", date.format = "%Y-%m-%d", verbose = FALSE)
  # wl.eke <- analyze.wavelet(ts.eke, "eke_s", loess.span = 0, dt = 1,
  #                            dj = 1/50, lowerPeriod = 6, make.pval = TRUE, n.sim = 30,
  #                            method = "Fourier.rand", date.format = "%Y-%m-%d", verbose = FALSE)

  # coherence plots
  wl.dir <- "/Users/ajsmit/Dropbox/R/WBCs/figures"
  pdf(file = paste0(wl.dir, "/", region, "wavelets_non-whitened.pdf"), width = 6 * 1.33, height = 5.25 * 1.33,
      pointsize = 7.5)
  par(mfrow = c(3, 1),
      mar = c(4, 4.5, 3, 1))
  wc.image(wl.w, main = "Cross-wavelet power spectrum, exceedence over EKE",
           legend.params = list(lab = "Cross-wavelet power level"),
           periodlab = "Period (days)", show.date = TRUE, graphics.reset = FALSE,
           which.arrow.sig = "wt", color.key = "interval",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)

  wt.image(wl.w, my.series = "ex_s", siglvl = 0.05, col.contour = "black", color.key = "quantile",
           legend.params = list(lab = "Wavelet power level", label.digits = 2, shrink = 1.0),
           show.date = TRUE, date.format = "%Y-%m-%d",
           timelab = "Year", periodlab = "Period (days)", lwd = 1, graphics.reset = FALSE,
           main = "Exceedence",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)
  wt.image(wl.w, my.series = "eke_s", siglvl = 0.05, col.contour = "black", color.key = "quantile",
           legend.params = list(lab = "Wavelet power level", label.digits = 2, shrink = 1.0),
           show.date = TRUE, date.format = "%Y-%m-%d",
           timelab = "Year", periodlab = "Period (days)", lwd = 1, main = "Kinetic energy",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)
  dev.off()

  # coherence - non-whitened time series
  wl.nl <- analyze.coherency(ts, c("ex_s", "eke_s"),
                             loess.span = 0,
                             dt = 1, dj = 1/50,
                             # window.size.t = 1, window.size.s = 1/2,
                             lowerPeriod = 6,
                             make.pval = TRUE, n.sim = 30,
                             date.format = "%Y-%m-%d")

  pdf(file = paste0(wl.dir, "/", region, "-wavelets_non-whitened.pdf"), width = 6 * 1.33, height = 5.25 * 1.33,
      pointsize = 7.5)
  par(mfrow = c(3, 1),
      mar = c(4, 4.5, 3, 1))
  wc.image(wl.nl, main = "Cross-wavelet power spectrum, exceedence over EKE",
           legend.params = list(lab = "Cross-wavelet power level"),
           periodlab = "Period (days)", show.date = TRUE, graphics.reset = FALSE,
           which.arrow.sig = "wt", color.key = "interval",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)

  wt.image(wl.nl, my.series = "ex_s", siglvl = 0.05, col.contour = "black", color.key = "quantile",
           legend.params = list(lab = "Wavelet power level", label.digits = 2, shrink = 1.0),
           show.date = TRUE, date.format = "%Y-%m-%d",
           timelab = "Year", periodlab = "Period (days)", lwd = 1, graphics.reset = FALSE,
           main = "Exceedence",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)

  wt.image(wl.nl, my.series = "eke_s", siglvl = 0.05, col.contour = "black", color.key = "quantile",
           legend.params = list(lab = "Wavelet power level", label.digits = 2, shrink = 1.0),
           show.date = TRUE, date.format = "%Y-%m-%d",
           timelab = "Year", periodlab = "Period (days)", lwd = 1, main = "Kinetic energy",
           spec.time.axis = list(at = tlabs,
                                 labels = c("1", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10")),
           timetcl = NULL, timetck = 1)
  dev.off()

}

wl_plts(AC.ts, "AC")
wl_plts(BC.ts, "BC")
wl_plts(EAC.ts, "EAC")
wl_plts(GS.ts, "GS")
wl_plts(KC.ts, "KC")


# Detect heat waves -------------------------------------------------------

library(heatwaveR)

AC.ev <- detect_event(data = AC.ts[, c("date", "temp", "seas", "thresh")], x = date)
BC.ev <- detect_event(data = BC.ts[, c("date", "temp", "seas", "thresh")], x = date)
EAC.ev <- detect_event(data = EAC.ts[, c("date", "temp", "seas", "thresh")], x = date)
GS.ev <- detect_event(data = GS.ts[, c("date", "temp", "seas", "thresh")], x = date)
KC.ev <- detect_event(data = KC.ts[, c("date", "temp", "seas", "thresh")], x = date)


c(top.n[, "date_start"])
