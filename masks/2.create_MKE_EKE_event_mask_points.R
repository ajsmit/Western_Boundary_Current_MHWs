library(tidyverse)
library(ncdf4)
library(RcppRoll)
library(lubridate)
library(data.table)
library(colorspace)
library(raster)


# Setup -------------------------------------------------------------------

source("setup/regionDefinition.R")
# for converting julian date to calendar date, see:
# https://stackoverflow.com/questions/20789612/convert-julian-date-to-calendar-dates-within-a-data-frame-in-r


# The function ------------------------------------------------------------

# region <- "EAC"
find_current_points <- function(region) {
  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(4)

  # Prepare the data.table
  dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                         lat = nc$dim$lat$vals,
                         t = nc$dim$time$vals)

  ke <- as.data.table(melt(ugos, value.name = "ugos"))
  ke[, c("vgos", "ugosa", "vgosa") := .( # slow
    # as.Date(t, origin = "1950-01-01 00:00:00"), # not necessary to convert as no date calcs needed
    as.vector(ncvar_get(nc, varid = "vgos")) * 100,
    as.vector(ncvar_get(nc, varid = "ugosa")) * 100,
    as.vector(ncvar_get(nc, varid = "vgosa")) * 100
  )
  ]
  nc_close(nc)

  geo <- ke %>%
    na.omit() %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))) %>%
    dplyr::filter(eke > -999) %>%
    dplyr::summarise(mke90 = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke90 = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mke70 = mke90,
                  eke70 = eke90) %>%
    data.table()

  # Read in OISST data
  csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
  events <- fread(paste0(csvDir, "/", region, "-avhrr-only-v2.19810901-20180930_events.csv"))
  int <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(int90 = mean(intensity_mean, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(int70 = int90) %>%
    data.table()

  # Join them
  out <- na.omit(int[geo, on=c("lon==lon", "lat==lat")])

  # return(out)
  # define the regions to mask based on MKE and EKE thresholds
  out$int70[which(out$int70 >= quantile(out$int70, probs = 0.7))] <- -999
  out$mke70[which(out$mke70 >= quantile(out$mke70, probs = 0.7))] <- -999
  out$eke70[which(out$eke70 >= quantile(out$eke70, probs = 0.7))] <- -999
  out$int90[which(out$int90 >= quantile(out$int90, probs = 0.9))] <- -999
  out$mke90[which(out$mke90 >= quantile(out$mke90, probs = 0.9))] <- -999
  out$eke90[which(out$eke90 >= quantile(out$eke90, probs = 0.9))] <- -999

  coordinates(out) <- ~ lat + lon
  gridded(out) <- TRUE
  proj4string(out) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  int70 <- raster(out, layer = "int70")
  mke70 <- raster(out, layer = "mke70")
  eke70 <- raster(out, layer = "eke70")
  int90 <- raster(out, layer = "int90")
  mke90 <- raster(out, layer = "mke90")
  eke90 <- raster(out, layer = "eke90")

  # mask <- rasterToPolygons(mke, fun = function(mke) { mke == 9999 }, dissolve = TRUE)
  int70.mask <- rasterToPolygons(int70, fun = function(int70) { int70 == -999 }, dissolve = TRUE)
  mke70.mask <- rasterToPolygons(mke70, fun = function(mke70) { mke70 == -999 }, dissolve = TRUE)
  eke70.mask <- rasterToPolygons(eke70, fun = function(eke70) { eke70 == -999 }, dissolve = TRUE)
  int90.mask <- rasterToPolygons(int90, fun = function(int90) { int90 == -999 }, dissolve = TRUE)
  mke90.mask <- rasterToPolygons(mke90, fun = function(mke90) { mke90 == -999 }, dissolve = TRUE)
  eke90.mask <- rasterToPolygons(eke90, fun = function(eke90) { eke90 == -999 }, dissolve = TRUE)

  # apply 'over' function, convert coordinates to data frame
  int90.prid <- over(out, int90.mask)
  int90.ptid <- na.omit(int90.prid)
  int90.pt.poly <- out[as.numeric(as.character(row.names(int90.ptid))), ]
  int90.out <- as.data.frame(int90.pt.poly@coords)

  mke90.prid <- over(out, mke90.mask)
  mke90.ptid <- na.omit(mke90.prid)
  mke90.pt.poly <- out[as.numeric(as.character(row.names(mke90.ptid))), ]
  mke90.out <- as.data.frame(mke90.pt.poly@coords)

  eke90.prid <- over(out, eke90.mask)
  eke90.ptid <- na.omit(eke90.prid)
  eke90.pt.poly <- out[as.numeric(as.character(row.names(eke90.ptid))), ]
  eke90.out <- as.data.frame(eke90.pt.poly@coords)

  mask.pts <- list(mke.pt = mke90.out,
                   eke.pt = eke90.out,
                   ev.pt = int90.out)

  save(mask.pts, file = paste0("masks/", region, "-mask_points.Rdata"))

  # check plot
  # pdf(file = paste0(region, "-MKE-mask.pdf"),
  #     width = 5, height = 5, pointsize = 10)
  # plot(to_mask, col = "grey80", cex = 0.5)
  # plot(mask, add = TRUE)
  # plot(pt.poly, add = TRUE, col = "turquoise", cex = 0.5)
  # dev.off()
  #
  return(mask.pts)
}

system.time(AC.pts <- find_current_points("AC"))
system.time(BC.pts <- find_current_points("BC"))
system.time(EAC.pts <- find_current_points("EAC"))
system.time(GS.pts <- find_current_points("GS"))
system.time(KC.pts <- find_current_points("KC"))
