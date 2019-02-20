library(ncdf4)
library(tidyverse)
library(lubridate)
library(data.table)
library(colorspace)
library(raster)


# Setup -------------------------------------------------------------------

source("setup/regionDefinition.R")
# for converting julian date to calendar date, see:
# https://stackoverflow.com/questions/20789612/convert-julian-date-to-calendar-dates-within-a-data-frame-in-r


# Extract eddy tracks -----------------------------------------------------

# one file to support all WBCs
ncFile <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/value-added/eddy-trajectory/eddy_trajectory_2.0exp_19930101_20180118.nc"
nc <- nc_open(ncFile)
eddies <- tibble(lat = round(ncvar_get(nc, varid = "latitude"), 4), # observation longitude
                 lon = round(ncvar_get(nc, varid = "longitude"), 4), # observation latitude
                 # observation sequence number, days from eddy start:
                 observation_number = ncvar_get(nc, varid = "observation_number"))
nc_close(nc); rm(nc)


# The function ------------------------------------------------------------

# region <- "BC"
find_current_points <- function(eddies, region) {
  # create unique lons/lats for the eddy database
  to_mask <- eddies %>%
    dplyr::filter(lat >= bbox[1, region] & lat <= bbox[2, region]) %>%
    dplyr::filter(lon >= bbox[3, region] & lon <= bbox[4, region]) %>%
    dplyr::filter(observation_number == 0) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    unique()
  coordinates(to_mask) <- ~ lon + lat
  proj4string(to_mask) <- CRS("+proj=utm +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # create a polygon of the main AC current
  aviso.dir <- "~/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") %>%
    round(4)
  dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                         lat = nc$dim$lat$vals,
                         time = nc$dim$time$vals)
  mke <- # slow!
    as_tibble(melt(ugos, value.name = "ugos"), row.names = NULL) %>%
    dplyr::mutate(vgos = as.vector(ncvar_get(nc, varid = "vgos")),
                  time = as.Date(time, origin = "1950-01-01 00:00:00"),
                  year = lubridate::year(time),
                  month = lubridate::month(time)) %>%
    dplyr::mutate(ugos = ugos * 100,
                  vgos = vgos * 100) %>%
    dplyr::group_by(lat, lon) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2)) %>%
    dplyr::ungroup() %>%
    na.omit()
  nc_close(nc)

  # define the region to mask
  mke$mke[which(mke$mke >= quantile(mke$mke, probs = 0.9))] <- 9999

  coordinates(mke) <- ~ lon + lat
  gridded(mke) <- TRUE
  proj4string(mke) <- CRS("+proj=utm +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  mke <- raster(mke)

  # mask <- rasterToPolygons(mke, fun = function(mke) { mke == quantile(mke, probs = 0.9) })
  mask <- rasterToPolygons(mke, fun = function(mke) { mke == 9999 }, dissolve = TRUE)


  # apply over function, convert coordinates to data frame
  prid <- over(to_mask, mask)
  ptid <- na.omit(prid)
  pt.poly <- to_mask[as.numeric(as.character(row.names(ptid))), ]
  out <- pt.poly@coords

  fwrite(out, file = paste0(region, "-eddy_pts.csv"), col.names = TRUE)

  # check plot
  pdf(file = paste0(region, "-mask.pdf"),
      width = 5, height = 5, pointsize = 10)
  plot(to_mask, col = "grey80", cex = 0.5)
  plot(mask, add = TRUE)
  plot(pt.poly, add = TRUE, col = "turquoise", cex = 0.5)
  dev.off()

  return(out)
}

system.time(AC_pts <- find_current_points(eddies, "AC"))
system.time(BC_pts <- find_current_points(eddies, "BC"))
system.tim(EAC_pts <- find_current_points(eddies, "EAC"))
system.time(GS_pts <- find_current_points(eddies, "GS"))
system.time(KC_pts <- find_current_points(eddies, "KC"))
