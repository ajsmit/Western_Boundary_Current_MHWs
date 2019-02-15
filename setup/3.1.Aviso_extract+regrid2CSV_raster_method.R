# 3.0.MODIS_extract+regrid2CSV_raster_method.R

# This function regrids the 9km MODIS Aqua data to the coarser OISST 0.25°
# grid resolution. Unlike the other processing approaches, here I use a
# method within the raster package (i.e. I do not read and process the netCDF
# files directly, although I am sure the netCDF package(s) gets used
# internally by raster).

# NOTE: As with all function that need to read a large number of netCDF files
# from disk, the function executes a lot faster if the files reside on the local
# disk; the I/O bottleneck between external hard drives and the CPU is severe.

library(ncdf4)
library(raster)
library(data.table)
library(tidyverse)


# Read MODIS --------------------------------------------------------------

# setup MODIS
source("setup/regionDefinition.R")
MODIS.dir <- "/Volumes/Benguela/OceanData/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm"
MODIS.csv.dir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/8-day"


# the list of files -------------------------------------------------------

ncList <- list.files(path = MODIS.dir, pattern = "*.nc",
                     full.names = TRUE, include.dirs = TRUE)


# output directory --------------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/regridded"


# HELPER FUNCTIONS --------------------------------------------------------

# prepare time vector
makeTimeVec <- function(nc.list) {
  timeStamp <- substr(basename(nc.list), 2, 8)
  origin <- paste0(substr(timeStamp, 1, 4), "-01-01")
  doy <- as.numeric(substr(timeStamp, 5, 7))
  date <- as.Date(as.numeric(substr(timeStamp, 5, 7)), origin)
  return(date)
}

# create an extent object for use in cropping
makeLims <- function(region) {
  ex <- bbox[, region]
  out <- c(ex[3], ex[4], ex[1], ex[2])
  return(out)
}

# make the output file name
makeFileName <- function(nc.list) {
  nc <- nc_open(nc.list[1])
  instrument <- ncatt_get(nc, 0, "instrument")$value
  platform <- ncatt_get(nc, 0, "platform")$value
  product_name <- ncatt_get(nc, 0, "product_name")$value
  fNameStem <- substr(product_name, 17, 34)
  timeStampStart <- substr(product_name, 2, 8)
  originStart <- paste0(substr(timeStampStart, 1, 4), "-01-01")
  strtDate <- as.Date(as.numeric(substr(timeStampStart, 5, 7)), originStart)
  nc_close(nc)
  nc <- nc_open(nc.list[length(nc.list)])
  product_name <- ncatt_get(nc, 0, "product_name")$value
  timeStampEnd <- substr(product_name, 2, 8)
  originEnd <- paste0(substr(timeStampEnd, 1, 4), "-01-01")
  endDate <- as.Date(as.numeric(substr(timeStampEnd, 5, 7)), originEnd)
  nc_close(nc)
  fileName <- paste0(".", instrument, ".", platform, ".", fNameStem, "_regrid_",
                     strtDate, "-", endDate, ".csv")
  return(fileName)
}


# MAIN FUNCTION -----------------------------------------------------------

# region <- "AC"
# nc.list <- ncList
regridChl <- function(nc.list, region, csvDir, negDeg = FALSE) {
  # create a raster stack
  # FYI: read about the difference between a raster stack and a raster brick...
  chl.in <- do.call(stack,
                    lapply(nc.list,
                           function(nc)
                             raster(nc, varname = 'chlor_a', level = 1)))

  # crop
  lims <- makeLims(region)
  if (negDeg) {
    reg.adj <- c(360, 360, 0, 0)
    lims <- lims - reg.adj
    }
  chl.hr <- crop(chl.in, extent(lims))

  # find and setup pre-made low res (OISST) grid
  grid.path <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/grids"
  new.grid <- fread(paste0(grid.path, "/", region, "_year_grid.csv")) %>%
    dplyr::select(lon, lat) %>%
    unique()
  coordinates(new.grid) <- ~ lon + lat
  gridded(new.grid) <- TRUE
  new.grid <- raster(new.grid)
  if (negDeg)
    new.grid <- raster::rotate(new.grid)

  # resample
  chl.lr <- resample(chl.hr, new.grid)

  # set time values
  chl.lr <- setZ(chl.lr, makeTimeVec(nc.list))

  # reset lons to 0:360
  if (negDeg)
    chl.lr <- raster::shift(raster::rotate(raster::shift(chl.lr, 180)), 180)

  # create long data.frame
  out <- as.data.frame(chl.lr, xy = TRUE, long = TRUE)
  names(out) <- c("lon", "lat", "time", "chlor_a")

  # output file
  fileName <- makeFileName(nc.list)
  fwrite(out,
         file = paste0(csvDir, "/", region, fileName),
         append = FALSE, col.names = TRUE)
  # return(out)
}


# execute the function ----------------------------------------------------

regridChl(ncList, "AC", csvDir)
regridChl(ncList, "BC", csvDir, negDeg = TRUE)
regridChl(ncList, "EAC", csvDir)
regridChl(ncList, "GS", csvDir, negDeg = TRUE)
regridChl(ncList, "KC", csvDir)
