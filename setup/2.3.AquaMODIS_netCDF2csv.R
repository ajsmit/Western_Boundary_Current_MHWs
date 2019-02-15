# 2.3.AquaMODIS_netCDF2csv.R

library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(plyr)
library(lubridate)
library(stringr)
library(doMC); doMC::registerDoMC(cores = 4)


# Read MODIS --------------------------------------------------------------

# setup MODIS
source("setup/regionDefinition.R")
# MODIS.dir <- "/Volumes/Benguela/OceanData/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm"
MODIS.dir <- "/Users/ajsmit/spatial/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm"
# MODIS.csv.dir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/8-day"
MODIS.csv.dir <- "/Users/ajsmit/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/8-day"

# the list of files
ncList <- list.files(path = MODIS.dir, pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
# ncList <- ncList[1:10]
strtDate <- str_sub(basename(ncList[1]), start = 2, end = 8)
endDate <- str_sub(basename(ncList[length(ncList)]), start = 2, end = 8)

# ncFile <- ncList[1]
# region <- "AC"

# function to extract the dims and data from MODIS netCDFs
read_nc <- function(ncFile, region = region, csvDir = csvDir, negDeg = FALSE) {
  coords <- bbox[, region]
  reg.adj <- c(0, 0, 360, 360)
  if (negDeg)
    coords <- coords - reg.adj
  nc <- nc_open(ncFile)
  instrument <- ncatt_get(nc, 0, "instrument")$value
  platform <- ncatt_get(nc, 0, "platform")$value
  product_name <- ncatt_get(nc, 0, "product_name")$value
  fNameStem <- substr(product_name, 17, 38)
  timeStamp <- substr(product_name, 2, 8)
  origin <- paste0(substr(timeStamp, 1, 4), "-01-01")
  date <- as.Date(as.numeric(substr(timeStamp, 5, 7)), origin)
  # create grid (uneven lat/lon step are extracted for some reason...)
  # LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  # LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  northernmost_latitude <- ncatt_get(nc, 0, "northernmost_latitude")$value
  southernmost_latitude <- ncatt_get(nc, 0, "southernmost_latitude")$value
  westernmost_longitude <- ncatt_get(nc, 0, "westernmost_longitude")$value
  easternmost_longitude <- ncatt_get(nc, 0, "easternmost_longitude")$value
  latitude_step <- ncatt_get(nc, 0, "latitude_step")$value
  longitude_step <- ncatt_get(nc, 0, "longitude_step")$value
  number_of_lines <- ncatt_get(nc, 0, "number_of_lines")$value # check
  number_of_columns <- ncatt_get(nc, 0, "number_of_columns")$value # check
  fullLons <- seq(from = westernmost_longitude,
                  to = easternmost_longitude,
                  by = longitude_step)
  fullLats <- seq(from = southernmost_latitude,
                  to = northernmost_latitude,
                  by = latitude_step)
  LonIdx <- which(fullLons > coords[3] & nc$dim$lon$vals < coords[4])
  LatIdx <- which(fullLats > coords[1] & nc$dim$lat$vals < coords[2])
  chl <- ncvar_get(nc,
                   varid = "chlor_a",
                   start = c(LonIdx[1], LatIdx[1]),
                   count = c(length(LonIdx), length(LatIdx))) %>%
    round(4)
  dimnames(chl) <- list(lon = fullLons[LonIdx],
                        lat = fullLats[LatIdx])
  nc_close(nc)
  chl <-
    as.data.table(melt(chl, value.name = "chl"), row.names = NULL) %>%
    mutate(t = date) %>%
    na.omit()
  fwrite(chl,
         file = paste(csvDir, "/", region, "-", instrument, ".",platform, ".",
                      fNameStem, "-", strtDate, "-", endDate, ".csv", sep = ""),
         append = TRUE, col.names = FALSE)
  rm(chl)
}

# note: chl-a to be displayed on a log scale

# apply the function
system.time(llply(ncList, read_nc, region = "AC", csvDir = MODIS.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "BC", csvDir = MODIS.csv.dir,
                  negDeg = TRUE, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "EAC", csvDir = MODIS.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "KC", csvDir = MODIS.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "GS", csvDir = MODIS.csv.dir,
                  negDeg = TRUE, .parallel = TRUE))
