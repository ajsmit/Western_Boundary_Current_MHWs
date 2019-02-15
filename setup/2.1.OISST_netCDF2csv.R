# 2.1.OISST_netCDF2csv.R

library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(plyr)
library(lubridate)
library(stringr)
library(doMC); doMC::registerDoMC(cores = 4)


# Read OISST --------------------------------------------------------------

# setup OISST
source("setup/regionDefinition.R")
OISST.dir <- "/Volumes/Benguela/OceanData/OISSTv2/daily/netCDF"
OISST.csv.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"

# function to extract the dims and data from OISST netCDFs
read_nc <- function(ncFile, region = region, csvDir = csvDir) {
  coords <- bbox[, region]
  nc <- nc_open(ncFile)
  pathLen <- nchar(OISST.dir) + 1 # to account for the "/" that needs to be inserted
  fNameStem <-
    substr(ncFile, pathLen + 1, pathLen + 13)
  fDate <- substr(ncFile, pathLen + 15, pathLen + 22)
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  sst <- ncvar_get(nc,
                   varid = "sst",
                   start = c(LonIdx[1], LatIdx[1], 1, 1),
                   count = c(length(LonIdx), length(LatIdx), 1, 1)) %>%
    round(4)
  dimnames(sst) <- list(lon = nc$dim$lon$vals[LonIdx],
                        lat = nc$dim$lat$vals[LatIdx])
  nc_close(nc)
  sst <-
    as.data.table(melt(sst, value.name = "temp"), row.names = NULL) %>%
    mutate(t = ymd(fDate)) %>%
    na.omit()
  fwrite(sst,
         file = paste(csvDir, "/", region, "-", fNameStem, ".", strtDate, "-", endDate, ".csv", sep = ""),
         append = TRUE, col.names = FALSE)
  rm(sst)
}

# the list of files
ncList <- list.files(path = OISST.dir, pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
strtDate <- str_sub(ncList[1], start = 64, end = 71)
endDate <- str_sub(ncList[length(ncList)], start = 64, end = 71)

# apply the function
system.time(llply(ncList, read_nc, region = "AC", csvDir = OISST.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "BC", csvDir = OISST.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "EAC", csvDir = OISST.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "KC", csvDir = OISST.csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, region = "GS", csvDir = OISST.csv.dir, .parallel = TRUE))
