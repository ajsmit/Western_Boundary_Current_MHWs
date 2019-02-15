# 2.2.Aviso+_netCDF2csv.R

library(ncdf4)
library(plyr)
library(data.table)
library(tidyverse)
library(reshape2)
# library(lubridate)
# library(stringr)
library(doMC); doMC::registerDoMC(cores = 1)


# Read Aviso+ function ----------------------------------------------------

ncDir <- "/Volumes/Benguela/Aviso+/global/delayed-time/grids/msl/all-sat-merged"
csvDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"

# ncFile <- ncList[8]
# region <- "GS"
# function to extract the dims and data from CMC netCDFs
read_nc <- function(ncDir = ncDir, region = region, csvDir = csvDir) {
  ncList <- list.files(path = paste0(ncDir, "/", region), pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
  ncFun <- function(ncFile = ncFile, region = region, csvDir = csvDir) {
    nc <- nc_open(ncFile)
    fNameStem <- "dataset-duacs-rep-global-merged-allsat-phy-l4-v3"
    ugos <- ncvar_get(nc, varid = "ugos") %>%
      round(4)
    dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                           lat = nc$dim$lat$vals,
                           time = nc$dim$time$vals)
    sla <-
      as.tibble(melt(ugos, value.name = "ugos"), row.names = NULL) %>%
      dplyr::mutate(vgos = as.vector(ncvar_get(nc, varid = "vgos")),
                    ugosa = as.vector(ncvar_get(nc, varid = "ugosa")),
                    vgosa = as.vector(ncvar_get(nc, varid = "vgosa")),
                    sla = as.vector(ncvar_get(nc, varid = "sla")),
                    adt = as.vector(ncvar_get(nc, varid = "adt"))) %>%
      dplyr::mutate(time = as.Date(time, origin = "1950-01-01 00:00:00")) %>%
      na.omit()
    nc_close(nc)
    fwrite(sla,
           file = paste0(csvDir, "/", region, "-", fNameStem, ".csv"),
           append = TRUE, col.names = FALSE)
    rm(sla)
  }
  out <- llply(ncList, ncFun, region = region, csvDir = csvDir, .parallel = FALSE) # single core when writing
  # to preserve date order
  return(out)
}

# headers not saved, but they are lon, lat, time, ugos, vgos, ugosa, vgosa, sla, adt
read_nc(ncDir, "AC", csvDir)
read_nc(ncDir, "BC", csvDir)
# read_nc(ncDir, "EAC", csvDir)
read_nc(ncDir, "KC", csvDir)
read_nc(ncDir, "GS", csvDir)


# Regrid them -------------------------------------------------------------

regDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily_regridded"
