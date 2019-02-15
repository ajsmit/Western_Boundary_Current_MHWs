# 1.1.makeBathy.R

library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)


# Setup -------------------------------------------------------------------

source("setup/regionDefinition.R")
ncDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid"
csvDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
bathy.nc <- paste0(ncDir, "/GEBCO_2014_2D.nc")


# Read bathy function -----------------------------------------------------

# function to extract the dims and data from OISST netCDFs
bathy_nc <- function(ncFile, region = region, csvDir = csvDir, negDeg = FALSE) {
  coords <- bbox[, region]
  reg.adj <- c(0, 0, 360, 360)
  if (negDeg)
    coords <- coords - reg.adj
  nc <- nc_open(ncFile)
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  z <- ncvar_get(nc,
                 varid = "elevation",
                 start = c(LonIdx[1], LatIdx[1]),
                 count = c(length(LonIdx), length(LatIdx)))
  dimnames(z) <- list(lon = round(nc$dim$lon$vals[LonIdx], 4),
                      lat = round(nc$dim$lat$vals[LatIdx], 4))
  nc_close(nc)
  z <-
    as.tibble(melt(z, value.name = "z"), row.names = NULL) %>%
    na.omit()
  if (negDeg)
    z$lon <- z$lon + 360
  fwrite(z,
         file = paste0(csvDir, "/", region, "_bathy.csv"),
         append = FALSE)
  rm(z)
}


# Read bathy --------------------------------------------------------------

bathy_nc(ncFile = bathy.nc, region = "AC", csvDir = csvDir)
bathy_nc(ncFile = bathy.nc, region = "BC", csvDir = csvDir, negDeg = TRUE)
bathy_nc(ncFile = bathy.nc, region = "EAC", csvDir = csvDir)
bathy_nc(ncFile = bathy.nc, region = "GS", csvDir = csvDir, negDeg = TRUE)
bathy_nc(ncFile = bathy.nc, region = "KC", csvDir = csvDir)
