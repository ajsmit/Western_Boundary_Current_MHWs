# 2.4.OC_cci_netCDF2csv.R

library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(plyr)
library(lubridate)
library(stringr)
library(readr)
library(doMC); doMC::registerDoMC(cores = 4)


# Read OC_cci --------------------------------------------------------------

# setup OC_cci
source("setup/regionDefinition.R")
OC_cci.dir <- "/Volumes/Benguela/OceanData/Ocean_Colour_cci_Version_3.1 "
OC_cci.csv.dir <- "/Volumes/Benguela/spatial/processed/Ocean_Colour_cci_Version_3.1/WBCs/8-day"

# function to extract the dims and data from OC_cci netCDFs
# hard-code file name info into function
read_nc <- function(ncDir, region = region, csvDir = csvDir, negDeg = FALSE) {
  coords <- bbox[, region]
  if (negDeg) {
    reg.adj <- c(0, 0, 360, 360)
    coords <- coords - reg.adj
    }
  nc <- nc_open(paste0(ncDir, "/", region, "_CCI_ALL-v3.1-8DAY.nc"))
  fNameStem <- "_CCI_ALL-v3.1-8DAY"
  chl <- ncvar_get(nc, varid = "chlor_a") %>%
    round(4)
  dimnames(chl) <- list(lon = nc$dim$lon$vals,
                        lat = nc$dim$lat$vals,
                        t = nc$dim$time$vals)
  nc_close(nc)
  chl <-
    as.data.table(melt(chl, value.name = "chl"), row.names = NULL) %>%
    mutate(t = as.Date(t, origin = "1970-01-01 00:00:00")) %>%
    na.omit()
  strtDate <- head(chl$t, 1)
  endDate <- tail(chl$t, 1)
  fwrite(chl,
         file = paste(csvDir, "/", region, fNameStem, "-", strtDate, "-", endDate, ".csv", sep = ""),
         append = FALSE, col.names = TRUE, nThread = 2) # 3 threads do not seem to work
  rm(chl); rm(strtDate); rm(endDate)
}

# note: chl-a to be displayed on a log scale

# apply the function
system.time(llply(OC_cci.dir, read_nc, "AC", OC_cci.csv.dir, .parallel = TRUE))
system.time(llply(OC_cci.dir, read_nc, "BC", OC_cci.csv.dir,
                  negDeg = TRUE, .parallel = TRUE))
system.time(llply(OC_cci.dir, read_nc, "EAC", OC_cci.csv.dir, .parallel = TRUE))
system.time(llply(OC_cci.dir, read_nc, "KC", OC_cci.csv.dir, .parallel = TRUE))
system.time(llply(OC_cci.dir, read_nc, "GS", OC_cci.csv.dir,
                  negDeg = TRUE, .parallel = TRUE))
