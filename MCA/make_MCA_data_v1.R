library(tidyverse)
library(ncdf4)
library(data.table)
library(RcppRoll)
doParallel::registerDoParallel(cores = 3)
source("/Users/ajsmit/Dropbox/R/misc_R/anom.detr.R")
source("setup/functions.R")

region <- "EAC"


# INITIAL SETUP -----------------------------------------------------------

# Lanczos filter weights
# I ended up not using these...
source("/Users/ajsmit/Dropbox/R/misc_R/lanczos_weights.R")
lw1 <- lanczos_weights(nwt = 101, srt = 1, ihp = "lowpass", fca = 1 / (1 * 365))
lw2 <- lanczos_weights(nwt = 201, srt = 1, ihp = "lowpass", fca = 1 / (2 * 365))
lw3 <- lanczos_weights(nwt = 301, srt = 1, ihp = "lowpass", fca = 1 / (3 * 365))

# start and end dates for MCA output data
start_date <- "2006-01-01"
end_date <- "2015-12-31"

# climatology for EKE data
climatologyPeriod <- c(as.Date("1994-01-01"), as.Date("2017-12-31"))

# the output
out_dir <- "/Volumes/Benguela/spatial/processed/WBC/misc_results"
out_base <- "_MCA_df_"


# FUNCTION TO PREPARE EXCEEDANCE DATA -------------------------------------

load_exceedance <- function(region) {
  ex_dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
  ex_base <- "-avhrr-only-v2.19810901-20180930_climatology.csv"
  ex <- fread(paste0(ex_dir, "/", region, ex_base))
  ex[, t := fastDate(t)]
  ex[, ex := roll_mean((temp - thresh), n = 5, align = "center", fill = c(NA, NA, NA)), by = .(lon, lat)]
  # ex[, ex := ifelse(ex <= 0, 0, 1)]
  ex[, c("doy", "temp", "seas", "thresh") := NULL]
  ex_df <- ex %>%
    as_tibble() %>%
    dplyr::filter(t >= start_date & t <= end_date)
  rm(ex)
  colnames(ex_df) <- c("lon", "lat", "t", "z")
  return(ex_df)
}


# FUNCTION TO PREPARE THE EKE DATA ----------------------------------------

load_geostrophic <- function(region) {
  aviso_dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso_dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso_dir, "/", region, "/", aviso_dat))
  ugosa <- ncvar_get(nc, varid = "ugosa") * 100 %>%
    round(4)
  dimnames(ugosa) <- list(lon = nc$dim$lon$vals,
                          lat = nc$dim$lat$vals,
                          t = nc$dim$time$vals)
  ke <- as.data.table(melt(ugosa, value.name = "ugosa"))
  ke[, c("t", "vgosa") := .(
    as.Date(t, origin = "1950-01-01 00:00:00"),
    as.vector(ncvar_get(nc, varid = "vgosa")) * 100
  )
  ]
  nc_close(nc)

  out <- ke %>%
    na.omit() %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::mutate(eke = 0.5 * (vgosa^2 + ugosa^2)) %>%
    dplyr::select(-ugosa, -vgosa, -n)
  rm(ke)
  return(out)
}

# lw <- NULL
# tag <- "unfiltered"

# ASSEMBLE FUNCTION -------------------------------------------------------

# anomalise (this step takes long, hence data are saved after calc)...
# calculate the anomaly (zero-centered), remove linear trend,
# calculate exceedance, and apply 30-day running mean smooth
assemble <- function(lw, tag) {
  eke_df <- ddply(out, .(lon, lat),
                  anom.detr, x = t, y = eke, climatologyPeriod = climatologyPeriod, lw,
                  .parallel = TRUE)

  # merge the exceedance and EKE data
  eke_df <- eke_df %>%
    as_tibble() %>%
    dplyr::filter(x >= start_date & x <= end_date)
  colnames(eke_df) <- c("lon", "lat", "t", "z")

  # merge the data
  data_full <- eke_df %>%
    dplyr::inner_join(ex_df, c("t", "lon", "lat"))
  colnames(data_full)[c(4, 5)] <- c("eke", "ex")
  rm(eke_df)

  # these data are used in python in the MCA
  fwrite(data_full, file = paste0(out_dir, "/", region, out_base, "-",
                                  start_date, "-", end_date, "-", tag, ".csv"))
  rm(data_full)
}


# GET READY: SET REGION, LOAD DATA, AND GO! -------------------------------


# Agulhas Current ---------------------------------------------------------

region <- "AC"
ex_df <- load_exceedance(region)
out <- load_geostrophic(region)
assemble(NULL, "unfiltered")
# assemble(lw2, "unfiltered")
# assemble(lw3, "unfiltered")
rm(ex_df); rm(out)


# Brazil Current ----------------------------------------------------------

region <- "BC"
ex_df <- load_exceedance(region)
out <- load_geostrophic(region)
assemble(NULL, "unfiltered")
# assemble(lw2, "unfiltered")
# assemble(lw3, "unfiltered")
rm(ex_df); rm(out)


# East Australian Current -------------------------------------------------

region <- "EAC"
ex_df <- load_exceedance(region)
out <- load_geostrophic(region)
assemble(NULL, "unfiltered_bin")
# assemble(NULL, "unfiltered")
# assemble(NULL, "unfiltered")
rm(ex_df); rm(out)


# Gulf Stream -------------------------------------------------------------

region <- "GS"
ex_df <- load_exceedance(region)
out <- load_geostrophic(region)
assemble(NULL, "unfiltered")
# assemble(lw2, "unfiltered")
# assemble(lw3, "unfiltered")
rm(ex_df); rm(out)


# Kuroshio Current --------------------------------------------------------

region <- "KC"
ex_df <- load_exceedance(region)
out <- load_geostrophic(region)
assemble(NULL, "unfiltered")
# assemble(lw2, "unfiltered")
# assemble(lw3, "unfiltered")
rm(ex_df); rm(out)

