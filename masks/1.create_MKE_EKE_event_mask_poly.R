# Calculate the mean (long-term) MKE, EKE and mean intensity for masking

library(tidyverse)
library(data.table) # still required for fread() and fwrite()
library(raster)
library(sf)
library(rgeos)
library(smoothr)
library(ncdf4)
library(RcppRoll)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 3)

# region <- "EAC"

mask.fun <- function(region) {
  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ke <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(2)

  dimnames(ke) <- list(lon = nc$dim$lon$vals,
                       lat = nc$dim$lat$vals,
                       t = nc$dim$time$vals)

    # Using data.table...
  # ke <- reshape2::melt(ke, value.name = "ugos")
  # ke <- as.data.table(ke)
  #
  # ke[, c("vgos", "ugosa", "vgosa") := .( # slow
  #   # as.Date(t, origin = "1950-01-01 00:00:00"), # not necessary to convert as no date calcs needed
  #   round(as.vector(ncvar_get(nc, varid = "vgos")) * 100, 2),
  #   round(as.vector(ncvar_get(nc, varid = "ugosa")) * 100, 2),
  #   round(as.vector(ncvar_get(nc, varid = "vgosa")) * 100, 2)
  # )
  # ]

  # ke[, eke := 0.5 * (vgosa^2 + ugosa^2)]
  # ke[, c("ugosa", "vgosa") := NULL]
  # ke[, eke := roll_mean(eke, n = 30, align = "center", na.rm = TRUE, fill = c(NA, NA, NA)),
  #    by = .(lon, lat)]
  #
  # ke <- ke[, .(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
  #              eke = mean(eke, na.rm = TRUE)),
  #          by = .(lon, lat)
  #          ]
  # ke <- na.omit(ke)

  ke <- reshape2::melt(ke, value.name = "ugos") # slower, but more memory efficient than data.table?
  ke <- ke %>%
    dplyr::mutate(vgos = round(as.vector(ncvar_get(nc, varid = "vgos")) * 100, 2),
                  ugosa = round(as.vector(ncvar_get(nc, varid = "ugosa")) * 100, 2),
                  vgosa = round(as.vector(ncvar_get(nc, varid = "vgosa")) * 100, 2)) %>%
    dplyr::mutate(eke = 0.5 * (vgosa^2 + ugosa^2)) %>%
    dplyr::select(-ugosa, -vgosa) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", na.rm = TRUE, fill = c(NA, NA, NA))) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    na.omit()
  nc_close(nc)

  # Read in OISST data
  csvDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
  col_types <- list(
    lon = col_double(),
    lat = col_double(),
    event_no = col_integer(),
    index_start = col_integer(),
    index_peak = col_integer(),
    index_end = col_integer(),
    duration = col_integer(),
    date_start = col_date(),
    date_peak = col_date(),
    date_end = col_date(),
    intensity_mean = col_double(),
    intensity_max = col_double(),
    intensity_var = col_double(),
    intensity_cumulative = col_double(),
    intensity_mean_relThresh = col_double(),
    intensity_max_relThresh = col_double(),
    intensity_var_relThresh = col_double(),
    ntensity_cumulative_relThresh = col_double(),
    intensity_mean_abs = col_double(),
    intensity_max_abs = col_double(),
    intensity_var_abs = col_double(),
    intensity_cumulative_abs = col_double(),
    rate_onset = col_double(),
    rate_decline = col_double()
  )
  events <- vroom(
    file       = paste0(csvDir, "/", region, "-avhrr-only-v2.19810901-20180930_events.csv"),
    skip       = 1,
    # col_select = c(lon, lat, intensity_mean),
    delim      = ",",
    col_names  = names(col_types),
    col_types  = col_types,
    na         = c("", "NA", "NULL"))
  int <- events %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(int = mean(intensity_mean, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    na.omit()

  # Join them
  # out <- na.omit(int[ke, on=c("lon==lon", "lat==lat")])

  out <- int %>%
    dplyr::left_join(ke, by = c("lon", "lat")) %>%
    na.omit()
  rm(int); rm(ke)

  # define the regions to mask based on MKE and EKE thresholds
  out$int[which(out$int >= quantile(out$int, probs = 0.9))] <- -999
  out$mke[which(out$mke >= quantile(out$mke, probs = 0.9))] <- -999
  out$eke[which(out$eke >= quantile(out$eke, probs = 0.9))] <- -999

  coordinates(out) <- ~ lon + lat
  gridded(out) <- TRUE
  proj4string(out) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  # out.st <- stack(out)
  int <- raster(out, layer = "int")
  mke <- raster(out, layer = "mke")
  eke <- raster(out, layer = "eke")
  rm(out)

  mask.list <- list(int = smooth(rasterToPolygons(int, fun = function(int) { int == -999 },
                                                  dissolve = TRUE), method = "ksmooth", smoothness = 3),
                    mke = smooth(rasterToPolygons(mke, fun = function(mke) { mke == -999 },
                                                  dissolve = TRUE), method = "ksmooth", smoothness = 3),
                    eke = smooth(rasterToPolygons(eke, fun = function(eke) { eke == -999 },
                                                  dissolve = TRUE), method = "ksmooth", smoothness = 3))
  save(mask.list, file = paste0("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/", region, "-mask_polys.RData"))
  return(mask.list)
  rm(mask.list); gc()
}

AC.mask <- mask.fun("AC")
BC.mask <- mask.fun("BC")
EAC.mask <- mask.fun("EAC")
GS.mask <- mask.fun("GS")
KC.mask <- mask.fun("KC")


# Make some plots ---------------------------------------------------------

# read in the region-specific components
source("setup/plot.layers.R")

# plot function
# plot.parameters <- EAC.layers
# bathy <- AC_bathy
# data <- mask.list
# str(data$int75)
# region <- "EAC"

poly.plot <- function(region, plot.parameters) {
  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))
  maskDir <- "/Volumes/GoogleDrive/My Drive/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke <- fortify(mask.list$mke)
  eke <- fortify(mask.list$eke)
  int <- fortify(mask.list$int)
  plt <- ggplot() +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = alpha("#c83c7e", 1.0), colour = NA) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    geom_polygon(data = eke, aes(long, lat, group = group),
                 fill = alpha("gray50", 0.5), colour = "navy", size = 0.3) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    guides(alpha = "none",
           fill = "none",
           size = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters

  return(plt)
}

AC.poly.plt <- poly.plot("AC", AC.layers)
BC.poly.plt <- poly.plot("BC", BC.layers)
EAC.poly.plt <- poly.plot("EAC", EAC.layers)
GS.poly.plt <- poly.plot("GS", GS.layers)
KC.poly.plt <- poly.plot("KC", KC.layers)

bbox <- data.frame("AC" = c(-45, -35, 10, 25),
                   "BC" = c(-45, -35, 305, 310),
                   "EAC" = c(-40, -32.5, 150, 157.5),
                   "GS" = c(35, 45, 290, 300),
                   "KC" = c(32.5, 42.5, 140, 150),
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))

fig_masks <- ggarrange(AC.poly.plt + annotate("rect", xmin = bbox["lonmin", "AC"], xmax = bbox["lonmax", "AC"],
                                              ymin = bbox["latmin", "AC"], ymax = bbox["latmax", "AC"],
                                              col = "black", size = 0.4, fill = NA),
                       BC.poly.plt + annotate("rect", xmin = bbox["lonmin", "BC"], xmax = bbox["lonmax", "BC"],
                                              ymin = bbox["latmin", "BC"], ymax = bbox["latmax", "BC"],
                                              col = "black", size = 0.4, fill = NA),
                       EAC.poly.plt + annotate("rect", xmin = bbox["lonmin", "EAC"], xmax = bbox["lonmax", "EAC"],
                                               ymin = bbox["latmin", "EAC"], ymax = bbox["latmax", "EAC"],
                                               col = "black", size = 0.4, fill = NA),
                       GS.poly.plt + annotate("rect", xmin = bbox["lonmin", "GS"], xmax = bbox["lonmax", "GS"],
                                              ymin = bbox["latmin", "GS"], ymax = bbox["latmax", "GS"],
                                              col = "black", size = 0.4, fill = NA),
                       KC.poly.plt + annotate("rect", xmin = bbox["lonmin", "KC"], xmax = bbox["lonmax", "KC"],
                                              ymin = bbox["latmin", "KC"], ymax = bbox["latmax", "KC"],
                                              col = "black", size = 0.4, fill = NA),
                       ncol = 1, nrow = 5, labels = "AUTO")

ggplot2::ggsave("publ_plots/Figure_2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
