# Get and prepare data OISST event data -----------------------------------

# Get and prepare Aviso data ----------------------------------------------

# Regridding performed with NCO
# First, make time into a record dimension (from a fixed dimension);
# do this for each of the the downloaded Aviso+ files:
# > ncks --mk_rec_dmn time infile.nc -O outfile.nc
# (the above are in the shell script in the Aviso+ netCDF directories, e.g.
# /Volumes/Benguela/Aviso/global/delayed-time/grids/msl/all-sat-merged/AC/unlimDim.sh)

# Now concatenate these files along the record dimension (i.e. time):
# > ncrcat infile????.nc outfile.nc
# (see the shell script for additional details)

# Make a 'template' file from the OISST data, which will be used as the grid to
# interpolate to. Since these data are global, they need to be subsetted first;
# this is termed creating a hyperslab in NCO speak. Thus, create a hyperslab that
# corresponds to each WBC region's spatial extent, e.g. for the Agulhas Current:
# > ncks -F -d lon,6.25,45. -d lat,-45.,-20. avhrr-only-v2.20180930.nc AC_template.nc
#
# The template file is just a daily netCDF files of the region of interest.

# Now interpolate the Aviso file to the OISST template file:
# > ncremap -i outfile.nc -d AC_template.nc -o outfile_regridded.nc

library(tidyverse)
library(data.table)
library(ncdf4)
library(lubridate)
library(colorspace)
library(scales)
library(RcppRoll)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# read in the region-specific components
source("setup/plot.layers.R")


# Bathy data --------------------------------------------------------------

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC_bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC_bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC_bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS_bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC_bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

# Some functions ----------------------------------------------------------

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# region <- "EAC"

# Exceedence above threshold per month:
# this function calculates for each month and each pixel the
# mean MKE and EKE (from Aviso+) and mean mean MHW event intensity
uber.fun <- function(region) {
  clim.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
  clim.base <- "-avhrr-only-v2.19810901-20180930_climatology.csv"
  clim.dat <- paste0(region, clim.base)
  clim <- fread(paste0(clim.dir, "/", clim.dat))

  clim[, t := fastDate(t)]
  ex <- clim[t >= "1993-01-01" & t <= "2018-06-10", ] # trim to match start date of Aviso data
  ex[, ex := roll_mean((temp - thresh), n = 30, align = "center"), by = .(lon, lat)]
  ex[, c("doy", "temp", "seas", "thresh") := NULL]
  rm(clim)

  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(4)

  dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                         lat = nc$dim$lat$vals,
                         t = nc$dim$time$vals)
  ke <- as.data.table(melt(ugos, value.name = "ugos"))
  ke[, c("t", "vgos", "ugosa", "vgosa") := .( # slow
    as.Date(t, origin = "1950-01-01 00:00:00"),
    as.vector(ncvar_get(nc, varid = "vgos")) * 100,
    as.vector(ncvar_get(nc, varid = "ugosa")) * 100,
    as.vector(ncvar_get(nc, varid = "vgosa")) * 100
  )
  ]
  nc_close(nc)

  ke[, eke := (0.5 * ((vgosa)^2 + (ugosa)^2))]
  ke[, c("mke_roll", "eke_roll1", "eke_roll2") := .(
    0.5 * (roll_mean(vgos, n = 30, align = "center")^2 + roll_mean(ugos, n = 30, align = "center")^2),
    0.5 * ((vgos - roll_mean(vgos, n = 30, align = "center"))^2 + (ugos - roll_mean(ugos, n = 30, align = "center"))^2),
    roll_mean(eke, n = 30, align = "center")
  ), by = .(lon, lat)
  ]
  ke <- na.omit(ke)
  ke[, c("ugos", "vgos", "ugosa", "vgosa") := NULL]

  # Join the data and do correlation
  # Join the tables
  setkey(ex, lon, lat, t)
  setkey(ke, lon, lat, t)
  dat <- ex[ke, nomatch = 0]
  dat <- dat[eke_roll1 > 0, ]
  dat <- na.omit(dat)

  mke.cor.fun <- function(df) cor(df$mke_roll, df$ex)
  eke1.cor.fun <- function(df) cor(df$eke_roll1, df$ex)
  eke2.cor.fun <- function(df) cor(df$eke_roll2, df$ex)
  dat.r <- dat %>%
    dplyr::group_by(lon, lat) %>%
    nest() %>%
    dplyr::mutate(r.mke = map(data, mke.cor.fun)) %>%
    dplyr::mutate(r.eke1 = map(data, eke1.cor.fun)) %>%
    dplyr::mutate(r.eke2 = map(data, eke2.cor.fun)) %>%
    unnest(r.mke, r.eke1, r.eke2)
  return(dat.r)
}

# Apply the function
system.time(AC.cor <- uber.fun("AC"))
save(AC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/AC_correlations_v3.RData")
rm(AC.cor)
system.time(BC.cor <- uber.fun("BC"))
save(BC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/BC_correlations_v3.RData")
rm(BC.cor)
system.time(EAC.cor <- uber.fun("EAC"))
save(EAC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/EAC_correlations_v3.RData")
rm(EAC.cor)
system.time(GS.cor <- uber.fun("GS"))
save(GS.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/GS_correlations_v3.RData")
rm(GS.cor)
system.time(KC.cor <- uber.fun("KC"))
save(KC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/KC_correlations_v3.RData")
rm(KC.cor)

# summary(EAC.cor)
#
# save(AC.cor, BC.cor, EAC.cor, GS.cor, KC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/correlations_v3.RData")

load(paste0(corDir, "/", "AC_correlations_v3.RData"))
load(paste0(corDir, "/", "BC_correlations_v3.RData"))
load(paste0(corDir, "/", "EAC_correlations_v3.RData"))
load(paste0(corDir, "/", "KC_correlations_v3.RData"))
load(paste0(corDir, "/", "GS_correlations_v3.RData"))

region <- "AC"

# Make a graph
gv.plot <- function(data, plot.parameters, bathy, region) {
  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-masks.RData"))
  mke <- fortify(mask.list$mke90)
  eke <- fortify(mask.list$eke90)
  int <- fortify(mask.list$int90)

  fig.mke <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = r.mke)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE) +
    # scale_fill_gradientn(colours = col1,
    #                      values = scales::rescale(c(min(data$r.mke), 0, max(data$r.mke)))) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = NA, colour = "purple", size = 0.3) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 fill = NA, colour = "red3", size = 0.3) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "r",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(40, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters

  fig.eke2 <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = r.eke2)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE) +
    # scale_fill_gradientn(colours = col1,
    #                      values = scales::rescale(c(min(data$r.eke2), 0, max(data$r.eke2)))) +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = int, aes(long, lat, group = group),
                 fill = NA, colour = "purple", size = 0.3) +
    geom_polygon(data = eke, aes(long, lat, group = group),
                 fill = NA, colour = "navy", size = 0.3) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "r",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(40, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters


  fig <- ggarrange(fig.mke,
                   fig.eke2,
                   ncol = 2, nrow = 1)
  return(fig)
}

AC.fig0x <- gv.plot(AC.cor, AC.layers, AC_bathy, "AC")
BC.fig0x <- gv.plot(BC.cor, BC.layers, BC_bathy, "BC")
EAC.fig0x <- gv.plot(EAC.cor, EAC.layers, EAC_bathy, "EAC")
GS.fig0x <- gv.plot(GS.cor, GS.layers, GS_bathy, "GS")
KC.fig0x <- gv.plot(KC.cor, KC.layers, KC_bathy, "KC")

fig0x <- ggarrange(AC.fig0x,
                   BC.fig0x,
                   EAC.fig0x,
                   GS.fig0x,
                   KC.fig0x,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("publ_plots/Fig0x_correlations_v4.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)
