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

# Now interpolate the Aviso file to the OISST template file:
# > ncremap -i outfile.nc -d AC_template.nc -o outfile_regridded.nc

library(data.table)
library(ncdf4)
library(lubridate)
library(colorspace)
library(scales)
library(ggpubr)
library(tidyverse)
library(doMC); doMC::registerDoMC(cores = 4)

# read in the region-specific components
source("setup/plot.layers.R")

# Some functions ----------------------------------------------------------

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# region <- "AC"

# Exceedence above threshold per month:
# this function calculates for each month and each pixel the
# mean MKE and EKE (from Aviso+) and mean mean MHW event intensity
uber.fun <- function(region) {
  clim.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
  clim.base <- "-avhrr-only-v2.19810901-20180930_climatology.csv"
  clim.dat <- paste0(region, clim.base)
  clim <- as_tibble(fread(paste0(clim.dir, "/", clim.dat)))

  ex <- clim %>%
    dplyr::mutate(t = fastDate(t)) %>%
    dplyr::filter(t >= "1993-01-01") %>%
    dplyr::mutate(time = floor_date(t, unit = "month")) %>%
    dplyr::group_by(lon, lat, time) %>%
    dplyr::summarise(ex = mean(temp - thresh, na.rm = TRUE)) %>%
    dplyr::ungroup()
  rm(clim)

  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") %>%
    round(4)

  dimnames(ugos) <- list(lon = nc$dim$lon$vals,
                         lat = nc$dim$lat$vals,
                         time = nc$dim$time$vals)
  sla <- # slow!
    as_tibble(melt(ugos, value.name = "ugos"), row.names = NULL) %>%
    dplyr::mutate(vgos = as.vector(ncvar_get(nc, varid = "vgos")),
                  ugosa = as.vector(ncvar_get(nc, varid = "ugosa")),
                  vgosa = as.vector(ncvar_get(nc, varid = "vgosa")),
                  time = as.Date(time, origin = "1950-01-01 00:00:00")) %>%
    na.omit()
  nc_close(nc)

  # Mean and Eddy Kinetic Energy (MKE & EKE) function
  # per each month of the year
  ke <- sla %>%
    dplyr::mutate(time = floor_date(time, unit = "month"),
                  ugos = ugos * 100,
                  vgos = vgos * 100,
                  ugosa = ugosa * 100,
                  vgosa = vgosa * 100,
                  eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lat, lon, time) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2),
                     eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup()
  rm(sla)

  # Join the data and do correlation
  eke.cor.fun <- function(df) cor(df$eke, df$ex)
  mke.cor.fun <- function(df) cor(df$mke, df$ex)
  dat <- ex %>%
    dplyr::left_join(ke) %>%
    dplyr::filter(eke > 0) %>%
    na.omit() %>%
    dplyr::group_by(lon, lat) %>%
    nest() %>%
    dplyr::mutate(r.eke = map(data, eke.cor.fun)) %>%
    dplyr::mutate(r.mke = map(data, mke.cor.fun)) %>%
    unnest(r.eke, r.mke)

  return(dat)
}

# Apply the function
system.time(AC.cor <- uber.fun("AC"))
system.time(BC.cor <- uber.fun("BC"))
system.time(EAC.cor <- uber.fun("EAC"))
system.time(GS.cor <- uber.fun("GS"))
system.time(KC.cor <- uber.fun("KC"))

save(AC.cor, BC.cor, EAC.cor, GS.cor, KC.cor, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/correlations_v2.RData")

# Make a graph
gv.plot <- function(data, plot.parameters) {
  fig.a <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = r.eke)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
    theme_map() +
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
    plot.parameters

  fig.b <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = r.mke)) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
    theme_map() +
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
    plot.parameters

  fig <- ggarrange(fig.a,
                   fig.b,
                   ncol = 2, nrow = 1)
  return(fig)
}

AC.fig0x <- gv.plot(AC.cor, AC.layers)
BC.fig0x <- gv.plot(BC.cor, BC.layers)
EAC.fig0x <- gv.plot(EAC.cor, EAC.layers)
GS.fig0x <- gv.plot(GS.cor, GS.layers)
KC.fig0x <- gv.plot(KC.cor, KC.layers)

fig06 <- ggarrange(AC.fig0x,
                   BC.fig0x,
                   EAC.fig0x,
                   GS.fig0x,
                   KC.fig0x,
                   ncol = 1, nrow = 5)
ggplot2::ggsave("setup/publ_plots/Fig0x_correlations_v2.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)
