# Get and prepare data OISST event data -----------------------------------

# Get and prepare Aviso data ----------------------------------------------

# Regridding performed with NCO
# First, make time into a record dimension (from a fixed dimension);
# do this for each of the the downloaded Aviso+ files:
# > ncks --mk_rec_dmn time infile.nc -O outfile.nc
# (the above are in the shell script in the Aviso+ netCDF directories, e.g.
# /Volumes/Benguela/Aviso/global/delayed-time/grids/msl/all-sat-merged/BC/unlimDim.sh)

# Now concatenate these files along the record dimension (i.e. time):
# > ncrcat infile????.nc outfile.nc
# (see the shell script for additional details)

# Make a 'template' file from the OISST data, which will be used as the grid to
# interpolate to. Since these data are global, they need to be subsetted first;
# this is termed creating a hyperslab in NCO speak. Thus, create a hyperslab that
# corresponds to each WBC region's spatial extent, e.g. for the Agulhas Current:
# > ncks -F -d lon,6.25,45. -d lat,-45.,-20. avhrr-only-v2.20180930.nc BC_template.nc
#
# The template file is just a daily netCDF files of the region of interest.

# Now interpolate the Aviso file to the OISST template file:
# > ncremap -i outfile.nc -d BC_template.nc -o outfile_regridded.nc

library(tidyverse)
library(data.table)
library(ncdf4)
library(lubridate)
library(colorspace)
library(scales)
library(RcppRoll)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 3)

# read in the region-specific components
source("setup/functions.R")
source("setup/plot_theme.R")

# Exceedence above threshold per month:
# this function calculates for each month and each pixel the
# mean MKE and EKE (from Aviso+) and mean mean MHW event intensity

region <- "BC"

clim.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
clim.base <- "-avhrr-only-v2.19810901-20180930_climatology.csv"
clim.dat <- paste0(region, clim.base)
clim <- fread(paste0(clim.dir, "/", clim.dat))
clim <- clim[lon >= 305 & lon <= 310, ] # trim to a section of lon that falls across area of maximal MHW activity
clim <- clim[lat >= -45 & lat <= -35, ]

clim[, t := fastDate(t)]
ex <- clim[t >= "1993-01-01" & t <= "2018-06-10", ] # trim to match start date of Aviso data
ex[, ex := roll_mean((temp - thresh), n = 30, align = "center", fill = c(-999, -999, -999)), by = .(lon, lat)] # calculate zonal means
ex2 <- ex[, list(ex = mean(temp - thresh)), by = list(lon = lon, t = t)]
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
ke <- as.data.table(reshape2::melt(ugos, value.name = "ugos"))
rm(ugos)
ke[, c("t", "vgos", "ugosa", "vgosa") := .( # slow
  as.Date(t, origin = "1950-01-01 00:00:00"),
  as.vector(ncvar_get(nc, varid = "vgos")) * 100,
  as.vector(ncvar_get(nc, varid = "ugosa")) * 100,
  as.vector(ncvar_get(nc, varid = "vgosa")) * 100
)
]
nc_close(nc)

ke <- ke[lon >= 305 & lon <= 310, ]
ke <- ke[lat >= -45 & lat <= -35, ]
ke[, eke := (0.5 * ((vgosa)^2 + (ugosa)^2))]
ke[, c("mke_roll", "eke_roll") := .(
  0.5 * ((vgos - roll_mean(vgos, n = 30, align = "center", fill = c(-999, -999, -999)))^2 + (ugos - roll_mean(ugos, n = 30, align = "center", fill = c(-999, -999, -999)))^2),
  roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))
), by = .(lon, lat)
]
ke <- na.omit(ke)
ke[, c("ugos", "vgos", "ugosa", "vgosa") := NULL]

ke2 <- ke[, list(eke = mean(eke),
                 mke_roll = mean(mke_roll),
                 eke_roll = mean(eke_roll)), by = list(lon = lon, t = t)]

# Join the data and do correlation
# Join the tables
setkey(ex2, lon, t)
setkey(ke2, lon, t)
dat <- ex2[ke2, nomatch = 0]
dat <- dat[eke_roll > 0, ]
dat <- na.omit(dat)

mke.cor.fun <- function(df) cor(df$mke_roll, df$ex)
eke.cor.fun <- function(df) cor(df$eke_roll, df$ex)
dat.r <- dat %>%
  dplyr::group_by(lon) %>%
  nest() %>%
  dplyr::mutate(MKE = map(data, mke.cor.fun)) %>%
  dplyr::mutate(EKE = map(data, eke.cor.fun)) %>%
  unnest(c(MKE, EKE)) %>%
  dplyr::select(-data)

dat.r.long <- gather(dat.r, key = "ke.type", value = "value", MKE:EKE)

save(dat.r.long, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/BC_correlations_v4.RData")

# Make a graph
BC.corr.plt <- ggplot(data = dat.r.long, aes(x = lon, group = ke.type)) +
  geom_line(aes(y = value, colour = ke.type), size = 0.7) +
  scale_color_manual(values = c("navy", "red3")) +
  labs(x = "Longitude (°E)",
       y = expression(paste("Correlation (", italic(r), ")")),
       colour = "Kinetic\nenergy") +
  theme_ajs()
BC.corr.plt

save(BC.corr.plt, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/BC_corr_plt_v1.RData")

ggplot2::ggsave("publ_plots/Fig0x_BC_meridional.jpg",
                width = 7.0, height = 2.5, scale = 2.5, units = "cm")
