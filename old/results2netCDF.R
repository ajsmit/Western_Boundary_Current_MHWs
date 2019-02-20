library(ncdf4)

# get the data and turn implicit missing values into explicit missing values
# latitude is the fastest varying dimension
aviso_mean <- mAC.sl %>%
  complete(lon, lat)

# get x and y vectors (dimensions)
lon_vec <- unique(aviso_mean$lon); nlon <- length(lon_vec)
lat_vec <- unique(aviso_mean$lat); nlat <- length(lat_vec)

# define the dimensions
dimX <- ncdim_def("Longitude", "degrees_east", lon_vec)
dimY <- ncdim_def("Latitude", "degrees_north", lat_vec)

# define missing value
mv = -9999

# define the data to be written
# it is important that the fastest varying dimension is listed first in 'dim'
velocity_mean_var <- ncvar_def(name = "velocity_mean",
                               units="m/s",
                               dim = list(dimY, dimX),
                               missval = mv,
                               prec = "double",
                               longname = "Mean Aviso+ sea surface velocity (1993-01-01 to 2017-05-15)",
                               shuffle = TRUE,
                               compression = 5,
                               verbose = TRUE)
velocity_trend_var <- ncvar_def(name = "velocity_trend",
                                units = "m/s/decade",
                                dim = list(dimY, dimX),
                                missval = mv,
                                prec = "double",
                                longname = "Decadal trend in Aviso+ sea surface velocity (1993-01-01 to 2017-05-15)",
                                shuffle = TRUE,
                                compression = 5,
                                verbose = TRUE)
eke_trend_var <- ncvar_def(name = "eke_trend",
                           units = "m^2/s^2/decade",
                           dim = list(dimY, dimX),
                           missval = mv,
                           prec = "double",
                           longname = "Decadal trend in Aviso+ eddy kinetic energy (1993-01-01 to 2017-05-15)",
                           shuffle = TRUE,
                           compression = 5,
                           verbose = TRUE)
OISST_trend_var <- ncvar_def(name = "SST_trend",
                             units = "°C/decade",
                             dim = list(dimY, dimX),
                             missval = mv,
                             prec = "double",
                             longname = "Decadal trend in Optimally Interpolated SST (1981-09-01 to 2017-11-07)",
                             shuffle = TRUE,
                             compression = 5,
                             verbose = TRUE)

# create the NetCDF
# if you want a NetCDF4, explicitly add force_v4=T
nc <- nc_create("writevals.nc", vars = list(velocity_mean_var, velocity_trend_var, eke_trend_var, OISST_trend_var))
nc <- nc_create("ncfile.nc", vars = list(velocity_mean_var))

# although it is probably better to explicitely create a 2D array, it seems to work fine
# if this is not done; so, I just supply the data directly as (e.g.) `aviso_mean$velocity`
vel_dat <- array(aviso_mean$velocity, dim=c(nlon, nlat))

# write data to the NetCDF
ncvar_put(nc, velocity_mean_var, aviso_mean$velocity)
ncvar_put(nc, velocity_mean_var, vel_dat)
ncvar_put(nc, velocity_trend_var, data2d)
ncvar_put(nc, eke_trend_var, data2d)
ncvar_put(nc, OISST_trend_var, data2d)

# close the new file to finish writing
nc_close(nc)
