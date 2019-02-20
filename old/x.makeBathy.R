bathy <- as.matrix(read.table("/Users/ajsmit/spatial/bathy/GEBCO/gebco_sa.asc", sep = "", skip = 6))
info <- read.table("/Users/ajsmit/spatial/bathy/GEBCO/gebco_sa.asc", sep = "", nrows = 6) # Load in plot size parameters

xrng <- seq(from = info[3,2], by = info[5,2], length.out = info[1,2]) # Create x range integers
yrng <- seq(from = info[4,2], by = info[5,2], length.out = info[2,2]) # Create y range integers
xyrng <- expand.grid(xrng, yrng)
ht(xyrng)

bathy[bathy == info[6,2]] <- NA # Renames error values with NA
bathy[bathy > 0] <- NA # Removes any values above the surface of the water
bathy <- t(bathy[nrow(bathy):1,]) # Adds a column giving row names

library(ncdf4)
gebco.nc <- nc_open("/Users/ajsmit/spatial/bathy/GEBCO/gridone.nc")

library(reshape2)
library(tidyverse)
bathy.long <- as.tibble(melt(bathy))
ht(bathy.long)
length(xrng) * length(yrng)
nrow(bathy.long)
