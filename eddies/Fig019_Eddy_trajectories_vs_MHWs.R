library(tidyverse)
library(ggpubr)
library(lubridate)
library(data.table)
library(colorspace)
library(OneR) # for binning


# Miscellaneous things to load --------------------------------------------

# Read in the region-specific components and functions
source("setup/functions.R")
source("setup/regionDefinition.R")
source("setup/plot.layers.R")
source("setup/colour_palettes.R")

# region <- "AC"
# plot.parameters <- AC.layers

# The Gulf Stream and Brazil don't plot here because it gets the
# axes coords from plot.layers, which doesn't match up with the
# way you use the East/West coords in your data...

eddy_plot <- function(eddies, region, plot.parameters, ev.date) {
  # bathy data
  bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
  bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))
  load(paste0("masks/", region, "-mask_polys.RData"))
  load(paste0("eddies/eddy_vs_MHW_",region,".Rdata"))

  # prep
  mke <- fortify(mask.list$mke) %>%
    dplyr::mutate(long = ifelse(long > 180, long - 360, long))
  eddy_MHW_res_daily <- eddy_MHW_res$daily %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lon_eddy = ifelse(lon_eddy > 180, lon_eddy - 360, lon_eddy))
  coords <- bbox %>%
    dplyr::select({{region}}) %>%
    dplyr::mutate(coord = row.names(.)) %>%
    spread(key = coord, value = {{region}}) %>%
    dplyr::mutate(lonmin = ifelse(lonmin > 180, lonmin - 360, lonmin),
           lonmax = ifelse(lonmax > 180, lonmax - 360, lonmax))

  eddy_plt <- ggplot(data = eddy_MHW_res_daily, aes(x = lon_eddy, y = lat_eddy)) +
    geom_point(aes(colour = pixel_no_MKE_MHW_prop, alpha = pixel_no_MKE_MHW_prop), size = 0.1) +
    scale_color_gradient(low = "white", high = "navy") +
    stat_contour(data = bathy, aes(x = lon, y = lat, z = z),
                 col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
    geom_polygon(data = mke, aes(long, lat, group = group),
                 colour = "red3", size = 0.4, fill = NA) +
    scale_alpha_continuous(range = c(0.1, 1.0)) +
    guides(alpha = "none",
           size = "none",
           colour = "none") +
    theme_map() +
    plot.parameters

  # ggplot2::ggsave("eddies/Fig018_Eddy_trajectories_AC.jpg",
  #                 width = 3.5 * (1/3), height = 2.6 * (1/3), scale = 3.7)

  return(eddy_plt)
}

AC.Fig019 <- eddy_plot(eddies, "AC", AC.layers, AC.ev)
BC.Fig019 <- eddy_plot(eddies, "BC", BC.layers, BC.ev)
EAC.Fig019 <- eddy_plot(eddies, "EAC", EAC.layers, EAC.ev)
GS.Fig019 <- eddy_plot(eddies, "GS", GS.layers, GS.ev)
KC.Fig019 <- eddy_plot(eddies, "KC", KC.layers, KC.ev)

Fig019 <- ggarrange(AC.Fig019,
                    BC.Fig019,
                    EAC.Fig019,
                    GS.Fig019,
                    KC.Fig019,
                    ncol = 1, nrow = 5, labels = "AUTO")
ggplot2::ggsave("eddies/Fig019_Eddy_trajectories_MHW_vs_eddies.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
