# Plots of MKE and EKE ----------------------------------------------------

library(tidyverse)
library(RcppRoll)
library(ncdf4)
library(lubridate)
library(data.table)
library(colorspace)
library(scales)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 3)

# Functions: velocities, EKE, etc.
source("setup/functions.R")
source("setup/colour_palettes.R")

# Read in the region-specific components
source("setup/plot.layers.R")
csvInDir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/daily"

# Bathy data
bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC.bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC.bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC.bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS.bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC.bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))

# region <- "EAC"

ke.fun <- function(region) {
  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ke <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(2)

  # Prepare the data.table
  dimnames(ke) <- list(lon = nc$dim$lon$vals,
                       lat = nc$dim$lat$vals,
                       t = nc$dim$time$vals)
  ke <- as.data.table(reshape2::melt(ke, value.name = "ugos"))
  ke[, c("t", "vgos", "ugosa", "vgosa") := .( # slow
    as.Date(t, origin = "1950-01-01 00:00:00"),
    round(as.vector(ncvar_get(nc, varid = "vgos")) * 100, 2),
    round(as.vector(ncvar_get(nc, varid = "ugosa")) * 100, 2),
    round(as.vector(ncvar_get(nc, varid = "vgosa")) * 100, 2)
  )
  ]
  nc_close(nc)

  # calculate EKE relative to the 30-day running mean ("eke") and relative to
  # long-term mean ("eke_lt)
  out <- ke %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = 0.5 * ((ugos - roll_mean(ugos, n = 30, align = "center", fill = c(-999, -999, -999)))^2 +
                                 (vgos - roll_mean(vgos, n = 30, align = "center", fill = c(-999, -999, -999)))^2)) %>%
    na.omit() %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos)^2 + mean(ugos)^2),
                     eke = mean(eke)) %>%
    dplyr::ungroup()

  return(out)
}

# No need to run each time as the data have been saved already... load below...
AC.dat <- ke.fun("AC")
save(AC.dat, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/AC_ke_v4.RData")
rm(AC.dat)
BC.dat <- ke.fun("BC")
save(BC.dat, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/BC_ke_v4.RData")
rm(BC.dat)
EAC.dat <- ke.fun("EAC")
save(EAC.dat, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/EAC_ke_v4.RData")
rm(EAC.dat)
GS.dat <- ke.fun("GS")
save(GS.dat, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/GS_ke_v4.RData")
rm(GS.dat)
KC.dat <- ke.fun("KC")
save(KC.dat, file = "/Volumes/Benguela/spatial/processed/WBC/misc_results/KC_ke_v4.RData")
rm(KC.dat)


# Load data prepared above ------------------------------------------------

keDir <- "/Volumes/Benguela/spatial/processed/WBC/misc_results"
load(paste0(keDir, "/", "AC_ke_v4.RData"))
load(paste0(keDir, "/", "BC_ke_v4.RData"))
load(paste0(keDir, "/", "EAC_ke_v4.RData"))
load(paste0(keDir, "/", "KC_ke_v4.RData"))
load(paste0(keDir, "/", "GS_ke_v4.RData"))


# Agulhas Current ---------------------------------------------------------

load("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/AC-mask_polys.RData")
mke <- fortify(mask.list$mke)
eke <- fortify(mask.list$eke)
int <- fortify(mask.list$int)
AC.Fig00a <- ggplot(dplyr::filter(AC.dat, mke <= 7700), aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 3500, 7000)) +
  # scale_fill_continuous_sequential(palette = "Viridis", na.value = "#011789",
  #                                  rev = TRUE, breaks = c(500, 4100, 7700)) +
  stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + AC.layers

AC.Fig00b <- ggplot(dplyr::filter(AC.dat, eke_lt <= 3500), aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke_lt)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789",
                       breaks = c(0, 1400, 2800)) +
  # scale_fill_continuous_sequential(palette = "Viridis", rev = TRUE) +
  stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + AC.layers

AC.Fig00c <- ggplot(dplyr::filter(AC.dat, eke < 5000), aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789") +
  # scale_fill_continuous_sequential(palette = "Viridis", rev = TRUE) +
  stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + AC.layers


# Brazil Current ----------------------------------------------------------

load("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/BC-mask_polys.RData")
mke <- fortify(mask.list$mke)
eke <- fortify(mask.list$eke)
int <- fortify(mask.list$int)
BC.Fig00a <- ggplot(BC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 750, 1500)) +
  stat_contour(data = BC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + BC.layers

BC.Fig00b <- ggplot(BC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke_lt)) +
  scale_fill_gradientn(colours = col4, breaks = c(0, 800, 1600)) +
  stat_contour(data = BC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + BC.layers


# East Australian Current -------------------------------------------------

load("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/EAC-mask_polys.RData")
mke <- fortify(mask.list$mke)
eke <- fortify(mask.list$eke)
int <- fortify(mask.list$int)
EAC.Fig00a <- ggplot(EAC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 900, 1800)) +
  stat_contour(data = EAC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + EAC.layers

EAC.Fig00b <- ggplot(EAC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke_lt)) +
  scale_fill_gradientn(colours = col4, breaks = c(0, 2000, 4000)) +
  stat_contour(data = EAC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + EAC.layers

EAC.Fig00c <- ggplot(EAC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 900, 1800)) +
  stat_contour(data = EAC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + EAC.layers


# Gulf Stream -------------------------------------------------------------

load("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/GS-mask_polys.RData")
mke <- fortify(mask.list$mke)
eke <- fortify(mask.list$eke)
int <- fortify(mask.list$int)
GS.Fig00a <- ggplot(GS.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 4000, 8000)) +
  stat_contour(data = GS.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + GS.layers

GS.Fig00b <- ggplot(GS.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke_lt)) +
  scale_fill_gradientn(colours = col4, breaks = c(0, 1750, 3500)) +
  stat_contour(data = GS.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + GS.layers


# Kuroshio Current --------------------------------------------------------

load("/Volumes/GoogleDrive/My Drive/R/WBCs/masks/KC-mask_polys.RData")
mke <- fortify(mask.list$mke)
eke <- fortify(mask.list$eke)
int <- fortify(mask.list$int)
KC.Fig00a <- ggplot(KC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_gradientn(colours = col4, na.value = "#011789", breaks = c(0, 2500, 5000)) +
  stat_contour(data = KC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               colour = "red3", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + KC.layers

KC.Fig00b <- ggplot(KC.dat, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke_lt)) +
  scale_fill_gradientn(colours = col4, breaks = c(0, 1400, 2800)) +
  stat_contour(data = KC.bathy, aes(x = lon, y = lat, z = z),
               col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               colour = "navy", size = 0.4, fill = NA) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = expression((cm^2~s^{-2})),
                                frame.colour = "black",
                                frame.linewidth = 0.4,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(26, units = "mm"),
                                draw.ulim = F,
                                title.position = 'left',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  theme_map() + labs(x = NULL, y = NULL) + KC.layers

Fig00.mke <- ggarrange(AC.Fig00a,
                       BC.Fig00a,
                       EAC.Fig00a,
                       GS.Fig00a,
                       KC.Fig00a,
                       ncol = 1, nrow = 5, labels = "AUTO")
ggplot2::ggsave("publ_plots/Fig00_Aviso+_MKE_v2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)

Fig00.eke <- ggarrange(AC.Fig00b,
                       BC.Fig00b,
                       EAC.Fig00b,
                       GS.Fig00b,
                       KC.Fig00b,
                       ncol = 1, nrow = 5, labels = "AUTO")
ggplot2::ggsave("publ_plots/Fig00_Aviso+_EKE_v2.jpg",
                width = 3.5 * (1/3), height = 13 * (1/3), scale = 3.7)
