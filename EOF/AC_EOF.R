library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(heatwaveR)
library(grid)
library(vegan)
library(doMC); doMC::registerDoMC(cores = 3)


# Load stuff --------------------------------------------------------------

source("setup/plot.layers.R")
source("setup/functions.R")

# Bathy data

bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
AC.bathy <- fread(paste0(bathyDir, "/", "AC", "_bathy.csv"))
BC.bathy <- fread(paste0(bathyDir, "/", "BC", "_bathy.csv"))
EAC.bathy <- fread(paste0(bathyDir, "/", "EAC", "_bathy.csv"))
GS.bathy <- fread(paste0(bathyDir, "/", "GS", "_bathy.csv"))
KC.bathy <- fread(paste0(bathyDir, "/", "KC", "_bathy.csv"))


# Load kinetic energy data ------------------------------------------------

region <- "AC"

ke.fun <- function(region) {
  # Read in the Aviso data
  aviso.dir <- "/Volumes/Benguela/lustre/Aviso/global/delayed-time/grids/msl/all-sat-merged"
  aviso.dat <- "outfile_regridded.nc"
  nc <- nc_open(paste0(aviso.dir, "/", region, "/", aviso.dat))
  ugos <- ncvar_get(nc, varid = "ugos") * 100 %>%
    round(4)

  # Prepare the data.table
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

  out <- ke %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::mutate(eke = roll_mean(eke, n = 30, align = "center", fill = c(-999, -999, -999))) %>%
    dplyr::ungroup() %>%
    dplyr::select(lon, lat, t, eke) %>%
    dplyr::filter(eke > -999) %>%
    dplyr::filter(t >= "2015-07-01" & t <= "2018-06-30") %>%
    na.omit() %>%
    tidyr::spread(t, eke)

  return(out)
}

AC.spread <- ke.fun("AC")
AC.df <- na.omit(AC.spread)
AC.dates <- data.frame(date = fastDate(colnames(AC.df[, 3:ncol(AC.df)])))
AC.coords <- AC.df[, 1:2]
AC.mat <- AC.df[, 3:ncol(AC.df)]

AC.pca <- rda(AC.mat)
summary(AC.pca)

# Proportion of variation explained
(4.435e+08 / 1.053e+09) * 100 # 1 Eddy flux
(5.393e+07 / 1.053e+09) * 100 # 2
(3.881e+07 / 1.053e+09) * 100 # 3
(3.567e+07 / 1.053e+09) * 100 # 4 Meanders

AC.sites <- cbind(AC.coords, scores(AC.pca, choices = c(1:4), display = "sites"))
AC.times <- cbind(AC.dates, scores(AC.pca, choices = c(1:4), display = "species"))

ggplot(data = AC.times, aes(x = date, y = PC4)) +
  geom_line() + theme_bw()



load("/Users/ajsmit/Dropbox/R/WBCs/masks/AC-mask_polys.RData")
mke <- fortify(mask.list$mke90)
eke <- fortify(mask.list$eke90)
int <- fortify(mask.list$int90)

AC.Fig.PC1 <- ggplot(AC.sites, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = PC1)) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
  # stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
  #              col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = int, aes(long, lat, group = group),
               fill = NA, colour = "purple", size = 0.3) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               fill = NA, colour = "red3", size = 0.3) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               fill = NA, colour = "navy", size = 0.3) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = "PC1",
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

AC.Fig.PC2 <- ggplot(AC.sites, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = PC2)) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
  # stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
  #              col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = int, aes(long, lat, group = group),
               fill = NA, colour = "purple", size = 0.3) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               fill = NA, colour = "red3", size = 0.3) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               fill = NA, colour = "navy", size = 0.3) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = "PC2",
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

AC.Fig.PC3 <- ggplot(AC.sites, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = PC3)) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
  # stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
  #              col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = int, aes(long, lat, group = group),
               fill = NA, colour = "purple", size = 0.3) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               fill = NA, colour = "red3", size = 0.3) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               fill = NA, colour = "navy", size = 0.3) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = "PC3",
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

AC.Fig.PC4 <- ggplot(AC.sites, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = PC4)) +
  scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
  # stat_contour(data = AC.bathy, aes(x = lon, y = lat, z = z),
  #              col = "black", size = 0.15, breaks = c(-500, -1000, -2000)) +
  geom_polygon(data = int, aes(long, lat, group = group),
               fill = NA, colour = "purple", size = 0.3) +
  geom_polygon(data = mke, aes(long, lat, group = group),
               fill = NA, colour = "red3", size = 0.3) +
  geom_polygon(data = eke, aes(long, lat, group = group),
               fill = NA, colour = "navy", size = 0.3) +
  guides(alpha = "none",
         size = "none",
         fill = guide_colourbar(title = "PC4",
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

Fig00.AC.EOF <- ggarrange(
  AC.Fig.PC1,
  AC.Fig.PC2,
  AC.Fig.PC3,
  AC.Fig.PC4,
  ncol = 2,
  nrow = 2,
  labels = "AUTO"
)
ggplot2::ggsave("figures/Figxx_AC_EOF.jpg",
                width = 7.0 * (1/3), height = 1.733333, scale = 3.7)

# Remove seasonal cycle first

