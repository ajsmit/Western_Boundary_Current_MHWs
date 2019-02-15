# Plots of MKE and EKE ----------------------------------------------------

library(tidyverse)
library(ggplot2)
library(lubridate)
library(data.table)
library(colorspace)
library(scales)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# Read in the region-specific components
source("setup/plot.layers.R")
csvInDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"

# Mean Kinetic Energy (MKE) function
mke.fun <- function(dat) {
  mke <- dat %>%
    dplyr::select(-ugosa, -vgosa) %>%
    dplyr::mutate(ugos = ugos * 100,
                  vgos = vgos * 100,
                  year = lubridate::year(fastDate(time)),
                  month = lubridate::month(fastDate(time))) %>%
    dplyr::group_by(lat, lon) %>%
    dplyr::summarise(mke = 0.5 * (mean(vgos, na.rm = TRUE)^2 + mean(ugos, na.rm = TRUE)^2)) %>%
    dplyr::ungroup()
  return(mke)
}

# Eddy Kinetic Energy (EKE) function
eke.fun <- function(dat) {
  eke <- dat %>%
    dplyr::mutate(ugos = ugos * 100,
                  vgos = vgos * 100,
                  ugosa = ugosa * 100,
                  vgosa = vgosa * 100,
                  year = lubridate::year(fastDate(time)),
                  month = lubridate::month(fastDate(time))) %>%
    dplyr::mutate(eke = 0.5 * ((vgosa)^2 + (ugosa)^2)) %>%
    dplyr::group_by(lon, lat, month) %>%
    dplyr::summarise(eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(eke = mean(eke, na.rm = TRUE)) %>%
    dplyr::ungroup()
  return(eke)
}


# Agulhas Current ---------------------------------------------------------

AC.file <- "AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
AC.sl <- as.tibble(fread(paste0(csvInDir, "/", AC.file), verbose = FALSE))
colnames(AC.sl) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")

AC.mke <- mke.fun(AC.sl)
AC.eke <- eke.fun(AC.sl)
rm(AC.sl)

plt1a <- ggplot(AC.mke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(500, 4100, 7700)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

plt1b <- ggplot(AC.eke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_continuous_sequential(palette = "Inferno", rev = TRUE, breaks = c(200, 1500, 3000)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

BC.file <- "BC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
BC.sl <- as.tibble(fread(paste0(csvInDir, "/", BC.file), verbose = FALSE))
colnames(BC.sl) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")

BC.mke <- mke.fun(BC.sl)
BC.eke <- eke.fun(BC.sl)
rm(BC.sl)

plt2a <- ggplot(BC.mke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(200, 800, 1600)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

plt2b <- ggplot(BC.eke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_continuous_sequential(palette = "Inferno", rev = TRUE, breaks = c(200, 800, 1600)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

EAC.file <- "EAC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
EAC.sl <- as.tibble(fread(paste0(csvInDir, "/", EAC.file), verbose = FALSE))
colnames(EAC.sl) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")
EAC.sl <- EAC.sl %>%
  dplyr::filter(lat <= -15)

EAC.mke <- mke.fun(EAC.sl)
EAC.eke <- eke.fun(EAC.sl)
rm(EAC.sl)

plt3a <- ggplot(EAC.mke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(200, 900, 1900)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

plt3b <- ggplot(EAC.eke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_continuous_sequential(palette = "Inferno", rev = TRUE, breaks = c(400, 2100, 4200)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

GS.file <- "GS-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
GS.sl <- as.tibble(fread(paste0(csvInDir, "/", GS.file), verbose = FALSE))
colnames(GS.sl) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")

GS.mke <- mke.fun(GS.sl)
GS.eke <- eke.fun(GS.sl)
rm(GS.sl)

plt4a <- ggplot(GS.mke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(1000, 4900, 9800)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

plt4b <- ggplot(GS.eke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_continuous_sequential(palette = "Inferno", rev = TRUE, breaks = c(200, 1800, 3600)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

KC.file <- "KC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
KC.sl <- as.tibble(fread(paste0(csvInDir, "/", KC.file), verbose = FALSE))
colnames(KC.sl) <- c("lon", "lat", "time", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt")

KC.mke <- mke.fun(KC.sl)
KC.eke <- eke.fun(KC.sl)
rm(KC.sl)

plt5a <- ggplot(KC.mke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = mke)) +
  scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE, breaks = c(300, 3000, 6000)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

plt5b <- ggplot(KC.eke, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = eke)) +
  scale_fill_continuous_sequential(palette = "Inferno", rev = TRUE, breaks = c(200, 1600, 3000)) +
  guides(fill = guide_colourbar(title = expression((cm^2~s^{-2})),
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

fig.ke <- ggarrange(plt1a, plt1b,
                    plt2a, plt2b,
                    plt3a, plt3b,
                    plt4a, plt4b,
                    plt5a, plt5b,
                    ncol = 2, nrow = 5)
ggplot2::ggsave("setup/publ_plots/Figure02_mke_eke.jpg",
                width = 7.0 * (1/3), height = 13 * (1/3), scale = 3.7)
