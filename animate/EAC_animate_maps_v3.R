# devtools::install_github("tidyverse/dplyr@rc_0.8.0") # for new dplyr function, group_split
library(tidyverse)
library(ggpubr)
library(data.table)
library(colorspace)
library(RcppRoll)


# Setup -------------------------------------------------------------------

source("setup/functions.R")
source("setup/plot.layers.R")
source("setup/regionDefinition.R")


# Define session ----------------------------------------------------------

region <- "EAC"
plot.parameters <- EAC.layers


# Load data function ------------------------------------------------------

load_data <- function(region) {

  # Load SST climatology
  sst_dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
  clim <- fread(paste0(sst_dir, "/", region, "-avhrr-only-v2.19810901-20180930_climatology.csv"),
                colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7),
                col.names = c("lon", "lat", "doy", "t", "temp", "seas", "thresh"))
  clim[, c("t", "ex", "anom") := .(
    fastDate(t),
    temp - thresh,
    temp - seas
  )]
  clim[, c("doy", "temp", "seas", "thresh") := NULL] # drop unused columns
  clim[, c("ex", "anom") := .(
    roll_mean(ex, n = 5, align = "center"),
    roll_mean(anom, n = 5, align = "center")
    ),
    by = c("lon", "lat")]
  clim[ex > 0, ev := 1]
  clim[ex <= 0, ev := 0]

  # Load Aviso+ data
  aviso_dir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/daily"
  geo <- fread(paste0(aviso_dir, "/", region, "-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
               colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:9),
               col.names = c("lon", "lat", "t", "ugos", "vgos", "ugosa", "vgosa", "sla", "adt"))
  geo[, c('t', 'velocity', 'eke', 'arrow_size') := .(
    fastDate(t),
    sqrt((vgos^2) + (ugos^2)),
    0.5 * (ugos^2 + vgos^2),
    ((abs(ugos * vgos) / max(abs(ugos * vgos))) + 0.3) / 6
  )]
  geo[, c("ugos", "vgos", "ugosa", "vgosa", "sla", "adt") := NULL]

  # Join the tables
  setkey(clim, lon, lat, t)
  setkey(geo, lon, lat, t)
  dat <- clim[geo, nomatch = 0]
  dat[, no := seq(1:.N), by = c("lon", "lat")]
  rm(clim); rm(geo)

  return(dat)
}

dat <- load_data(region)
plt.dat <- dat[t >= "2015-07-01" & t <= "2018-06-30", ]
plt.dat[, cum := cumsum(ev), by = c("lon", "lat")]
summary(plt.dat)
test.dat <- dat[t >= "2015-07-01" & t <= "2015-07-31", ]
test.dat[, cum := cumsum(ev), by = c("lon", "lat")]
summary(test.dat)


# The plot function -------------------------------------------------------

combo_plot <- function(data, plot.parameters, region) {
  maskDir <- "/Users/ajsmit/Dropbox/R/WBCs/masks/"
  load(paste0(maskDir, region, "-mask_polys.RData"))
  mke90 <- fortify(mask.list$mke90)
  eke90 <- fortify(mask.list$eke90)
  int90 <- fortify(mask.list$int90)
  top_right <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = ex)) +
    # geom_polygon(data = mke90, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    # geom_polygon(data = int90, aes(long, lat, group = group),
    #              fill = alpha("#c83c7e", 1.0), colour = NA) +
    # geom_polygon(data = eke75, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "navy", size = 0.3,
    #              linetype = "dashed") +
    geom_polygon(data = eke90, aes(long, lat, group = group),
                 fill = NA, colour = "navy", size = 0.2) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE) +
    guides(fill = guide_colourbar(title = "[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(31, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() +
    plot.parameters + labs(subtitle = "SST exceedence over threshold", title = " ")

  top_left <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = eke)) +
    # geom_polygon(data = mke90, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    # geom_polygon(data = int90, aes(long, lat, group = group),
    #              fill = alpha("#c83c7e", 1.0), colour = NA) +
    # geom_polygon(data = eke75, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "navy", size = 0.3,
    #              linetype = "dashed") +
    geom_polygon(data = eke90, aes(long, lat, group = group),
                 fill = NA, colour = "navy", size = 0.2) +
    scale_fill_continuous_sequential(palette = "Inferno", na.value = "#011789", rev = TRUE) +
    guides(fill = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters + labs(subtitle = "Gestrophic velocity", title = as.character(unique(data$t)))

  # llim <- min(data$cum)
  # ulim <- max(data$cum)

  bottom_left <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = anom)) +
    # geom_polygon(data = mke90, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "red3", size = 0.3) +
    # geom_polygon(data = int90, aes(long, lat, group = group),
    #              fill = alpha("#c83c7e", 1.0), colour = NA) +
    # geom_polygon(data = eke75, aes(long, lat, group = group),
    #              fill = alpha("gray50", 0.5), colour = "navy", size = 0.3,
    #              linetype = "dashed") +
    geom_polygon(data = eke90, aes(long, lat, group = group),
                 fill = NA, colour = "navy", size = 0.2) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", na.value = "#011789", rev = FALSE,
                                    limits = c(-4, 5.5)) +
    guides(fill = guide_colourbar(title = "[°C]",
                                  frame.colour = "black",
                                  frame.linewidth = 0.4,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(31, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'left',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters + labs(subtitle = "SST anomaly", title = " ")

  bottom_right <- ggplot(data, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = cum)) +
    scale_fill_continuous_sequential(palette = "Purples 3", na.value = "#011789", rev = TRUE,
                                    limits = c(0, 665)) +
    guides(fill = "none") +
    theme_map() + labs(x = NULL, y = NULL) +
    plot.parameters + labs(subtitle = "Event count", title = " ")

  fig <- ggarrange(top_left, top_right, bottom_left, bottom_right, ncol = 2, nrow = 2)

  ggsave(filename = paste0("~/temp/", region, "/", data[1, "no"], "-", region, "-",
                           as.character(unique(data$t)), "-anim_sequence.jpg"),
         plot = fig, width = 2.0626, height = 0.858 * 2, scale = 3.7)
}

# Apply the function! -----------------------------------------------------

# fun1 <- function(data, plot.parameters, region) {
#   combo_plot(data, plot.parameters, region)
# }

plt.dat %>%
  dplyr::group_split(t) %>%
  purrr::map(.f = combo_plot, plot.parameters, region)


# Animate... --------------------------------------------------------------

# animate in the terminal using ffmpeg:
ffmpeg \
-framerate 5 \
-pattern_type glob -i '*.jpg' \
-vf scale=1280:-2 \
EAC_out.mp4 \
;
