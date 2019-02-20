# Questions and notes
# 1. how do we determine how often the current loops and pinches off from the average current path?
# Southern hemisphere eddies:


# setup -------------------------------------------------------------------

library(data.table)
library(fasttime)
library(heatwaveR)
library(plyr)
library(tidyverse)
library(lubridate)
library(scales)
library(ggpubr)
library(doMC); doMC::registerDoMC(cores = 4)

# Functions: velocities, EKE, etc.
source("setup/functions.R")

# Read in the region-specific components
source("setup/plot.layers.R")

region  <- "AC"

# The bathy data
source("setup/regionDefinition.R")
bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
bathy <- fread(paste0(bathyDir, "/", region, "_bathy.csv"))

# Various paths
subDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/subregions"
OISSTDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/climatology"
evDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/MHWs"
dailyDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
AvisoDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"


# recode factors ----------------------------------------------------------

# recode factors and create mean of replicate 'in' and 'out' regions within each current, etc.
subs.in <- as.tibble(fread(paste0(subDir, "/all-sub-avhrr-only-v2-19810904-20180928.csv")))
subs.raw <- subs.in %>%
  mutate(t = fastDate(t)) %>%
  mutate(location = fct_recode(ID,
                               "current" = "sub1",
                               "current" = "sub2",
                               "current" = "sub3",
                               "reference" = "sub4",
                               "reference" = "sub5",
                               "reference" = "sub6",
                               "event" = "sub7",
                               "event" = "sub8",
                               "event" = "sub9")) %>%
  group_by(t, region, location) %>%
  dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(region == "AC")
rm(subs.in)


# make climatology and detect events --------------------------------------

# events and climatologies: create climatology and create events for each
# of the regions and 'in' and 'out' groups

# CLIMATOLOGY function
ts2clm.grid <- function(data) {
  out <- ts2clm(data, x = t, y = temp,
                climatologyPeriod = c("1983-01-01", "2012-12-31"),
                robust = FALSE, maxPadLength = 3, windowHalfWidth = 5,
                pctile = 90, smoothPercentile = TRUE,
                smoothPercentileWidth = 31, clmOnly = FALSE)
  return(out)
}

climatology <- as.tibble(ddply(subs.raw, .(location), ts2clm.grid, .parallel = TRUE, .progress = "text"))

# DETECT FUNCTION using heatwaveR
detect.grid <- function(dat) {
  out <- heatwaveR::detect_event(dat)
  event <- out$event
  return(event)
}

events <- as.tibble(ddply(climatology, .(location), detect.grid, .parallel = TRUE, .progress = "text"))

detect.grid.long <- function(dat) {
  out <- heatwaveR::detect_event(dat, protoEvents = TRUE)
  return(out)
}

events.long <- as.tibble(ddply(climatology, .(location), detect.grid.long, .parallel = TRUE, .progress = "text"))


# find the largest events for later plotting and analysis -----------------

# largest Agulhas Current events
top.ev <- events %>%
  filter(location == "event") %>%
  dplyr::arrange(desc(intensity_mean)) %>%
  head(5) %>%
  dplyr::select(date_start:date_end) %>%
  as.matrix()
# focus on events 2, 4, and 5 (all > 2000-01-01)

# load gridded OISST data -------------------------------------------------

# load the temp and climatology (this is a long climatology, i.e. includes
# daily SST for the full time period)
sst.clim <- fread(paste0(OISSTDir, "/AC-avhrr-only-v2.19810901-20180930_climatology.csv"),
                  colClasses = list("numeric" = 1:3, "myDate" = 4, "numeric" = 5:7))

# load the events prepared from the above climatology
ev <- as.tibble(fread(paste0(evDir, "/AC-avhrr-only-v2.19810901-20180930_events.csv")))


# plot climatology and daily SST sections ---------------------------------

# extract data along -40°S (10°E to 25°E)
# event 2
sst.clim.ev2.ts <- sst.clim %>%
  dplyr::filter(t == top.ev[2,2]) %>%
  dplyr::filter(lon >= 10 & lon <= 25) %>%
  dplyr::filter(lat == -40.125)

# event 4
sst.clim.ev4.ts <- sst.clim %>%
  dplyr::filter(t == top.ev[4,2]) %>%
  dplyr::filter(lon >= 10 & lon <= 25) %>%
  dplyr::filter(lat == -40.125)

# event 5
sst.clim.ev5.ts <- sst.clim %>%
  dplyr::filter(t == top.ev[5,2]) %>%
  dplyr::filter(lon >= 10 & lon <= 25) %>%
  dplyr::filter(lat == -40.125)

# plot SST transect
ggplot(data = sst.clim.ev2.ts, aes(x = lon, y = temp - seas)) +
  geom_line(size = 0.8, colour = "blue") +
  geom_line(data = sst.clim.ev4.ts, aes(x = lon, y = temp - seas), size = 0.8, colour = "blue3") +
  geom_line(data = sst.clim.ev5.ts, aes(x = lon, y = temp - seas), size = 0.8, colour = "navy") +
  geom_line(aes(y = seas - mean(seas)), colour = "green3", size = 0.8)
# make a similar plot as above, but outside of heatwave events?


# make flame plots --------------------------------------------------------

plt1 <- flame.plot(climatology, "AC", "reference")
plt2 <- flame.plot(climatology, "AC", "current")
plt3 <- flame.plot(climatology, "AC", "event")

ggarrange(plt1,
          plt2,
          plt3,
          ncol = 1, nrow = 3)
# ggplot2::ggsave("setup/regional_analyses/Agulhas_Current/AC_OISST_flame_2010-2018.png",
#                 width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# get the Aviso+ data -----------------------------------------------------

vel_AC <- as.tibble(fread(paste0(AvisoDir, "/AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                          col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))


# plot current/EKE on peak of selected events -----------------------------

# events limited to the Aviso+ era
# plot current/EKE on peak of second most intense mean event
AC_ev2 <- vel_AC %>%
  filter(time == top.ev[2, 2]) %>%
  mutate(velocity = sqrt((v^2) + (u^2)),
         eke = 0.5 * (u^2 + v^2),
         arrow_size = ((abs(u * v) / max(abs(u * v))) + 0.3) / 6)

# plot current/EKE on peak of fourth most intense mean event
AC_ev4 <- vel_AC %>%
  filter(time == top.ev[4, 2]) %>%
  mutate(velocity = sqrt((v^2) + (u^2)),
         eke = 0.5 * (u^2 + v^2),
         arrow_size = ((abs(u * v) / max(abs(u * v))) + 0.3) / 6)

# plot current/EKE on peak of fifth most intense mean event
AC_ev5 <- vel_AC %>%
  filter(time == top.ev[5, 2]) %>%
  mutate(velocity = sqrt((v^2) + (u^2)),
         eke = 0.5 * (u^2 + v^2),
         arrow_size = ((abs(u * v) / max(abs(u * v))) + 0.3) / 6)

vel.plot(AC_ev2, bathy, AC.layers) + labs(subtitle = paste(top.ev[2,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_peak_mean_intensity_map_", top.ev[2,2], ".png"),
#                 width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

vel.plot(AC_ev4, bathy, AC.layers) + labs(subtitle = paste(top.ev[4,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_peak_mean_intensity_map_", top.ev[4,2], ".png"),
#                 width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

vel.plot(AC_ev5, bathy, AC.layers) + labs(subtitle = paste(top.ev[5,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_peak_mean_intensity_map_", top.ev[5,2], ".png"),
#                 width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# setup and make Hovmöller plots ------------------------------------------

hov.ev2 <- sst.clim %>%
  dplyr::filter(lon >= 10 & lon <= 22.5) %>%
  dplyr::filter(lat >= -41 & lat <= -40) %>%
  dplyr::filter(t >= as.Date(top.ev[2, 1]) %m-% months(1) &
                  t <= as.Date(top.ev[2, 3]) %m+% months(1)) %>%
  dplyr::mutate(t = fastDate(t)) %>%
  dplyr::group_by(lon, t) %>%
  dplyr::summarise(temp = mean(temp, na.rm = TRUE),
                   seas = mean(seas, na.rm = TRUE),
                   thresh = mean(thresh, na.rm = TRUE)) %>%
  dplyr::ungroup()

hov.ev4 <- sst.clim %>%
  dplyr::filter(lon >= 10 & lon <= 22.5) %>%
  dplyr::filter(lat >= -41 & lat <= -40) %>%
  dplyr::filter(t >= as.Date(top.ev[4, 1]) %m-% months(1) &
                  t <= as.Date(top.ev[4, 3]) %m+% months(1)) %>%
  dplyr::mutate(t = fastDate(t)) %>%
  dplyr::group_by(lon, t) %>%
  dplyr::summarise(temp = mean(temp, na.rm = TRUE),
                   seas = mean(seas, na.rm = TRUE),
                   thresh = mean(thresh, na.rm = TRUE)) %>%
  dplyr::ungroup()

hov.ev5 <- sst.clim %>%
  dplyr::filter(lon >= 10 & lon <= 22.5) %>%
  dplyr::filter(lat >= -41 & lat <= -40) %>%
  dplyr::filter(t >= as.Date(top.ev[5, 1]) %m-% months(1) &
                  t <= as.Date(top.ev[5, 3]) %m+% months(1)) %>%
  dplyr::mutate(t = fastDate(t)) %>%
  dplyr::group_by(lon, t) %>%
  dplyr::summarise(temp = mean(temp, na.rm = TRUE),
                   seas = mean(seas, na.rm = TRUE),
                   thresh = mean(thresh, na.rm = TRUE)) %>%
  dplyr::ungroup()

hov.plt <- function(data, start, peak, end) {
  ggplot(data, aes(x = lon, y = t)) +
  geom_raster(aes(fill = temp - seas), interpolate = FALSE) +
  scale_fill_gradientn(name = "SST anomaly (°C)",
                       colours = col1, limits = c(-8.5, 8.5), position = "bottom") +
  guides(fill = guide_colourbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                barheight = unit(2, units = "mm"),
                                barwidth = unit(100, units = "mm"),
                                draw.ulim = F,
                                title.position = 'top',
                                title.hjust = 0.5,
                                label.hjust = 0.5)) +
  geom_hline(aes(yintercept = as.Date(start)), linetype = "dashed") +
  geom_hline(aes(yintercept = as.Date(end)), linetype = "dashed") +
  geom_hline(aes(yintercept = as.Date(peak))) +
  geom_vline(aes(xintercept = 15), alpha = 0.8, linetype = "dotted") +
  geom_vline(aes(xintercept = 16), alpha = 0.8, linetype = "dotted") +
  theme(legend.position = "bottom") +
  labs(title = "Agulhas Current", subtitle = paste0(start, " to ", end),
       y = NULL, x = "Longitude (°E)")
}

# 2014-03-01 2014-03-24 2014-04-29
hov.plt(hov.ev2, top.ev[2, 1], top.ev[2, 2], top.ev[2, 3]) + labs(subtitle = paste(top.ev[2,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_Hovmöller_", top.ev[2,2], ".png"),
#                 width = 5, height = 3.5, units = "cm", dpi = "retina", scale = 2.8)

# 2014-11-23 2014-12-21 2015-02-24
hov.plt(hov.ev4, top.ev[4, 1], top.ev[4, 2], top.ev[4, 3]) + labs(subtitle = paste(top.ev[4,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_Hovmöller_", top.ev[4,2], ".png"),
#                 width = 5, height = 3.5, units = "cm", dpi = "retina", scale = 2.8)

# 2018-03-05 2018-04-05 2018-05-02
hov.plt(hov.ev5, top.ev[5, 1], top.ev[5, 2], top.ev[5, 3]) + labs(subtitle = paste(top.ev[5,2]))
# ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_Hovmöller_", top.ev[5,2], ".png"),
#                 width = 5, height = 3.5, units = "cm", dpi = "retina", scale = 2.8)


# plot OISST anomalies on top event days ----------------------------------

sst.ev2 <- sst.clim %>%
  dplyr::filter(t == top.ev[2, 2])

sst.ev4 <- sst.clim %>%
  dplyr::filter(t == top.ev[4, 2])

sst.ev5 <- sst.clim %>%
  dplyr::filter(t == top.ev[5, 2])

sst.plot <- function(sst.dat, bathy.dat, plot.parameters) {
  fig <- ggplot(sst.dat, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = (temp - seas)), interpolate = FALSE) +
    scale_fill_gradientn(colours = col1, limits = c(-12.4, 12.4)) +
    stat_contour(data = bathy.dat, aes(x = lon, y = lat, z = z, alpha = ..level..),
                 col = "black", size = 0.15, binwidth = 500) +
    scale_alpha_continuous(range = c(0.2, 0.9)) +
    guides(alpha = "none",
           fill = guide_colourbar(title = "SST anomaly (°C)",
                                  frame.colour = "black",
                                  frame.linewidth = 0.5,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(100, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'top',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5)) +
    theme_map() +
    plot.parameters
  return(fig)
}

sst.plot(sst.ev2, bathy, AC.layers) + labs(subtitle = paste(top.ev[2,2]))
ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_SST_anom_map_", top.ev[2,2], ".png"),
                width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

sst.plot(sst.ev4, bathy, AC.layers) + labs(subtitle = paste(top.ev[4,2]))
ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_SST_anom_map_", top.ev[4,2], ".png"),
                width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

sst.plot(sst.ev5, bathy, AC.layers) + labs(subtitle = paste(top.ev[5,2]))
ggplot2::ggsave(paste0("setup/regional_analyses/Agulhas_Current/AC_SST_anom_map_", top.ev[5,2], ".png"),
                width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# read in monthly MODIS chl-a ---------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/subregions"
chl <- fread(paste0(csvDir, "/all-sub-MODIS-Aqua-20020705-20181118.csv"))

cols <- rep(c("red3", "blue3", "green4"), 3)

# recode the ID variable
chlor_a <- chl %>%
  dplyr::filter(region == "AC") %>%
  dplyr::filter(time >= "2014-01-01") %>%
  mutate(time = fastDate(time)) %>%
  mutate(location = fct_recode(ID,
                               "current" = "sub1",
                               "current" = "sub2",
                               "current" = "sub3",
                               "reference" = "sub4",
                               "reference" = "sub5",
                               "reference" = "sub6",
                               "event" = "sub7",
                               "event" = "sub8",
                               "event" = "sub9")) %>%
  dplyr::select(-region)

ggplot(chlor_a, aes(x = time, y = chlor_a, group = location)) +
  geom_line(aes(colour = ID, group = ID), alpha = 0.7, size = 0.4) +
  scale_color_manual(values = cols) +
  scale_y_continuous(trans = 'log10') +
  guides(alpha = "none", colour = "none") +
  labs(title = "MODIS Aqua 8-day chlorophyll-a",
       subtitle = "2014 to 2018",
       x = NULL,
       y = "Chlorophyll-a [mg/m3]") +
  facet_wrap(~ location, nrow = 3) + theme_minimal()

# from: https://stackoverflow.com/questions/9500114/find-which-season-a-particular-date-belongs-to
# make lookup table for finding a season that a particular date belongs to;
# accounting for leap years...
d <- function(month_day) which(lut$month_day == month_day)
lut <- data.frame(all_dates = as.POSIXct("2012-1-1") + ((0:365) * 3600 * 24),
                  season = NA)
lut <- within(lut, { month_day = strftime(all_dates, "%b-%d") })
lut[c(d("Jan-01"):d("Mar-20"), d("Dec-21"):d("Dec-31")), "season"] = "winter"
lut[c(d("Mar-21"):d("Jun-20")), "season"] = "spring"
lut[c(d("Jun-21"):d("Sep-20")), "season"] = "summer"
lut[c(d("Sep-21"):d("Dec-20")), "season"] = "autumn"
rownames(lut) = lut$month_day

chlor_a <- within(chlor_a, {
  season = lut[strftime(time, "%b-%d"), "season"]
 })

# ANOVA comparing the chl-a concs. between times with events and times
# without events, separately for 'event', 'current' and 'reference' areas...
events.present <- events.long %>%
  dplyr::filter(event == TRUE) %>%
  dplyr::select(t)
events.absent <- events.long %>%
  dplyr::filter(event == FALSE) %>%
  dplyr::select(t)

aov.fun <- function(data) {
  require(nlme)
  A <- data %>%
    dplyr::filter(time %in% events.present$t) %>%
    dplyr::mutate(presence = "yes")
  B <- data %>%
    dplyr::filter(time %in% events.absent$t) %>%
    dplyr::mutate(presence = "no")
  dat <- rbind(A, B)
  # return(t.test(A$chlor_a, B$chlor_a, na.rm = TRUE))
  mod <- aov(chlor_a ~ presence * season, data = dat)
  out <- summary(mod)
  return(out)
}

dlply(chlor_a, .(location), aov.fun)

# for "current", "event", and "reference" areas, extract the chl-a concentrations during times corresponding to MHWs and compare graphically with chl-a concentrations during times not coinciding with MHWs...
box.fun <- function(data) {
  require(nlme)
  A <- data %>%
    dplyr::filter(time %in% events.present$t) %>%
    dplyr::mutate(presence = "yes")
  B <- data %>%
    dplyr::filter(time %in% events.absent$t) %>%
    dplyr::mutate(presence = "no")
  dat <- rbind(A, B)

  library(ggpubr)
  out <- ggplot(dat, aes(x = season, y = chlor_a)) +
    geom_boxplot(aes(colour = presence), notch = TRUE) +
    facet_wrap(~location) +
    scale_y_continuous(trans = 'log10') + theme_minimal()
  return(out)
}

box.fun(chlor_a)

# compare for events 2, 4 and 5 only (i.e. the top events):
# 1. chl-a conc during event with chl-a conc prior to event
# 2. chl-a conc during event with chl-a conc post event
#
# compare for events of increasingly longer durations:
# 1. chl-a conc during event with chl-a conc prior to event
# 2. chl-a conc during event with chl-a conc post event
#
# since both cyclonic and anticyclonic eddies are being formed from the WBCs, compare the chl-a content of cyclonic vs. anticyclonic eddies:
# 1. chl-a conc during event with chl-a conc prior to event
# 2. chl-a conc during event with chl-a conc post event

# some plots of eddies
eddies <- fread("/Users/ajsmit/Downloads/tracks_AVISO_DT2014_daily_web.csv")

eddies2 <- eddies %>%
  dplyr::filter(lat >= bbox$AC[1] & lat <= bbox$AC[2]) %>%
  dplyr::filter(lon >= bbox$AC[3] + 360 & lon <= bbox$AC[4] + 360) %>%
  dplyr::mutate(lon = lon - 360,
                row_ID = seq(1:n())) %>%
  dplyr::group_by(track) %>%
  dplyr::mutate(dur = max(n)) %>%
  dplyr::ungroup()

# find only eddies that originate in the Agulhas retroflection region:
# 15 to 25°E and -45 to -35°S

track_no <- eddies2 %>%
  dplyr::filter(n == 1) %>%
  dplyr::filter(lon >= 15 & lon <= 25) %>%
  dplyr::filter(lat >= -45 & lat <= -35) %>%
  dplyr::select(track)

eddies3 <- dplyr::right_join(eddies2, track_no)

library(OneR)
eddies3$bin <- bin(eddies3$dur, method = "content")

eddies3 %>%
  # dplyr::filter(dur <= quantile(dur, 0.25)) %>%
  ggplot(aes(x = lon, y = lat)) +
  geom_path(aes(colour = cyc, group = track)) +
  stat_contour(data = bathy, aes(x = lon, y = lat, z = z, alpha = ..level..),
               col = "magenta", size = 0.15, breaks = c(-500, -1000, -2000)) +
  scale_alpha_continuous(range = c(0.5, 1.0)) +
  guides(alpha = "none") +
  facet_wrap(~ bin, nrow = 3) +
  theme_map() +
  AC.layers



