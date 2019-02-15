# 7.make_plot_DailyMonthlySubregions.R

library(RmarineHeatWaves)
library(data.table)
library(plyr)
library(tidyverse)
library(ggpubr)
library(lubridate)
library(mgcv)
library(broom)
library(doMC); doMC::registerDoMC(cores = 4)


# Functions: velocities, EKE, etc. ----------------------------------------

source("setup/functions.R")


# Set-up ------------------------------------------------------------------

# recode factors and create mean of replicate 'in' and 'out' regions within each current, etc.
subs.in <- as.tibble(fread("/Volumes/Benguela/spatial/processed/OISSTv2/WBC/subregions/all-sub-avhrr-only-v2-19810904-20180928.csv"))
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
  ungroup()

# daily
# make whole
subs.daily <- as.tibble(ddply(subs.raw, .(region, location), make_whole, .parallel = TRUE))


# Read in the region-specific components ----------------------------------

source("setup/plot.layers.R")


# The bathy data ----------------------------------------------------------

source("setup/regionDefinition.R")
bathyDir <- "/Users/ajsmit/spatial/GEBCO_2014_Grid/csv"
bathy_AC <- fread(paste0(bathyDir, "/AC_bathy.csv"))
bathy_BC <- fread(paste0(bathyDir, "/BC_bathy.csv"))
bathy_EAC <- fread(paste0(bathyDir, "/EAC_bathy.csv"))
bathy_GS <- fread(paste0(bathyDir, "/GS_bathy.csv"))
bathy_KC <- fread(paste0(bathyDir, "/KC_bathy.csv"))


# The OISST daily data ----------------------------------------------------

OISSTDir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
OISST_GS <- as.tibble(fread(paste0(OISSTDir, "/GS-avhrr-only-v2.19810901-20180930.csv")))
colnames(OISST_GS) <- c("lon", "lat", "temp", "time")

# calc SST anomaly and climatology
OISST_GS <- OISST_GS %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::mutate(anom = temp - mean(temp, na.rm = TRUE)) %>%
  dplyr::mutate(clim = mean(temp, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(time == "2014-03-24")
# rm(OISST_GS)

# extract data along 295°E (30°N to 43°N)
OISST_GS.ts <- OISST_GS %>%
  dplyr::filter(lon == 295.125) %>%
  dplyr::filter(lat >= 30 & lat <= 43)


# Events and climatologies ------------------------------------------------

# create climatology and create events for each of the regions and 'in' and 'out' groups
# events
makeEvents <- function(sst) {
  out <- RmarineHeatWaves::detect(sst[, 3:5], min_duration = 5, max_gap = 2, cold_spells = FALSE,
                climatology_start = "1983-01-01", climatology_end = "2012-12-31")
  event <- out$event
  return(event)
  }

events <- as.tibble(ddply(subs.daily, .(region, location), makeEvents, .parallel = TRUE, .progress = "text"))

# climatology
makeClimatology <- function(sst) {
  out <- RmarineHeatWaves::detect(sst[, 3:5], min_duration = 5, max_gap = 2, cold_spells = FALSE,
                                  climatology_start = "1983-01-01", climatology_end = "2012-12-31")
  clim <- out$clim[, 1:5]
  return(clim)
}

climatology <- as.tibble(ddply(subs.daily, .(region, location), makeClimatology, .parallel = TRUE, .progress = "text"))


# Agulhas Current ---------------------------------------------------------

# largest Agulhas Current events"
events %>%
  filter(region == "AC" & location == "event") %>%
  dplyr::arrange(desc(int_mean))

gv.plot <- function(climatology, region, location) {
  dat <- climatology
  reg <- region
  loc <- location
  fig <- dat %>%
    filter(region == reg & location == loc) %>%
    filter(t >= "2010-01-01" & t <= "2017-12-31") %>%
    ggplot(aes(x = t, y = temp)) +
    geom_hline(aes(yintercept = mean(seas_clim_year)), colour = "pink", size = 0.8) +
    geom_flame(aes(y2 = thresh_clim_year)) +
    geom_line(aes(y = temp), colour = "black") +
    geom_line(aes(y = seas_clim_year), colour = "green4") +
    geom_line(aes(y = thresh_clim_year), colour = "red4") +
    xlab("Date") + ylab(expression(paste("Temperature [", degree, "C]"))) +
    ggtitle(paste(reg, loc, sep = " - "))
  return(fig)
}

plt1 <- gv.plot(climatology, "AC", "reference")
plt2 <- gv.plot(climatology, "AC", "current")
plt3 <- gv.plot(climatology, "AC", "event")

ggarrange(plt1,
          plt2,
          plt3,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/AC_OISST_flame_2010-2018.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# Brazil Current ----------------------------------------------------------

# largest Brazil Current events"
events %>%
  filter(region == "BC" & location == "event") %>%
  dplyr::arrange(desc(int_mean))

plt4 <- gv.plot(climatology, "BC", "reference")
plt5 <- gv.plot(climatology, "BC", "current")
plt6 <- gv.plot(climatology, "BC", "event")

ggarrange(plt4,
          plt5,
          plt6,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/BC_OISST_flame_2010-2018.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# East Australian Current -------------------------------------------------

# largest East Australian Current events"
events %>%
  filter(region == "EAC" & location == "event") %>%
  dplyr::arrange(desc(int_mean))

plt7 <- gv.plot(climatology, "EAC", "reference")
plt8 <- gv.plot(climatology, "EAC", "current")
plt9 <- gv.plot(climatology, "EAC", "event")

ggarrange(plt7,
          plt8,
          plt9,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/EAC_OISST_flame_2010-2018.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# Gulf Stream -------------------------------------------------------------

# largest Gulf Stream events"
events %>%
  filter(region == "GS" & location == "event") %>%
  dplyr::arrange(desc(int_mean))

plt10 <- gv.plot(climatology, "GS", "reference")
plt11 <- gv.plot(climatology, "GS", "current")
plt12 <- gv.plot(climatology, "GS", "event")

ggarrange(plt10,
          plt11,
          plt12,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/GS_OISST_flame_2010-2018.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)


# Kuroshio Current --------------------------------------------------------

# largest Kuroshio Current events"
events %>%
  filter(region == "KC" & location == "event") %>%
  dplyr::arrange(desc(int_mean))

plt13 <- gv.plot(climatology, "KC", "reference")
plt14 <- gv.plot(climatology, "KC", "current")
plt15 <- gv.plot(climatology, "KC", "event")

ggarrange(plt13,
          plt14,
          plt15,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/KC_OISST_flame_2010-2018.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

# create anomaly time series and THEN detect

subs.anom <- subs.daily %>%
  group_by(region, location) %>%
  mutate(temp = temp - mean(temp, na.rm = TRUE)) %>%
  ungroup()

events.anom <- as.tibble(ddply(subs.anom, .(region, location), makeEvents, .parallel = TRUE, .progress = "text"))
clim.anom <- as.tibble(ddply(subs.anom, .(region, location), makeClimatology, .parallel = TRUE, .progress = "text"))

plt1.a <- gv.plot(clim.anom, "AC", "reference")
plt2.a <- gv.plot(clim.anom, "AC", "current")
plt3.a <- gv.plot(clim.anom, "AC", "event")

ggarrange(plt1.a,
          plt2.a,
          plt3.a,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/AC_OISST_flame_2010-2017-anom.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

plt4.a <- gv.plot(clim.anom, "BC", "reference")
plt5.a <- gv.plot(clim.anom, "BC", "current")
plt6.a <- gv.plot(clim.anom, "BC", "event")

ggarrange(plt4.a,
          plt5.a,
          plt6.a,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/BC_OISST_flame_2010-2017-anom.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

plt7.a <- gv.plot(clim.anom, "EAC", "reference")
plt8.a <- gv.plot(clim.anom, "EAC", "current")
plt9.a <- gv.plot(clim.anom, "EAC", "event")

ggarrange(plt7.a,
          plt8.a,
          plt9.a,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/EAC_OISST_flame_2010-2017-anom.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

plt10.a <- gv.plot(clim.anom, "GS", "reference")
plt11.a <- gv.plot(clim.anom, "GS", "current")
plt12.a <- gv.plot(clim.anom, "GS", "event")

ggarrange(plt10.a,
          plt11.a,
          plt12.a,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/GS_OISST_flame_2010-2017-anom.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)

plt13.a <- gv.plot(clim.anom, "KC", "reference")
plt14.a <- gv.plot(clim.anom, "KC", "current")
plt15.a <- gv.plot(clim.anom, "KC", "event")

ggarrange(plt13.a,
          plt14.a,
          plt15.a,
          ncol = 1, nrow = 3)

ggplot2::ggsave("figures/KC_OISST_flame_2010-2017-anom.pdf", width = 7.0 * (1/3), height = 5 * (1/3), scale = 3.5)



# Generalised additive models
# 1. omit GAM on 'raw' data
# 2. GAM on de-seasoned data (raw minus climatology)

# prepare data frame
climatology.1 <- climatology %>%
  mutate(de.seas = temp - seas_clim_year) %>%
  group_by(region, location) %>%
  mutate(anom = temp - mean(temp, na.rm = TRUE)) %>%
  nest()

# anomaly.gam <- gam(temp ~ s(c(1:length(temp))), data = anomaly)

# function for linear model
lm.fun <- function(df) {
  lm(de.seas ~ c(1:length(de.seas)), data = df)
}

# function for GAM
gam.fun <- function(df) {
  gam(de.seas ~ s(c(1:length(de.seas))), data = df)
}

# fit models and add into new columns
mods <- climatology.1 %>%
  mutate(mod.lm = map(data, lm.fun),
         mod.gam = map(data, gam.fun))

daily.mods <- mods %>%
  mutate(
    glance_lm = mod.lm %>% map(glance), # model summary: r-squared...
    rsq = glance_lm %>% map_dbl("r.squared"),  # extract r-squared
    tidy_lm = mod.lm %>% map(tidy), # model estimate: coeff...
    augment_lm = mod.lm %>% map(augment), # observation stats: reslocation,hat...
    res_lm = augment_lm %>% map(".resid"), # extract reslocation
    glance_gam = mod.gam %>% map(broom::glance), # GAM model summary
    AIC_gam = glance_gam %>% map_dbl("AIC"), # extract AIC
    augment_gam = mod.gam %>% map(augment),
    predict_gam = augment_gam %>% map(".fitted"),
    res_gam = augment_gam %>% map(".resid")
    )

daily.mods %>%
  unnest(data, predict_gam) %>%
  ggplot(aes(x = t, y = de.seas, group = location)) +
  geom_line(aes(colour = location), size = 0.2, alpha = 0.8) +
  geom_line(aes(y = predict_gam, colour = location), alpha = 1.0, size = 0.75) +
  facet_grid(region ~ .) +
  scale_colour_manual(values = c("green3", "blue3", "red4")) +
  labs(x = "Date", y = "Temperature anomaly (°C)",
       title = "Generalised Additive Model -- daily data")
ggplot2::ggsave("figures/Fig13_dailyGam.pdf", width = 8.25, height = 9, scale = 1.1)
# although the current and the event regions generally display the same trend, they differ
# i.t.o. their variability--the event regions show larger deviations from the trend line


# # calculate the monthly seasonal trend to remove from the monthly data below
# seas_anom <- daily.mods %>%
#   unnest(data) %>%
#   dplyr::mutate(date = floor_date(t, unit = "month")) %>%
#   dplyr::group_by(region, location, t) %>%
#   dplyr::summarise(seas_clim_year = mean(seas_clim_year, na.rm = TRUE)) %>%
#   dplyr::group_by(region, location) %>%
#   dplyr::mutate(seas_anom = seas_clim_year - mean(seas_clim_year, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::select(seas_anom) %>%
#   as.tibble()
#
# # monthly (subs.daily to subs.monthly)
# subs.monthly <- subs.daily %>%
#   dplyr::mutate(date = floor_date(t, unit = "month"),
#                 year = year(date)) %>%
#   dplyr::group_by(region, location, date, year) %>%
#   dplyr::summarise(temp = mean(temp, na.rm = TRUE)) %>%
#   dplyr::group_by(region, location) %>%
#   dplyr::mutate(anom = temp - mean(temp, na.rm = TRUE)) %>%
#   dplyr::mutate(num = seq(1:length(temp))) %>%
#   dplyr::ungroup() %>%
#   as.tibble()
#
# subs.monthly <- as.tibble(cbind(subs.monthly, seas_anom))
# subs.monthly <- subs.monthly %>%
#   # dplyr::mutate(de.seas = anom - seas_anom) %>%
#   dplyr::group_by(region, location) %>%
#   nest()
#
# # a gls function
# gls.fun <- function(df) {
#   gls(temp ~ num, correlation = corARMA(form = ~ 1 | year, p = 2),
#       method = "REML", data = df, na.action = na.exclude)
# }
#
# # fit a GAM
# monthly.mod1 <- subs.monthly %>%
#   mutate(mod.gam = map(data, gam.fun))
#
# monthly.mod1 <- monthly.mod1 %>%
#   mutate(
#     glance_gam = mod.gam %>% map(broom::glance), # GAM model summary
#     AIC_gam = glance_gam %>% map_dbl("AIC"), # extract AIC
#     augment_gam = mod.gam %>% map(augment),
#     predict_gam = augment_gam %>% map(".fitted"),
#     res_gam = augment_gam %>% map(".resid")
#   )
#
# monthly.mod1 %>%
#   unnest(data, predict_gam) %>%
#   ggplot(aes(x = date, y = de.seas)) +
#   geom_line(aes(colour = location), size = 0.4, alpha = 1) +
#   geom_line(aes(y = predict_gam, colour = location), alpha = 1.0, size = 0.75) +
#   facet_grid(region ~ .) +
#   scale_colour_manual(values = c("red4", "blue3")) +
#   labs(x = "Date", y = "Temperature anomaly (°C)",
#        title = "Generalised Additive Model -- montly means")
# ggplot2::ggsave("figures/Fig13_monthlyGAM.pdf", width = 8.25, height = 9, scale = 1.1)
#
# # try a gamm
# # a gamm function
# gamm.fun <- function(df) {
#   out <- gamm(de.seas ~ s(num), data = df, correlation = corARMA(form = ~ 1, p = 2))
#   return(out$gam)
# }
#
# monthly.mod2 <- subs.monthly %>%
#   mutate(mod.gamm = map(data, gamm.fun))
#
# # a custom augment function because broom doesn't have one
# augment.gamm <- function(x, data = stats::model.frame(x), ...) {
#   augment_columns(x, data)
# }
#
# # massage the thing into something useful
# monthly.mapped <- map(monthly.mod2$mod.gamm, augment.gamm)
# names(monthly.mapped) <- with(subs.monthly, str_c(region, location, sep = "_"))
# monthly.mod2 <- ldply(monthly.mapped, data.frame, .progress = "text")
# monthly.mod2 <- as.tibble(select(monthly.mod2, -c(Xr, X)))
# monthly.mod2 <- separate(monthly.mod2, .id, c("region", "location"), sep = "_")
#
# monthly.mod2 %>%
#   ggplot(aes(x = num, y = de.seas)) +
#   geom_line(aes(colour = location), size = 0.4, alpha = 1) +
#   geom_line(aes(y = .fitted, colour = location), alpha = 1.0, size = 0.75) +
#   facet_grid(region ~ .) +
#   scale_colour_manual(values = c("red4", "blue3")) +
#   labs(x = "Date", y = "Temperature anomaly (°C)",
#        title = "Generalised Additive Mixed Model (AR1) -- montly means")
# ggplot2::ggsave("figures/Fig13_monthlyGAMM.pdf", width = 8.25, height = 9, scale = 1.1)
