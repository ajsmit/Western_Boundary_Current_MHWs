# 8. ICOADS.R

library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(doMC); doMC::registerDoMC(cores = 4)
library(viridis)

# the GLS-AR2 -------------------------------------------------------------
lm_fun <- function(df) {
  # require(nlme)
  out <- tryCatch( {
    model <- lm(
      temp ~ num,
      data = df, na.action = na.exclude
    )
    stats <-
      data.frame(
        DT_model = round(as.numeric(coef(model)[2]) * 120, 3),
        se_trend = round(summary(model)[["tTable"]][[2, 2]], 10),
        sd_initial = round(sd(df$temp, na.rm = TRUE), 2),
        sd_residual = round(sd(model$residuals), 2),
        r = NA,
        R2 = NA,
        p_trend = summary(model)[["tTable"]][[2, 4]],
        p_seas = NA,
        length = length(df$temp),
        model = "lm"
      )
    return(stats)
  },
  error = function(cond) {
    stats <- data.frame(
      DT_model = NA,
      se_trend = NA,
      sd_initial = round(sd(df$temp, na.rm = TRUE), 2),
      sd_residual = NA,
      r = NA,
      R2 = NA,
      p_trend = NA,
      p_seas = NA,
      length = length(df$temp),
      model = "lm"
    )
    return(stats)
  })
  rm(model)
  return(out)
}

# setup ICOADS
ICOADS.AC.in <- fread("/Volumes/Benguela/OceanData/ICOADS/ICOADS_Release_3_Monthly_Summaries/subsetted/AC/MSG1.S.std.196001.201710_1",
                      skip = 1)
ICOADS.AC.in <- as_tibble(ICOADS.AC.in)

source("setup/regionDefinition.R")
source("setup/theme.R")
theme_set(bw_update) # Set theme

# Make world coastline
library(PBSmapping)
gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"
lats <- c(-90, 90)
lons <- c(0, 360)
shore <- importGSHHS(paste0(gshhsDir, "/gshhs_i.b"), xlim = lons, ylim = lats, maxLevel = 1, useWest = FALSE)

# AC = c(-45, -20, 6.25, 45)

library(zoo)
ICOADS.AC <- ICOADS.AC.in %>%
  dplyr::mutate(date = as.Date(as.yearmon(paste(ICOADS.AC.in$YEAR, ICOADS.AC.in$MON), "%Y %m"))) %>%
  dplyr::select(date, YEAR, MON, BLO, BLA, M) %>%
  dplyr::filter(date > "1981-08-01") %>%
  dplyr::rename(year = YEAR,
                month = MON,
                lon = BLO,
                lat = BLA,
                temp = M) %>%
  dplyr::filter(lon >= 6.25 & lon <= 45) %>%
  dplyr::filter(lat >= -45 & lat <= -20) %>%
  dplyr::group_by(date, lon, lat) %>%
  dplyr::summarize(temp = mean(temp, na.rm = TRUE)) %>%
  dplyr::group_by(lon, lat) %>%
  dplyr::mutate(num = as.numeric(date) / 1000) %>%
  dplyr::ungroup()

system.time(AC_gls <- dlply(ICOADS.AC, .(lon, lat), .progress = "text",
                             .parallel = FALSE, gls_fun))
AC_df <- as_tibble(ldply(AC_gls, data.frame, .progress = "text"))

system.time(AC_lm <- dlply(ICOADS.AC, .(lon, lat), .progress = "text",
                            .parallel = FALSE, gls_fun))
AC.lm_df <- as_tibble(ldply(AC_lm, data.frame, .progress = "text"))

ICOADS.AC %>%
  ggplot(aes(x = lon, y = lat, fill = temp)) +
  # stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE) +
  stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE, h = c(1, 1)) +
  scale_fill_viridis(name = "Total:\ncount") +
  geom_vline(xintercept = seq(10, 45, 5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = seq(-45, -20, 5), linetype = "dashed", size = 0.2) +
  geom_polygon(data = shore, aes(x = X, y = Y, group = PID),
               colour = "black", fill = "gray98", size = 0.2) +
  coord_fixed(ratio = 1, xlim = c(-1.875, 53.125), ylim = c(-52.5, -12.5), expand = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = scales::unit_format("°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0), labels = scales::unit_format("°N", sep = "")) +
  labs(title = "A", x = NULL, y = NULL)
