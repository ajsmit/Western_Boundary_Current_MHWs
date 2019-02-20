# 1.regional_setup.R


# Plotting functions ------------------------------------------------------

current_uv_scalar <- 2.25
vel.plot <- function(sl.dat, bathy.dat, plot.parameters) {
  vec <- oceFun2(sl.dat)
  fig <- ggplot(sl.dat, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = eke), interpolate = FALSE) +
    scale_fill_gradientn(colours = rev(rainbow(7, end = 4/6)),
                         space = "Lab", limits = c(0, 1.45)) +
    stat_contour(data = bathy.dat, aes(x = lon, y = lat, z = z, alpha = ..level..),
                 col = "black", size = 0.15, binwidth = 500) +
    scale_alpha_continuous(range = c(0.2, 0.9)) +
    geom_segment(data = vec,
                 aes(xend = lon + u * current_uv_scalar,
                     yend = lat + v * current_uv_scalar,
                     size = velocity / 2),
                 colour = "red",
                 arrow = arrow(angle = 17.5, length = unit(0.1, "cm"), type = "open"),
                 linejoin = "mitre") +
    scale_size(range = c(0.025, 0.25)) +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(EKE~(m^2~s^{-2})),
                                  frame.colour = "black",
                                  frame.linewidth = 0.5,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(100, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'top',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    theme_map() +
    plot.parameters
  return(fig)
}

sla.plot <- function(sl.dat, bathy.dat, plot.parameters) {
  vec <- oceFun2(sl.dat)
  fig <- ggplot(sl.dat, aes(x = lon, y = lat)) +
    geom_raster(aes(fill = sla), interpolate = FALSE) +
    scale_fill_gradientn(colours = col1, limits = c(-1.2, 1.2)) +
    stat_contour(data = bathy.dat, aes(x = lon, y = lat, z = z, alpha = ..level..),
                 col = "black", size = 0.15, binwidth = 500) +
    scale_alpha_continuous(range = c(0.2, 0.9)) +
    geom_segment(data = vec,
                 aes(xend = lon + u * current_uv_scalar,
                     yend = lat + v * current_uv_scalar,
                     size = velocity / 2),
                 colour = "red",
                 arrow = arrow(angle = 17.5, length = unit(0.1, "cm"), type = "open"),
                 linejoin = "mitre") +
    scale_size(range = c(0.025, 0.25)) +
    guides(alpha = "none",
           fill = guide_colourbar(title = expression(SLA~(m)),
                                  frame.colour = "black",
                                  frame.linewidth = 0.5,
                                  ticks.colour = "black",
                                  barheight = unit(2, units = "mm"),
                                  barwidth = unit(100, units = "mm"),
                                  draw.ulim = F,
                                  title.position = 'top',
                                  title.hjust = 0.5,
                                  label.hjust = 0.5),
           size = "none") +
    theme_map() +
    plot.parameters
  return(fig)
}

# plot the SST anomaly on the day of the max mean intensity (Peak: 2014-03-24)
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

flame.plot <- function(climatology, region, location) {
  dat <- climatology
  reg <- region
  loc <- location
  fig <- dat %>%
    filter(region == reg & location == loc) %>%
    filter(t >= "2010-01-01" & t <= "2018-09-30") %>%
    ggplot(aes(x = t, y = temp)) +
    geom_hline(aes(yintercept = mean(seas)), colour = "pink", size = 0.8) +
    geom_flame(aes(y2 = thresh)) +
    geom_line(aes(y = temp), colour = "black") +
    geom_line(aes(y = seas), colour = "green4") +
    geom_line(aes(y = thresh), colour = "red4") +
    xlab("Date") + ylab(expression(paste("Temperature [", degree, "C]"))) +
    ggtitle(paste(reg, loc, sep = " - ")) + theme_minimal()
  return(fig)
}
