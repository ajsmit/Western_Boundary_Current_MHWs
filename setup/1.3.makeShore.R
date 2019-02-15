# 1.3.makeShore.R

library(data.table)
library(ggplot2)

# setup
source("setup/regionDefinition.R")
source("setup/functions.R")
shoreDir <- "/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/shorelines"
gshhsDir <- "/Users/ajsmit/spatial/gshhg-bin-2.3.7"

# make the shores
makeShore(shoreDir, gshhsDir, "AC")
makeShore(shoreDir, gshhsDir, "BC")
makeShore(shoreDir, gshhsDir, "BC")
makeShore(shoreDir, gshhsDir, "EAC")
makeShore(shoreDir, gshhsDir, "KC")
# makeShore(shoreDir, gshhsDir, "GS") # doesn't work due to self-intersection; do differently...
lats <- c(bbox["latmin", "GS"], bbox["latmax", "GS"])
lons <- c(bbox["lonmin", "GS"], bbox["lonmax", "GS"])
shore <- Rgshhs(paste0(gshhsDir, "/gshhs_f.b"), ## run Rgshhs directly; returns list
                xlim = lons, ylim = lats, level = 1, no.clip = FALSE, checkPolygons = TRUE)
shore <- fortify(shore$SP)
save(shore, file = paste0(shoreDir, "/", "GS", "-shore.Rdata"))
remove(shore)

# plot the shorelines
baseMap(shoreDir, "AC-shore.Rdata")
baseMap(shoreDir, "BC-shore.Rdata") # ???
load(paste0(shoreDir, "/", "BC-shore.Rdata")) # the plotting of the BC requires a different lon specification
coords <- c(-45, -5, -60, -25)
ggplot(shore, aes(x = lon, y = lat)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               fill = "#929292", colour = "black", size = 0.1, show.legend = FALSE) +
  xlab("Longitude (°E)") +
  ylab("Latitude (°S)") +
  coord_fixed(ratio = 1,
              xlim = c(coords[3], coords[4]),
              ylim = c(coords[1], coords[2]), expand = FALSE) +
  theme_bw()
ggsave(paste0(shoreDir, "/", "BC", "-shore.pdf"), width = 8, height = 8)
baseMap(shoreDir, "EAC-shore.Rdata")
baseMap(shoreDir, "KC-shore.Rdata")
baseMap(shoreDir, "GS-shore.Rdata")
