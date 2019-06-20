library(rgdal)
library(maptools)

load("masks/AC-mask_polys.RData")
writeOGR(obj = mask.list$int90, dsn = "masks/shapedir", layer = "AC_int90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$eke90, dsn = "masks/shapedir", layer = "AC_eke90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$mke90, dsn = "masks/shapedir", layer = "AC_mke90", driver = "ESRI Shapefile")

load("masks/EAC-mask_polys.RData")
writeOGR(obj = mask.list$int90, dsn = "masks/shapedir", layer = "EAC_int90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$eke90, dsn = "masks/shapedir", layer = "EAC_eke90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$mke90, dsn = "masks/shapedir", layer = "EAC_mke90", driver = "ESRI Shapefile")

load("masks/KC-mask_polys.RData")
writeOGR(obj = mask.list$int90, dsn = "masks/shapedir", layer = "KC_int90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$eke90, dsn = "masks/shapedir", layer = "KC_eke90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$mke90, dsn = "masks/shapedir", layer = "KC_mke90", driver = "ESRI Shapefile")

# Gulf Stream and Brazil Current shapefiles must be shifted...

shift.xy <- c(-360, 0)

load("masks/GS-mask_polys.RData")
writeOGR(obj = mask.list$int90, dsn = "masks/shapedir", layer = "GS_int90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$eke90, dsn = "masks/shapedir", layer = "GS_eke90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$mke90, dsn = "masks/shapedir", layer = "GS_mke90", driver = "ESRI Shapefile")
xx <- readOGR("masks/shapedir/GS_int90.shp") # only shifting MHW intensity here...
## update the geometry with elide arguments
shifted <- elide(xx, shift = shift.xy)
## write out a new shapfile
writeOGR(shifted, "masks/shapedir/GS_int90_shifted.shp", layer = "int90", driver = "ESRI Shapefile")


load("masks/BC-mask_polys.RData")
writeOGR(obj = mask.list$int90, dsn = "masks/shapedir", layer = "BC_int90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$eke90, dsn = "masks/shapedir", layer = "BC_eke90", driver = "ESRI Shapefile")
writeOGR(obj = mask.list$mke90, dsn = "masks/shapedir", layer = "BC_mke90", driver = "ESRI Shapefile")
xx <- readOGR("masks/shapedir/BC_int90.shp") # only shifting MHW intensity here...
## update the geometry with elide arguments
shifted <- elide(xx, shift = shift.xy)
## write out a new shapfile
writeOGR(shifted, "masks/shapedir/BC_int90_shifted.shp", layer = "int90", driver = "ESRI Shapefile")
