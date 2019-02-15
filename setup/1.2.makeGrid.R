# 1.2.makeGrid.R

library(data.table)

# setup
source("setup/regionDefinition.R")
source("setup/functions.R")
inDir <- "/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/daily"
outDir <- "/Volumes/Benguela/spatial/OISSTv2/daily/csv/WBC/grids"

# make the grids
setupGrid(inDir, outDir, "AC-avhrr-only-v2.19810901-20171019.csv", byYear = TRUE)
setupGrid(inDir, outDir, "EAC-avhrr-only-v2.19810901-20171019.csv", byYear = TRUE)
setupGrid(inDir, outDir, "BC-avhrr-only-v2.19810901-20171019.csv", byYear = TRUE)
setupGrid(inDir, outDir, "KC-avhrr-only-v2.19810901-20171019.csv", byYear = TRUE)
setupGrid(inDir, outDir, "GS-avhrr-only-v2.19810901-20171019.csv", byYear = TRUE)
