# 6.make_daily_subsets.R

library(data.table)
library(tidyverse)
library(lubridate)

source("setup/functions.R")


# Load subregion definition -----------------------------------------------

subRegions <- as.tibble(fread("setup/subRegions.csv"))


# A function --------------------------------------------------------------

# region <- "KC"
# is.west <- TRUE

# a function to do the subsetting
# set is.west = TRUE when the lons for regions west of the date line are
# formatted as degrees west (i.e. -180 to 0); else, treated as 0 - 360
subExtract <- function(data, region, is.west = FALSE) {
  sub.def <- subRegions[subRegions$region == region, ]
  if (is.west) {
    sub.def$lonmin <- sub.def$lonmin - 360
    sub.def$lonmiax<- sub.def$lonmax - 360 }
  else {
    sub.def <- sub.def
  }
  fun1 <- function(subregion) {
    sub <- data %>%
      filter(lon >= as.numeric(sub.def[subregion, "lonmin"]) & lon <= as.numeric(sub.def[subregion, "lonmax"])) %>%
      filter(lat >= as.numeric(sub.def[subregion, "latmin"]) & lat <= as.numeric(sub.def[subregion, "latmax"])) %>%
      mutate(region = region) %>%
      mutate(ID = paste0("sub", subregion))
    return(sub)
  }
  subs <- bind_rows(fun1(1), fun1(2), fun1(3), fun1(4), fun1(5), fun1(6), fun1(7), fun1(8), fun1(9))
  return(subs)
}


# Subset OISST ------------------------------------------------------------

csv.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily"
sub.dir <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/subregions"

AC.sst <- as.tibble(fread(paste0(csv.dir, "/", "AC-avhrr-only-v2.19810901-20180930.csv"),
                          colClasses = list("numeric" = 1:3, "myDate" = 4),
                          col.names = c("lon", "lat", "temp", "t")))
AC.sub <- subExtract(AC.sst, "AC")
rm(AC.sst)

BC.sst <- as.tibble(fread(paste0(csv.dir, "/", "BC-avhrr-only-v2.19810901-20180930.csv"),
                          colClasses = list("numeric" = 1:3, "myDate" = 4),
                          col.names = c("lon", "lat", "temp", "t")))
BC.sub <- subExtract(BC.sst, "BC")
rm(BC.sst)

EAC.sst <- as.tibble(fread(paste0(csv.dir, "/", "EAC-avhrr-only-v2.19810901-20180930.csv"),
                           colClasses = list("numeric" = 1:3, "myDate" = 4),
                           col.names = c("lon", "lat", "temp", "t")))
EAC.sub <- subExtract(EAC.sst, "EAC")
rm(EAC.sst)

KC.sst <- as.tibble(fread(paste0(csv.dir, "/", "KC-avhrr-only-v2.19810901-20180930.csv"),
                          colClasses = list("numeric" = 1:3, "myDate" = 4),
                          col.names = c("lon", "lat", "temp", "t")))
KC.sub <- subExtract(KC.sst, "KC")
rm(KC.sst)

GS.sst <- as.tibble(fread(paste0(csv.dir, "/", "GS-avhrr-only-v2.19810901-20180930.csv"),
                          colClasses = list("numeric" = 1:3, "myDate" = 4),
                          col.names = c("lon", "lat", "temp", "t")))
GS.sub <- subExtract(GS.sst, "GS")
rm(GS.sst)

# combine and save them
strt <- gsub("-", "", AC.sub[1, 4])
end <- gsub("-", "", AC.sub[nrow(AC.sub), 4])

sub.sst <- bind_rows(AC.sub, BC.sub, EAC.sub, KC.sub, GS.sub)
fwrite(sub.sst,
       file = paste0(sub.dir, "/", "all-sub-avhrr-only-v2-", strt, "-", end, ".csv"))
rm(AC.sub); rm(BC.sub); rm(EAC.sub); rm(KC.sub); rm(GS.sub)


# Subset Aviso+ -----------------------------------------------------------

csvDir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/daily"
sub.dir <- "/Volumes/Benguela/spatial/processed/Aviso+/WBC/subregions"

AC.file <- "AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
AC.sl <- as.tibble(fread(paste0(csvDir, "/", AC.file),
                         colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                         col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
AC.sub <- subExtract(AC.sl, "AC", is.west = FALSE)
rm(AC.sl); rm(AC.file)

BC.file <- "BC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
BC.sl <- as.tibble(fread(paste0(csvDir, "/", BC.file),
                         colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                         col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
BC.sub <- subExtract(BC.sl, "BC", is.west = FALSE)
rm(BC.sl); rm(BC.file)

# EAC = c(-42.5, -15, 145, 160), # East Australian Current
EAC.file <- "EAC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
EAC.sl <- as.tibble(fread(paste0(csvDir, "/", EAC.file),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                          col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
EAC.sub <- subExtract(EAC.sl, "EAC", is.west = FALSE)
EAC.sl <- EAC.sl %>% dplyr::filter(lat <= -15)
rm(EAC.sl); rm(EAC.file)

KC.file <- "KC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
KC.sl <- as.tibble(fread(paste0(csvDir, "/", KC.file),
                         colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                         col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
KC.sub <- subExtract(KC.sl, "KC", is.west = FALSE)
rm(KC.sl); rm(KC.file)

GS.file <- "GS-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
GS.sl <- as.tibble(fread(paste0(csvDir, "/", GS.file),
                         colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4:7),
                         col.names = c("lon", "lat", "time", "u", "v", "sla", "adt")))
GS.sub <- subExtract(GS.sl, "GS", is.west = FALSE)
rm(GS.sl); rm(GS.file)

# and combine and save them
strt <- gsub("-", "", AC.sub[1, 3])
end <- gsub("-", "", AC.sub[nrow(AC.sub), 3])

sub.sl <- bind_rows(AC.sub, BC.sub, EAC.sub, KC.sub, GS.sub)
rm(AC.sub); rm(BC.sub); rm(EAC.sub); rm(KC.sub); rm(GS.sub)
fwrite(sub.sl,
       file = paste0(sub.dir, "/", "all-sub-Aviso+-", strt, "-", end, ".csv"))

sub.sl.grp <- sub.sl %>%
  dplyr::group_by(region, ID, time) %>%
  dplyr::summarise(ugos = round(mean(u, na.rm = TRUE), 3),
                   vgos = round(mean(v, na.rm = TRUE), 3),
                   sla = round(mean(sla, na.rm = TRUE), 3),
                   adt = round(mean(adt, na.rm = TRUE), 3)) %>%
  ungroup()

fwrite(sub.sl.grp,
       file = paste0(sub.dir, "/", "all-sub-summarised-Aviso+-", strt, "-", end, ".csv"))


# Subset MODIS Aqua -------------------------------------------------------

csv.dir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/regridded"
sub.dir <- "/Volumes/Benguela/spatial/processed/MODIS-Aqua_Chlorophyll_Concentration_OCI_Algorithm/WBC/subregions"

AC.chl <- as.tibble(fread(paste0(csv.dir, "/", "AC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4)))
AC.sub <- subExtract(AC.chl, "AC")
rm(AC.chl)

BC.chl <- as.tibble(fread(paste0(csv.dir, "/", "BC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4)))
BC.sub <- subExtract(BC.chl, "BC")
rm(BC.chl)

EAC.chl <- as.tibble(fread(paste0(csv.dir, "/", "EAC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4)))
EAC.sub <- subExtract(EAC.chl, "EAC")
rm(EAC.chl)

GS.chl <- as.tibble(fread(paste0(csv.dir, "/", "GS.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4)))
GS.sub <- subExtract(GS.chl, "GS")
rm(GS.chl)

KC.chl <- as.tibble(fread(paste0(csv.dir, "/", "KC.MODIS.Aqua.L3m_8D_CHL_chlor_a_regrid_2002-07-05-2018-11-18.csv"),
                          colClasses = list("numeric" = 1:2, "myDate" = 3, "numeric" = 4)))
KC.sub <- subExtract(KC.chl, "KC")
rm(KC.chl)

# and combine and save them
strt <- gsub("-", "", AC.sub[1, 3])
end <- gsub("-", "", AC.sub[nrow(AC.sub), 3])

sub.chl <- bind_rows(AC.sub, BC.sub, EAC.sub, KC.sub, GS.sub)
rm(AC.sub); rm(BC.sub); rm(EAC.sub); rm(KC.sub); rm(GS.sub)
fwrite(sub.chl,
       file = paste0(sub.dir, "/", "all-sub-MODIS-Aqua-", strt, "-", end, ".csv"))

sub.chl.grp <- sub.chl %>%
  dplyr::group_by(region, ID, time) %>%
  dplyr::summarise(chlor_a = round(mean(chlor_a, na.rm = TRUE), 3)) %>%
  ungroup()

fwrite(sub.chl.grp,
       file = paste0(sub.dir, "/", "all-sub-summarised-MODIS-Aqua-", strt, "-", end, ".csv"))


# UNUSED ------------------------------------------------------------------

# Subset MW_OI-REMSS-L4-GLOB-v4.0 -----------------------------------------

csv.dir <- "/Users/ajsmit/spatial/processed/MW_OI-REMSS-L4-GLOB-v4.0/WBC/daily"
sub.dir <- "/Users/ajsmit/spatial/processed/MW_OI-REMSS-L4-GLOB-v4.0/WBC/daily/subregions"
names <- c("lon", "lat", "temp", "t")

AC.sst <- as.tibble(fread(paste0(csv.dir, "/", "AC-MW_OI-REMSS-L4-GLOB-v4.0-19980101-20171231.csv")))
colnames(AC.sst) <- names
AC.sub <- subExtract(AC.sst, "AC")
rm(AC.sst)

fwrite(AC.sub,
       file = paste0(sub.dir, "/", "AC-sub-MW_OI-REMSS-L4-GLOB-v4.0-19980101-20171231.csv"))

AC.sub.grp <- AC.sub %>%
  dplyr::group_by(region, ID, t) %>%
  dplyr::summarise(temp = round(mean(temp, na.rm = TRUE), 3)) %>%
  dplyr::mutate(temp = temp - 273.15) %>%
  ungroup()

fwrite(AC.sub.grp,
       file = paste0(sub.dir, "/", "AC-sub-summarised-MW_OI-REMSS-L4-GLOB-v4.0-19980101-20171231.csv"))


# Subset CMC --------------------------------------------------------------

csv.dir <- "/Volumes/Challenger_Deep/spatial/processed/CMC0.2deg-CMC-L4-GLOB-v2.0/WBC/daily"
sub.dir <- "/Volumes/Challenger_Deep/spatial/processed/CMC0.2deg-CMC-L4-GLOB-v2.0/WBC/subregions"
names <- c("lon", "lat", "temp", "t")

AC.sst <- as.tibble(fread(paste0(csv.dir, "/", "AC-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv")))
colnames(AC.sst) <- names
AC.sub <- subExtract(AC.sst, "AC")
# fwrite(AC.sub,
#        file = paste0(sub.dir, "/", "AC-sub-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))
rm(AC.sst)

BC.sst <- as.tibble(fread(paste0(csv.dir, "/", "BC-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv")))
colnames(BC.sst) <- names
BC.sub <- subExtract(BC.sst, "BC", is.west = TRUE)
rm(BC.sst)

EAC.sst <- as.tibble(fread(paste0(csv.dir, "/", "EAC-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv")))
colnames(EAC.sst) <- names
EAC.sub <- subExtract(EAC.sst, "EAC")
# fwrite(EAC.sub,
#        file = paste0(sub.dir, "/", "EAC-sub-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))
rm(EAC.sst)

KC.sst <- as.tibble(fread(paste0(csv.dir, "/", "KC-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv")))
colnames(KC.sst) <- names
KC.sub <- subExtract(KC.sst, "KC")
# fwrite(KC.sub,
#        file = paste0(sub.dir, "/", "KC-sub-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))
rm(KC.sst)

GS.sst <- as.tibble(fread(paste0(csv.dir, "/", "GS-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv")))
colnames(GS.sst) <- names
GS.sub <- subExtract(GS.sst, "GS", is.west = TRUE)
# fwrite(GS.sub,
#        file = paste0(sub.dir, "/", "GS-sub-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))
rm(GS.sst)

# and combine and save them
sub.sst <- bind_rows(AC.sub, BC.sub, EAC.sub, KC.sub, GS.sub)
rm(AC.sub); rm(BC.sub); rm(EAC.sub); rm(KC.sub); rm(GS.sub)
fwrite(sub.sst,
       file = paste0(sub.dir, "/", "all-sub-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))

sub.sst.grp <- sub.sst %>%
  dplyr::mutate(temp = temp - 273.15) %>%
  dplyr::group_by(region, ID, t) %>%
  dplyr::summarise(temp = round(mean(temp, na.rm = TRUE), 3)) %>%
  ungroup()

fwrite(sub.sst.grp,
       file = paste0(sub.dir, "/", "all-sub-summarised-CMC0.2deg-GLOB-v02.0-19910901-20170317.csv"))
