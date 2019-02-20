# 4.4.Aviso+_overallMean.R

library(data.table)
library(tidyverse)


# Functions: velocities, EKE, etc. ----------------------------------------

source("setup/functions.R")

csvInDir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/daily"
csvOutDir <- "/Volumes/Benguela/spatial/processed/Aviso/WBC/mean"


# Calculate monthly Aviso+ data -------------------------------------------

## apply the function, and save the data -- do once only when new data are
## available, then load the saved data

AC.file <- "AC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
AC.sl <- as.tibble(fread(paste0(csvInDir, "/", AC.file), verbose = FALSE))
dateRng <- range(fastDate(AC.sl$V3))
dateStrt <- str_replace_all(dateRng[1], "-", "")
dateEnd <- str_replace_all(dateRng[2], "-", "")
mAC.sl <- oceFun1(AC.sl)
fwrite(mAC.sl, paste0(csvOutDir, "/mean_AC_sla-", dateStrt, "-", dateEnd, ".csv"))
rm(AC.sl); rm(mAC.sl)

BC.file <- "BC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
BC.sl <- as.tibble(fread(paste0(csvInDir, "/", BC.file), verbose = FALSE))
dateRng <- range(fastDate(BC.sl$V3))
dateStrt <- str_replace_all(dateRng[1], "-", "")
dateEnd <- str_replace_all(dateRng[2], "-", "")
mBC.sl <- oceFun1(BC.sl)
fwrite(mBC.sl, paste0(csvOutDir, "/mean_BC_sla-", dateStrt, "-", dateEnd, ".csv"))
rm(BC.sl); rm(mBC.sl)

# EAC = c(-42.5, -15, 145, 160), # East Australian Current
EAC.file <- "EAC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
EAC.sl <- as.tibble(fread(paste0(csvInDir, "/", EAC.file), verbose = FALSE))
dateRng <- range(fastDate(EAC.sl$V3))
dateStrt <- str_replace_all(dateRng[1], "-", "")
dateEnd <- str_replace_all(dateRng[2], "-", "")
mEAC.sl <- oceFun1(EAC.sl)
mEAC.sl <- mEAC.sl %>% dplyr::filter(lat <= -15)
fwrite(mEAC.sl, paste0(csvOutDir, "/mean_EAC_sla-", dateStrt, "-", dateEnd, ".csv"))
rm(EAC.sl); rm(mEAC.sl)

KC.file <- "KC-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
KC.sl <- as.tibble(fread(paste0(csvInDir, "/", KC.file), verbose = FALSE))
dateRng <- range(fastDate(KC.sl$V3))
dateStrt <- str_replace_all(dateRng[1], "-", "")
dateEnd <- str_replace_all(dateRng[2], "-", "")
mKC.sl <- oceFun1(KC.sl)
fwrite(mKC.sl, paste0(csvOutDir, "/mean_KC_sla-", dateStrt, "-", dateEnd, ".csv"))
rm(KC.sl); rm(mKC.sl)

GS.file <- "GS-dataset-duacs-rep-global-merged-allsat-phy-l4-v3.csv"
GS.sl <- as.tibble(fread(paste0(csvInDir, "/", GS.file), verbose = FALSE))
dateRng <- range(fastDate(GS.sl$V3))
dateStrt <- str_replace_all(dateRng[1], "-", "")
dateEnd <- str_replace_all(dateRng[2], "-", "")
mGS.sl <- oceFun1(GS.sl)
fwrite(mGS.sl, paste0(csvOutDir, "/mean_GS_sla-", dateStrt, "-", dateEnd, ".csv"))
rm(GS.sl); rm(mGS.sl)
