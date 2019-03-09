# The purpose of this script is to first create the AVISO and MHW databases
# used throughout this analysis because this script is to be run on the 
# tikoraluk server and the database files are too large to send via
# conventional digital means
# These files are created here by:
  # 1) Loading the bounding box definitions
  # 2) Subsetting from the MHW and AVISO databases on tikoraluk
    # 2.1) The grids for these data already match on tikoraluk 
      # so no converting is required

# With these databases recreated here, the next step is to re-calculate KE

# With all of the results ready, the next step is to use the 90th perc.
# EKE mask and the 75th perc. EKE mask to define the region within which 
# we want to compare poleward velocity (U+V) and mean intensity

# With the data masked the correlations may then be run


# Functions ---------------------------------------------------------------

source("setup/meanderFunctions.R")


# Subsetting AVISO --------------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Subset AVISO data
for(i in 1:ncol(bbox)){
  region <- colnames(bbox)[i]
  print(paste0("Began run on ",region," at ",Sys.time()))
  AVISO_sub_save(region)
  print(paste0("Finished run on ",region," at ",Sys.time()))
} #~4 -- 10 minutes per region


# Subsetting MHW results --------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Subset MHW data
# for(i in 1:ncol(bbox)){
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   MHW_sub_save(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# }


# Calculate EKE -----------------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

## NB: This is a bit too beefy
  ## It looks like R will need to be manually restarted after each run
  ## to effectively purge the memmory
# Calculate EKE
# for(i in 1:ncol(bbox)){
#   region <- colnames(bbox)[i]
  # print(paste0("Began run on ",region," at ",Sys.time()))
  # ke_region_save(region)
  # print(paste0("Finished run on ",region," at ",Sys.time()))
# }

## Individual runs
# ke_region_save(colnames(bbox)[1])
# ke_region_save(colnames(bbox)[2])
# ke_region_save(colnames(bbox)[3])
# ke_region_save(colnames(bbox)[4])
# ke_region_save(colnames(bbox)[5])
# ke_region_save(colnames(bbox)[6])


# EKE percentile masks ----------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)


# Filter 90th and 75th areas ----------------------------------------------


# Correlate pixels --------------------------------------------------------
# In this final step we take only the EKE and MHW pixels that are in the 
# 75th perc. mask but not in the 90th perc. mask
# We then take the days the pixels in the 75th to 90th perc. areas
# are themselves in the 90th perc. for EKE
# A count of how often this occurs and MHWs occur is made to find what
# the proportion of this occurence is
# The EKE and mean intensities on days when EKE is in the 90th perc. 
# are then correlated
# Thiw result will help to illustrate the potential relationship between
# meanders and MHWs

