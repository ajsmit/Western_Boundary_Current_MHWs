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


# Subet AVISO and calculate KE --------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

## NB: This is too beefy to run in a for loop
# Run sequentially
# for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   AVISO_KE_save(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
#   Sys.sleep(10) # Let the server catch its breath
# } #~10 -- 15 minutes per region
# AVISO_KE_save(colnames(bbox)[1])
# AVISO_KE_save(colnames(bbox)[2])
# AVISO_KE_save(colnames(bbox)[3])
AVISO_KE_save(colnames(bbox)[4])
# AVISO_KE_save(colnames(bbox)[5])


# EKE percentile masks ----------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Run sequentially
# for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   mke_masks(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# }


# Subsetting MHW results --------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Subset MHW data
# for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   MHW_sub_save(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# }


# Correlate pixels --------------------------------------------------------
# In this final step we take only the MKE and MHW pixels that are in the 
# 75th percentile mask but not in the 90th percentile mask.
# We then take the days the pixels in the 75th percentile areas
# are themselves in the 90th percentile for MKE.
# A count of how often this occurs when MHWs are occurring is made to find
# what the proportion of this occurence is.
# The MKE and mean intensities on days when MKE is in the 90th perc. 
# are then correlated.
# This result will help to illustrate the potential relationship between
# meanders and MHWs.


# Visualise results -------------------------------------------------------


