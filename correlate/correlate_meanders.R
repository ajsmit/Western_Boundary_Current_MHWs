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
# MKE mask and the 75th perc. MKE mask to define the region within which 
# we want to compare MKE and MHW intensity

# With the data masked the co-occurences and correlations may then be run

# NB: This script is intended to be run via an R terminal
  # Therefore each different section should be (un)commented when desired


# Functions ---------------------------------------------------------------

source("setup/meanderFunctions.R")


# Subet AVISO and calculate KE --------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Run sequentially
# for(i in 1:(ncol(bbox)-1)){
# 
#   # Determine region
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
# 
#   # Create the masks
#   AVISO_KE_save(colnames(bbox)[i])
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# 
#   # Clear up some RAM
#   gc()
# } # ~ 13 minute each


# Subsetting MHW results --------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Subset MHW data
# for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   MHW_sub_save(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
#   gc()
# } ~ 10 minutes each


# Percentile masks --------------------------------------------------------

# Set cores
library(doMC); doMC::registerDoMC(cores = 50)

# Run sequentially
# for(i in 1:(ncol(bbox)-1)){
# 
#   # Determine region
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
# 
#   # Load the desired files
#   AVISO_KE <- fread(paste0("../data/WBC/AVISO_KE_",region,".csv"))
#   MHW_event <- fread(paste0("../data/WBC/MHW_event_",region,".csv"))
#   print("Data loaded")
# 
#   # Create the masks
#   masks(AVISO = AVISO_KE, MHW = MHW_event)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# 
#   # Clear up some RAM
#   rm(AVISO_KE, MHW_event); gc()
# } # ~ 1 minute each


# Correlate pixels --------------------------------------------------------

# In this final step we take only the MKE and MHW pixels that are in the 
# 50th percentile mask but not in the 90th percentile mask.
# We also take the pixels that are in the 90th percentile max intensity mask.
# We then take the days the pixels in the 50th perc. MKE and 90th perc. max int.
# areas are themselves in the 90th percentile for MKE.
# A count of how often this occurs when MHWs are occurring is made to find
# what the proportion of this occurence is.
# The MKE and mean intensities on days when MKE is in the 90th perc. 
# are then correlated.
# This result will help to illustrate the potential relationship between
# meanders and MHWs.

# Set cores
library(doMC); doMC::registerDoMC(cores = 5)

# Calculate the results in serial
# for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   meander_co_calc(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
#   gc()
# } # ~xxx seconds each


# meander_co_calc(colnames(bbox)[1])
# meander_co_calc(colnames(bbox)[2])
# meander_co_calc(colnames(bbox)[3])
# meander_co_calc(colnames(bbox)[4])
# meander_co_calc(colnames(bbox)[5])


# Visualise results -------------------------------------------------------

# A for loop for ease of creating all desired figures
for(i in 1:(ncol(bbox)-1)){
  region <- colnames(bbox)[i]
  print(paste0("Began run on ",region," at ",Sys.time()))
  meander_vis(region)
  print(paste0("Finished run on ",region," at ",Sys.time()))
  gc()
} # ~xxx seconds each

# meander_vis(colnames(bbox)[1])
# meander_vis(colnames(bbox)[2])
# meander_vis(colnames(bbox)[3])
# meander_vis(colnames(bbox)[4])
# meander_vis(colnames(bbox)[5])

