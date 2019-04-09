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

# Then the eddy trajectory product from AVISO is used to create an eddy mask

# With all of the results ready, the next step is to use the 90th perc.
# MKE mask and the 50th perc. MKE mask to define the region within which 
# we want to compare MKE and MHW intensity

# We then mask these results with the eddy trajectory mask

# With the data masked in two progressive steps the co-occurences and correlations may then be run

# NB: This script is intended to be run via an R terminal
  # Therefore each different section should be (un)commented when desired


# Functions ---------------------------------------------------------------

source("setup/meanderFunctions.R")


# Subet AVISO and calculate KE --------------------------------------------

# Set cores
doMC::registerDoMC(cores = 50)

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
doMC::registerDoMC(cores = 50)

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
doMC::registerDoMC(cores = 50)

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
#   mke_masks(AVISO = AVISO_KE, MHW = MHW_event)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# 
#   # Clear up some RAM
#   rm(AVISO_KE, MHW_event); gc()
# } # ~ 1 minute each


# Eddy masks --------------------------------------------------------------

# In this step we run a function that looks at the AVISO eddy track database
# and finds which pixels in each bbox are within the radius and centre point
# of all of the recorded eddies from 1993-01-01 to 2019-
doMC::registerDoMC(cores = 50)
# eddy_masks()


# Mask regions ------------------------------------------------------------

# In this step we take only the MKE and MHW pixels that are in the 
# 50th percentile mask but not in the 90th percentile mask.
# We also take the pixels that are in the 90th percentile max intensity mask.
# Set cores
doMC::registerDoMC(cores = 50)

# Run sequentially
# for(i in 1:(ncol(bbox)-1)){
# 
#   # Determine region
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
# 
#   # Mask the regions
#      # NB: This runs the MKE and Eddy masks
#   mask_region(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
# 
#   # Clear up some RAM
#   gc()
# } # ~ 8 minutes each


# Correlate pixels --------------------------------------------------------

# Here we take the days the pixels in the 50th perc. MKE and 90th perc. max int.
# areas are themselves in the 90th percentile for MKE.
# A count of how often this occurs when MHWs are occurring is made to find
# what the proportion of this occurence is.
# The MKE and mean intensities on days when MKE is in the 90th perc. 
# are then correlated.
# This result helps to illustrate the potential relationship between
# meanders and MHWs.

# Set cores
doMC::registerDoMC(cores = 50)

# Calculate the results in serial
for(i in 1:(ncol(bbox)-1)){ # Not running for Benguela
  region <- colnames(bbox)[i]
  print(paste0("Began run on ",region," at ",Sys.time()))
  meander_co_calc(region)
  print(paste0("Finished run on ",region," at ",Sys.time()))
  gc()
} # ~30 seconds each


# Visualise results -------------------------------------------------------

# Here is the naming convention for the figures produced below:

# region_finalStatistic_pixelFilter_timeUnit.pdf

# region: the acronym for the study region of interest

# finalStatistic: The result being shown. Presently all of the figures I've uploaded
# will say "cooc" here. 
# This stands for co-occurrence and the figures show the rate of co-occurence 
# of the proportion of days when MWHs are occurring in a pixel against the total
# number of days in that pixel (0 = none, 1 = perfect) 

# pixelFilter: This will say either "50" or "max"
# 50 means that only pixels that are within the 50th to 90th percentile for 
# mean MKE in the region are being shown
# max means that only pixels that have had events whose max intensity is in
# the 90th percentile for the study region are being shown

# timeUnit: This shows if the figure is showing the "total" number of days
# or if it is showing "month" facets.
# Because there was no apparent monthly pattern I didn't save all of the 
# figures as I didn't want to create too much clutter
# The code to create them exists in "meander_vis()" in "setup/meanderFunctions.R"
# and must just be uncommented to run

# Within each figure are two panels labeled "cooc_90" and "cooc_flat"
# cooc_flat shows the co-occurrence rates for days with MHWs against 
# all days (~9000) in a given pixel. Generally this is around 1 - 3%
# cooc_90 shows the rates of co-occurrence between MHW days and days when
# the MKE is in the 90th percentile for the entire study area. These rates
# tend to be much higher, some of them are even 1.0

# A for loop for ease of creating all desired figures
# for(i in 1:(ncol(bbox)-1)){
#   region <- colnames(bbox)[i]
#   print(paste0("Began run on ",region," at ",Sys.time()))
#   meander_vis(region)
#   print(paste0("Finished run on ",region," at ",Sys.time()))
#   gc()
# } # ~1 second each

