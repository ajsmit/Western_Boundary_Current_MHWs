# The purpose of this script is to first create the AVISO and MHW databases
# used throughout this analysis because this script is to be run on the 
# tikoraluk server and the database files are too large to send via
# conventional digital means
# These files are created here by:
  # 1) Loading the bounding box definitions
  # 2) Loading the AVISO grid conversion mask
  # 3) Subsetting from the MHW and AVISO databases on tikoraluk

# With these databases recreated here, the next step is to re-calculate EKE,
# which can be done with the existing code in xxxxx.R

# With all of the results ready, the next step is to use the 90th perc.
# EKE mask and the 75th perc. EKE mask to define the region within which 
# we want to compare poleward velocity (U+V) and mean intensity

