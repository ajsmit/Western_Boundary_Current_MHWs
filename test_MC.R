library(data.table)
# devtools::install_github("hadley/multidplyr")
library(tidyverse)
library(multidplyr)
library(heatwaveR)
library(fasttime)

inFile <- "/Volumes/Benguela/spatial/processed/OISSTv2/WBC/daily/EAC-avhrr-only-v2.19810901-20180930.csv"

sst <- fread(inFile,
             colClasses = list("numeric" = 1:3, "Date" = 4),
             col.names = c("lon", "lat", "sst", "date"))
sst[,date := as.Date(fastPOSIXct(date, tz = NULL))]

cl <- future::makeClusterPSOCK(worker = "Darwin.local", 2L, timeout = 30, verbose = TRUE)
cl <- parallel::makeCluster(4L, timeout = 30, verbose = TRUE)

registerDoParallel(cores = 4)
getDoParWorkers()

sst %>%
  dplyr::filter(lon <= 151 & lat <= -41) %>%
  partition(lon, lat, cluster = cl) %>%
  dplyr::summarize(sst_mean = mean(sst))

library("future")
plan(sequential)
demo("mandelbrot", package = "future", ask = FALSE)

plan(multiprocess)
demo("mandelbrot", package = "future", ask = FALSE)

library("listenv")
plan(list(multiprocess, sequential))

htseq_align <- function(fq, chr) { chr }

fqs <- dir(pattern = "[.]fastq$")

bams <- listenv()
for (ss in seq_along(fqs)) {
  fq <- fqs[ss]
  bams[[ss]] %<-% {
    bams_ss <- listenv()
    for (cc in 1:24) {
      bams_ss[[cc]] %<-% htseq_align(fq, chr = cc)
    }
    as.list(bams_ss)
  }
}
bams <- as.list(bams)

# doParallel
library(doParallel)
# cl <- makeCluster(spec = 2, type = "PSOCK") # doesn't work
# registerDoParallel(cl, cores = 4) # doesn't work
numCores <- detectCores()
getDoParWorkers()
getDoParName()
getDoParVersion()
registerDoParallel(cores = numCores)
system.time(foreach(i=1:100000) %dopar% sum(tanh(1:i)))

# another example
getPrimeNumbers <- function(n) {
  n <- as.integer(n)
  if(n > 1e6) stop("n too large")
  primes <- rep(TRUE, n)
  primes[1] <- FALSE
  last.prime <- 2L
  for(i in last.prime:floor(sqrt(n)))
  {
    primes[seq.int(2L*last.prime, n, last.prime)] <- FALSE
    last.prime <- last.prime + min(which(primes[(last.prime+1):n]))
  }
  which(primes)
}

# below doesn work...
library(doParallel)
# with parLapply ... faster on Lengau
no_cores <- detectCores() - 1
registerDoParallel(cores = no_cores)
cl <- makeCluster(no_cores, type = "FORK")
system.time(result <- parLapply(cl, 10:10000, getPrimeNumbers)  )
stopCluster(cl)

# with foreach
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cl)
system.time(result <- foreach(i = 10:10000) %dopar% getPrimeNumbers(i))
stopCluster(cl)
