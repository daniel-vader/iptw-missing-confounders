.libPaths(c(.libPaths(), "rlibs"))

library(dplyr)
library(survival)
library(twang)
library(parallel)

source("functions/fanalyze.R")

set.seed(324769)

sdat.c <- readRDS("analyticdat/baseccdat.rds")
tsim.dat <- readRDS("plasdatm.Rdata")

teffb.all <- mclapply(tsim.dat, FUN=tsim, main = sdat.c, boosted = T, mc.cores = 3)

saveRDS(teffb.all, paste0("output/tr-effects-boosted_", Sys.Date(), ".rds"))