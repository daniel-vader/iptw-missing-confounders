# Set path to local library for running extra packages on the cluster
.libPaths(c(.libPaths(), "rlibs"))

library(tidyverse)
library(mice)
library(micemd)
library(survival)
library(pbapply)

################################################################################
# Run models using different missing data approaches: CC, MI, and PSC
################################################################################

# Load data
load("../analyticdat/simdat1000d.Rdata")

# Load function to run analyses across analytic data sets
source("../functions/fanalyze.R")

# Run MI and save
eff.mi <- md.analyze(missing.sim.d, approach="mi", grid=T)
save(eff.mi, file= paste0("../output/mi-effects-d_",Sys.Date(), ".Rdata"))
print("done")