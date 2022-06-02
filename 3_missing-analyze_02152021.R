# Analyze data with simulated missingness
# Author: Daniel Vader
# Date created: 12/7/2020

# Set path to local library for running extra packages on the cluster
.libPaths(c(.libPaths(), "rlibs"))

library(tidyverse)
library(mice)
library(micemd)
library(survival)
library(pbapply)
library(twang) # Used for gradient boosted ps model

################################################################################
# Generate KM Curve with complete data
################################################################################
library(survminer)
library(survival)
library(ggquickeda)

sdat.p <- readRDS("Y:/programs/dan/analyticdat/baseccdat.rds") %>%
  mutate(cmonth = ifelse(cmonth > 60, 60, cmonth),
         Treatment = factor(treat, levels=c(0,1), labels=c("Chemotherapy", "Immunotherapy"))) # 5 years

sfit <- survfit(Surv(cmonth, dead) ~ treat, data=sdat.p)


tiff("Y:/programs/dan/figures/survival_base-data.tiff", 
     width=5, height=2.5, units="in", res=300)
ggplot(sdat.p, aes(time=cmonth, 
                   status = dead, 
                   linetype=Treatment, 
                   #color=treatc,
                   group=Treatment)) +
  geom_km() + 
  #geom_kmband() +
  ylim(0,1) +
  xlab("Time (months)") +
  ylab("Survival probability") +
  theme_bw()
dev.off()

plot(sfit,
     lty=c("solid", "dashed"),
     xlab="Months",
     ylab="Survival probability")
legend("topright", c("Non-immunotherapy", "Immunotherapy"),
       lty=c("solid", "dashed"))

################################################################################
# Run models using different missing data approaches: CC, MI, and PSC
################################################################################

# Load data
load("Y:/programs/dan/analyticdat/simdat1000.Rdata")

# Load function to run analyses across analytic data sets
source("Y:/programs/dan/functions/fanalyze.R")


# No misclassification #########################################################

# Complete case
eff.cc <- md.analyze(missing.sim, approach="cc")
save(eff.cc, file= paste0("Y:/programs/dan/output/cc-effects_",Sys.Date(), ".Rdata"))

# Propensity score calibration
eff.psc <- md.analyze(missing.sim, approach="psc")
save(eff.psc, file= paste0("Y:/programs/dan/output/psc-effects_",Sys.Date(), ".Rdata"))

# Complete data (unadjusted)
eff.cd <- md.analyze(missing.sim, approach="nomiss")
save(eff.cd, file= paste0("Y:/programs/dan/output/cd-effects_",Sys.Date(), ".Rdata"))

# Complete data (adjusted)
eff.cdc <- md.analyze(missing.sim, approach="nomissc")
save(eff.cdc, file= paste0("Y:/programs/dan/output/cdc-effects_",Sys.Date(), ".Rdata"))

# Non-differential misclassification ###########################################

# Complete case
eff.cc.mc <- md.analyze(missing.sim, approach="cc", misclass=T)
save(eff.cc.mc, file= paste0("Y:/programs/dan/output/cc-effects-mc_",Sys.Date(), ".Rdata"))

# Propensity score calibration
eff.psc.mc <- md.analyze(missing.sim, approach="psc", misclass=T)
save(eff.psc.mc, file= paste0("Y:/programs/dan/output/psc-effects-mc_",Sys.Date(), ".Rdata"))

# Complete data (unadjusted)
eff.cd.mc <- md.analyze(missing.sim, approach="nomiss", misclass=T)
save(eff.cd.mc, file= paste0("Y:/programs/dan/output/cd-effects-mc_",Sys.Date(), ".Rdata"))

# Complete data (adjusted)
eff.cdc.mc <- md.analyze(missing.sim, approach="nomissc", misclass=T)
save(eff.cdc.mc, file= paste0("Y:/programs/dan/output/cdc-effects-mc_",Sys.Date(), ".Rdata"))

# Multiple imputation
# Best to run this one on the cluster, currently set to 5 cores.
# eff.mi <- md.analyze(f.list, approach="mi", misclass=T)
# save(eff.mi, file= paste0("Y:/programs/dan/output/mi-effects_",Sys.Date(), ".Rdata"))
# print("done")

# Differential misclassification ###########################################

# Complete case
eff.cc.dmc <- md.analyze(missing.sim, approach="cc", misclass=T, diff.miss=T)
save(eff.cc.dmc, file= paste0("Y:/programs/dan/output/cc-effects-dmc_",Sys.Date(), ".Rdata"))

# Propensity score calibration
eff.psc.dmc <- md.analyze(missing.sim, approach="psc", misclass=T, diff.miss=T)
save(eff.psc.dmc, file= paste0("Y:/programs/dan/output/psc-effects-dmc_",Sys.Date(), ".Rdata"))

# Complete data (unadjusted)
eff.cd.dmc <- md.analyze(missing.sim, approach="nomiss", misclass=T, diff.miss=T)
save(eff.cd.dmc, file= paste0("Y:/programs/dan/output/cd-effects-dmc_",Sys.Date(), ".Rdata"))

# Complete data (adjusted)
eff.cdc.dmc <- md.analyze(missing.sim, approach="nomissc", misclass=T, diff.miss=T)
save(eff.cdc.dmc, file= paste0("Y:/programs/dan/output/cdc-effects-dmc_",Sys.Date(), ".Rdata"))

# Multiple imputation
# Best to run this one on the cluster, currently set to 5 cores.
# eff.mi <- md.analyze(f.list, approach="mi")
# save(eff.mi, file= paste0("Y:/programs/dan/output/mi-effects_",Sys.Date(), ".Rdata"))
# print("done")

# Doubled sample size, no misclassification ####################################
load("Y:/programs/dan/analyticdat/simdat1000d.Rdata")

# Complete case
eff.cc.d <- md.analyze(missing.sim.d, approach="cc")
save(eff.cc.d, file= paste0("Y:/programs/dan/output/cc-effects-d_",Sys.Date(), ".Rdata"))

# Propensity score calibration
eff.psc.d <- md.analyze(missing.sim.d, approach="psc")
save(eff.psc.d, file= paste0("Y:/programs/dan/output/psc-effects-d_",Sys.Date(), ".Rdata"))

# Complete data (unadjusted)
eff.cd.d <- md.analyze(missing.sim.d, approach="nomiss")
save(eff.cd.d, file= paste0("Y:/programs/dan/output/cd-effects-d_",Sys.Date(), ".Rdata"))

# Complete data (adjusted)
eff.cdc.d <- md.analyze(missing.sim.d, approach="nomissc")
save(eff.cdc.d, file= paste0("Y:/programs/dan/output/cdc-effects-d_",Sys.Date(), ".Rdata"))


################################################################################
# Check runtimes
################################################################################

load("Y:/programs/dan/analyticdat/simdat1.Rdata")
eff.cc.t <- md.analyze(missing.sim.1, approach="cc", timer=T)
eff.psc.t <- md.analyze(missing.sim.1, approach="psc", timer=T)
eff.mi.t <- md.analyze(missing.sim.1, approach="mi", timer=T)

tlist <- list(eff.cc.t, eff.mi.t, eff.psc.t)
save(tlist, file=paste0("Y:/programs/dan/output/model-time-elapsed_",
                        Sys.Date(), 
                        ".Rdata"))

durtab <- delist_dur(tlist) %>% 
  filter(mech == "MAR") %>%
  select(-mech) %>%
  group_by(sim.hr, prop, approach) %>% 
  tidyr::pivot_wider(names_from=approach, 
                     values_from=c(dur)) %>%
  relocate(sim.hr, prop, CC, MI, PSC)

write.csv(durtab, file=paste0("Y:/programs/dan/tables/model-time-elapsed-tab_",
                              Sys.Date(), 
                              ".csv"))


################################################################################
# Simulate "truth." Run analyses with complete data from large simulation set.
################################################################################

# Simulate truth using logistic regression ps
sdat.c <- readRDS("Y:/programs/dan/analyticdat/baseccdat.rds")
tsim.dat <- readRDS("Y:/programs/dan/plasdatm.Rdata")

teff10 <- tsim(tsim.dat$sor10, sdat.c)
teff09 <- tsim(tsim.dat$sor09, sdat.c)
teff05 <- tsim(tsim.dat$sor05, sdat.c)
teff.all <- list()
teff.all[["sor10"]] <- teff10
teff.all[["sor09"]] <- teff09
teff.all[["sor05"]] <- teff05
saveRDS(teff.all, paste0("Y:/programs/dan/output/tr-effects_", Sys.Date(), ".rds"))

# Simulate truth with gradient boosted ps (use separate script to run on the grid)
# set.seed(324769)
# 
# teffb10 <- tsim(tsim.dat$sor10, sdat.c, boosted = T)
# teffb09 <- tsim(tsim.dat$sor09, sdat.c, boosted = T)
# teffb05 <- tsim(tsim.dat$sor05, sdat.c, boosted = T)
# teffb.all <- list()
# teffb.all[["sor10"]] <- teffb10
# teffb.all[["sor09"]] <- teffb09
# teffb.all[["sor05"]] <- teffb05
# saveRDS(teffb.all, paste0("output/tr-effects-boosted_", Sys.Date(), ".rds"))

################################################################################
# Summarize bias
################################################################################
#library(simhelpers) # https://cran.r-project.org/web/packages/simhelpers/vignettes/MCSE.html
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)

source("functions/fanalyze.R")

eff.tr <- readRDS("output/tr-effects_2021-11-29.rds")
eff.tr2 <- readRDS("output/tr-effects-boosted_2021-11-29.rds")

wdate <- "2022-01-12"
wdate2 <- "2021-11-22"
wdate3 <- "2021-11-29"

load(paste0("output/cc-effects_", wdate, ".Rdata"))
load(paste0("output/psc-effects_", wdate, ".Rdata"))
load(paste0("output/mi-effects_", "2021-11-27", ".Rdata"))
load(paste0("output/cd-effects_", wdate2, ".Rdata"))
load(paste0("output/cdc-effects_", wdate2, ".Rdata"))

eff.cc.d <- loadRData(paste0("output/cc-effects-d_", wdate, ".Rdata"))
eff.psc.d <- loadRData(paste0("output/psc-effects-d_", wdate, ".Rdata"))
eff.mi.d <- loadRData(paste0("output/mi-effects-d_", "2021-11-24", ".Rdata"))
eff.cd.d <- loadRData(paste0("output/cd-effects-d_", wdate3, ".Rdata"))
eff.cdc.d <- loadRData(paste0("output/cdc-effects-d_", wdate3, ".Rdata"))

load(paste0("output/cc-effects-mc_", wdate, ".Rdata"))
load(paste0("output/psc-effects-mc_", wdate, ".Rdata"))
eff.mi.mc <- loadRData(paste0("output/mi-effects-mc_", "2021-11-23", ".Rdata"))
load(paste0("output/cd-effects-mc_", wdate2, ".Rdata"))
load(paste0("output/cdc-effects-mc_", wdate2, ".Rdata"))

load(paste0("output/cc-effects-dmc_", "2022-02-08", ".Rdata"))
load(paste0("output/psc-effects-dmc_", "2022-02-08", ".Rdata"))
eff.mi.dmc <- loadRData(paste0("output/mi-effects-dmc_", "2022-02-12", ".Rdata"))
load(paste0("output/cd-effects-dmc_", "2022-02-09", ".Rdata"))
load(paste0("output/cdc-effects-dmc_", "2022-02-09", ".Rdata"))

# Transform data from lists into data frame (tibble) for ease of use ###########

# Primary analysis
# With truth estimated using regular PS model
effs.reg <- delist_all(eff.cc, eff.psc, eff.mi, eff.cd, eff.cdc, eff.tr)
# With truth estimated using boosted CART PS model
effs.boost <- delist_all(eff.cc, eff.psc, eff.mi, eff.cd, eff.cdc, eff.tr2)

# Double sample size
effs.reg.d <- delist_all(eff.cc.d, eff.psc.d, eff.mi.d, eff.cd.d, eff.cdc.d, eff.tr)
effs.boost.d <- delist_all(eff.cc.d, eff.psc.d, eff.mi.d, eff.cd.d, eff.cdc.d, eff.tr2)

# Non-differntial misclassification
effs.reg.mc <- delist_all(eff.cc.mc, eff.psc.mc, eff.mi.mc, eff.cd.mc, eff.cdc.mc, eff.tr)
effs.boost.mc <- delist_all(eff.cc.mc, eff.psc.mc, eff.mi.mc, eff.cd.mc, eff.cdc.mc, eff.tr2)

# Differential misclassification
effs.reg.dmc <- delist_all(eff.cc.dmc, eff.psc.dmc, eff.mi.dmc, eff.cd.dmc, eff.cdc.dmc, eff.tr)
effs.boost.dmc <- delist_all(eff.cc.dmc, eff.psc.dmc, eff.mi.dmc, eff.cd.dmc, eff.cdc.dmc, eff.tr2)

# Generate box plots of estimates ##############################################
library(viridis)
library(ggplot2)
fht <- 7 # figure height (in)
fwh <- 7 # figure width (in)
fres <- 300 # figure res (ppi)

# # Plot effects with GLM PS truth
# ep <- plot.effects(eff.all)
# tiff("figures/eff_by_model.tiff", width=fwh, height=fht, units="in", res=fres)
# ep
# dev.off()
# 
# # Plot bias with GLM PS truth
# bp <- plot.bias(eff.all)
# tiff("figures/bias_by_model.tiff", width=fwh, height=fht, units="in", res=fres)
# bp
# dev.off()

# # Plot effects with boosted CART PS truth
# ep.b <- plot.effects(eff.all.b)
# tiff("figures/eff_by_model_bt.tiff", width=fwh, height=fht, units="in", res=fres)
# ep.b
# dev.off()

# Plot bias with boosted CART PS truth
bp.b <- plot.bias(effs.boost$eff.all, alt=T)
tiff("figures/bias_by_modela.tiff", width=fwh, height=fht, units="in", res=fres)
bp.b
dev.off()

png("figures/bias_by_modela.png", width=5, height=8, units="in", res=fres)
bp.b + theme(legend.position="bottom", 
             legend.direction="vertical",
             legend.title=element_text(size=12), 
             legend.text=element_text(size=11))
dev.off()

# Plot bias with double sample size (boosted PS)
bp.b.d <- plot.bias(effs.boost.d$eff.all, alt=T)
tiff("figures/bias_by_model_da.tiff", width=fwh, height=fht, units="in", res=fres)
bp.b.d
dev.off()

# Plot bias with non-differential misclassificaiton
bp.b.mc <- plot.bias(effs.boost.mc$eff.all, alt=T)
tiff("figures/bias_by_model_mca.tiff", width=fwh, height=fht, units="in", res=fres)
bp.b.mc
dev.off()
png("figures/bias_by_model_mca.png", width=fwh, height=fht, units="in", res=fres)
bp.b.mc
dev.off()

# Plot bias with differential misclassification
bp.b.dmc <- plot.bias(effs.boost.dmc$eff.all, alt=T)
tiff("figures/bias_by_model_dmca.tiff", width=fwh, height=fht, units="in", res=fres)
bp.b.dmc
dev.off()
png("figures/bias_by_model_dmca.png", width=fwh, height=fht, units="in", res=fres)
bp.b.dmc
dev.off()


# Report bias, 95% CI coverage proportions, mse, 
    # empircal standard error, and model-based standard-error.
# bt <- tab.bias(eff.all.cd)
# write.csv(bt, "tables/model-stats.csv", na="")

bt <- tab.bias(effs.boost$eff.all.cd)
write.csv(bt, "tables/model-stats_b.csv", na="")

bt.d <- tab.bias(effs.boost.d$eff.all.cd)
write.csv(bt.d, "tables/model-stats_d.csv", na="")

bt.mc <- tab.bias(effs.boost.mc$eff.all.cd)
write.csv(bt.mc, "tables/model-stats_mc.csv", na="")

bt.dmc <- tab.bias(effs.boost.dmc$eff.all.cd)
write.csv(bt.dmc, "tables/model-stats_dmc.csv", na="")


# Don't need to repeat effects since they don't use truth model
#bt <- tab.eff(eff.all.cd)
#write.csv(bt, "tables/model-stats-eff.csv", na="")

bt.e <- tab.eff(effs.boost$eff.all)
write.csv(bt.e, "tables/model-stats-eff_b.csv", na="")

bt.e.d <- tab.eff(effs.boost.d$eff.all)
write.csv(bt.e.d, "tables/model-stats-eff_d.csv", na="")

bt.e.mc <- tab.eff(effs.boost.mc$eff.all)
write.csv(bt.e.mc, "tables/model-stats-eff_mc.csv", na="")

bt.e.dmc <- tab.eff(effs.boost.dmc$eff.all)
write.csv(bt.e.dmc, "tables/model-stats-eff_dmc.csv", na="")

# Plot 95% Coverage probs & MSE
mse95p <- plot.mse95(effs.boost)
tiff("figures/mse_by_model.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p$mse
dev.off()
tiff("figures/95cip_by_model.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p$cover95
dev.off()
png("figures/mse_by_model.png", width=fwh, height=fht, units="in", res=fres)
mse95p$mse
dev.off()
png("figures/95cip_by_model.png", width=6, height=5, units="in", res=fres)
mse95p$cover95 + coord_cartesian(ylim=c(0.7, 1.0))
dev.off()



mse95p.d <- plot.mse95(effs.boost.d, lim1=T)
tiff("figures/mse_by_model_d.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.d$mse
dev.off()
png("figures/mse_by_model_d.png", width=fwh, height=fht, units="in", res=fres)
mse95p.d$mse
dev.off()
tiff("figures/95cip_by_model_d.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.d$cover95
dev.off()
png("figures/95cip_by_model_d.png", width=fwh, height=fht, units="in", res=fres)
mse95p.d$cover95
dev.off()

mse95p.mc <- plot.mse95(effs.boost.mc)
tiff("figures/mse_by_model_mc.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.mc$mse
dev.off()
png("figures/mse_by_model_mc.png", width=fwh, height=fht, units="in", res=fres)
mse95p.mc$mse
dev.off()
tiff("figures/95cip_by_model_mc.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.mc$cover95
dev.off()
png("figures/95cip_by_model_mc.png", width=fwh, height=fht, units="in", res=fres)
mse95p.mc$cover95
dev.off()

mse95p.dmc <- plot.mse95(effs.boost.dmc)
tiff("figures/mse_by_model_dmc.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.dmc$mse
dev.off()
png("figures/mse_by_model_dmc.png", width=fwh, height=fht, units="in", res=fres)
mse95p.dmc$mse
dev.off()
tiff("figures/95cip_by_model_dmc.tiff", width=fwh, height=fht, units="in", res=fres)
mse95p.dmc$cover95
dev.off()
png("figures/95cip_by_model_dmc.png", width=fwh, height=fht, units="in", res=fres)
mse95p.dmc$cover95
dev.off()

