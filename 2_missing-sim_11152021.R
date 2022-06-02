# Generate simulated data sets for missing data project
# Author: Daniel Vader
# Date created: 08/21/2020

# load packages
library(Plasmode)
library(tidyverse)
library(survival)
library(naniar)
#library(splines)
library(mice)

# Optional stages for debuging / descriptive stats
assess.missing = F

# set number of simulations
nsim = 1000

# set seed
# sample(1:100000, 1) # 63269
set.seed(63269)

################################################################################
# Load data and select variables for analysis
  # Outcome is overall survival using treatment start -> death
  # Primary exposure is immunotherapy vs chemotherapy 
  # (should we do carboplatin vs immunotherapy?).
  # Covariables should be selected to build structure for the simulated exposure
  # + outcome and to be included in the propensity score adjustment.
################################################################################
sdat <- read.csv("Y:/programs/dan/analyticdat/deident_comorb_flatiron_ucc_2021-11-03.csv") %>%
  # what variables are we using? Other project recommends albumin, hemoglobin,
  # number of mets, and BMI
  dplyr::mutate(
    age = baseline_year-birthyear,
    mage = log(age) - min(log(age)), # potential age val for ampute weights
    # Change unknown to NA
    race.ethnicity = as.factor(ifelse(race.ethnicity == "Unknown", NA, race.ethnicity)),
    diseasegrade = as.factor(ifelse(diseasegrade == "Unknown/not documented", NA, diseasegrade)),
    denovomet = ifelse(groupstage %in% c("Stage IV", "Stage IVA", "Stage IVB"), 1,
                       ifelse(groupstage == "Unknown/not documented", NA, 0)),
    smokingstatus = as.factor(ifelse(smokingstatus == "Unknown", NA, smokingstatus)),
    #groupstage = as.factor(ifelse(groupstage == "Unknown/not documented", NA, groupstage)),
    # Inverse of dead for plasmode sim
    ndead = 1-dead,
    # Set treatment class to immuno vs non-immuno
    treat = ifelse(trtclass1 == "immuno", 1, 0),
    cmonth = lastcontactmonth + 1, # censoring month including death
    cmonth2 = censormonth + 1, # censoring month excluding death
    # Create dummy vars
    reth_black = ifelse(race.ethnicity == "Black", 1, 0),
    reth_hisp = ifelse(race.ethnicity == "Hispanic or Latino", 1, 0),
    reth_oth = ifelse(race.ethnicity == "Other", 1, 0),
    genderf = ifelse(gender == "F", 1, 0),
    #practypec = ifelse(practicetype == "COMMUNITY", 1, 0),
    smokey = ifelse(smokingstatus == "History of smoking", 1, 0),
    dgradeh = ifelse(diseasegrade == "High grade (G2/G3/G4)", 1, 0),
    surgery = ifelse(surgery == T, 1, 0),
    site_ureter = ifelse(primarysite == "Ureter", 1, 0), # ref = bladder
    site_renal = ifelse(primarysite == "Renal Pelvis", 1, 0),
    site_urethra = ifelse(primarysite == "Urethra", 1, 0),
    # ecog01 = ifelse(b.ecogvalue < 1, 1, 0), # ref = ecog < 1
    ecog12 = ifelse(b.ecogvalue >= 1 & b.ecogvalue < 2, 1, 0),
    ecog24 = ifelse(b.ecogvalue >= 2, 1, 0)
    
  ) %>%
  dplyr::select(patientid, genderf, 
                reth_black, reth_hisp, reth_oth,
                age, smokey,
                surgery, 
                # site_ureter, site_renal, site_urethra, 
                dead, ndead, cmonth,
                treat, ecog12, ecog24
                # elixhauser score vars
                #chf, carit, valv, pcd, pvd, hypunc, hypc, para, ond, cpd,
                #diabunc, diabc, hypothy, rf, ld, pud, aids, lymph, metacanc, 
                #solidtum, rheumd, coag, obes, wloss, fed, blane, dane, alcohol, 
                #drug, psycho, depre
  ) %>% relocate(patientid, dead, ndead, cmonth, treat)

# Remove incomplete data. N = 4361 (Previous N = 7309)
sdat.c <- sdat %>% na.omit()

# save
saveRDS(sdat.c, "Y:/programs/dan/analyticdat/baseccdat.rds")

################################################################################
# Generate baseline descriptive stats
################################################################################

# Create baseline data table ###################################################
btab.dat <- sdat %>% 
  mutate(gender = factor(genderf, levels=c(0,1), labels=c("Male", "Female")),
         reth = ifelse(reth_hisp == 1, 2, 
                       ifelse(reth_black == 1, 1,
                              ifelse(reth_oth == 1, 4, 3))),
         reth = factor(reth, levels=c(1,2,3,4), labels=c("Black", "Hispanic", "White", "Other")),
         ecog = ifelse(ecog12 == 1, 2, ifelse(ecog24 == 1, 3, 1)),
         ecog = factor(ecog, levels=c(1,2,3), labels=c("0", "1", "2-4")),
         #grade = factor(dgradeh, levels=c(0,1), labels=c("< G2", "High (G2/G3/G4)")),
         surgery = factor(surgery, levels=c(0,1), labels=c("No", "Yes")),
         eversmoke = factor(smokey, levels=c(0,1), labels=c("No", "Yes")),
         immunotherapy = factor(treat, levels=c(1,0), labels=c("Immunotherapy", "Other"))
  ) %>%
  select(immunotherapy, gender, reth, ecog, 
         #grade, 
         surgery, eversmoke, age, cmonth)

# Initiate table generating function
btab.gen <- function(dat, catvars, numvars, una="ifany"){
  tabdat <- tibble(var = character(),
                   cat = character(),
                   ntot = numeric(), 
                   ptot = numeric(),
                   qtot = numeric())
  
  for(i in 1:length(catvars)){
    v <- dat %>% pull(catvars[i])
    tn <- table(v, useNA = una)
    tp <- prop.table(tn)
    
    
    for(j in 1:length(tn)){
      tabdat <- tabdat %>% 
        add_row(var=catvars[i], cat=as.character(names(tn)[j]), 
                ntot=tn[j], ptot=tp[j], qtot=NA)
    }
    if(una == "no"){
      ttotna <- table(is.na(v))
      pttot <- prop.table(ttotna)
      
      if(sum(is.na(v)) > 0){
        tabdat <- tabdat %>% 
          add_row(var=catvars[i], cat="Missing", 
                  ntot=ttotna[2], ptot=NA,#pnat, 
                  qtot=NA)
      }
    }
  }
  
  for(i in 1:length(numvars)){
    v <- dat %>% pull(numvars[i])
    
    mtot <- median(v, na.rm=T)
    lqtot <- quantile(v, 0.25, na.rm=T)
    uqtot <- quantile(v, 0.75, na.rm=T)
    
    tabdat <- tabdat %>% 
      add_row(var=numvars[i], cat="median(IQR)", 
              ntot=mtot, ptot=lqtot, qtot=uqtot)
    
    ttot <- table(is.na(v))
    pttot <- prop.table(ttot)
    
    if(length(ttot) > 1){
      nat <- ttot[2]
      pnat <- pttot[2]
      
      tabdat <- tabdat %>% 
        add_row(var=numvars[i], cat="Missing", 
                ntot=nat, ptot=NA,#pnat, 
                qtot=NA)
    }
  }
  return(tabdat)
}

cvars <- c("immunotherapy", "gender", "reth", "ecog", "surgery", "eversmoke")
nvars <- c("age", "cmonth")

t1.all <- btab.gen(btab.dat, cvars, nvars, una="no")
write.csv(t1.all, "tables/t1-base-all.csv", na="")

t1.complete <- btab.gen(na.omit(btab.dat), cvars, nvars, una="always")
write.csv(t1.complete, "tables/t1-base-complete.csv", na="")


# glm(treat ~ genderf +
#       reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#       smokey + surgery + age, 
#     family=binomial(link="logit"),
#     data=dc)
# 
# # Generate baseline PS model ###################################################
# glm(treat ~ genderf +
#       reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#       smokey + surgery + age, 
#     family=binomial(link="logit"),
#     data=sdat.c)

################################################################################
# Check NA proportions 
################################################################################
# ECOG is the main reason for missing data. 
if(assess.missing){
  natab <- function(x){
    if(anyNA(x)){
      table(is.na(x))
    }else(NA)
  }
  na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
  apply(sdat, 2, natab) %>% na.omit.list()
  
  # Missing data visualization (generate and output to file)
  # 1778 missing due to ECOG alone, 384 missing due to disease grade alone,
  # 264 due to race/ethnicity alone, 465 missing due to ECOG with others,
  # # 34 due to race/eth + disease grade combo, 1 missing due to gender
  p1 <- gg_miss_upset(sdat)
  tiff("U:/projects/FLATIRON missing data sim/missing-data-pattern.tiff", units="in",
       width=7, height=6, res=300)
  p1
  dev.off()
  
  # Compare event rate in base data to incomplete data
  brate.b <- pyears(Surv(sdat$cmonth, sdat$dead) ~ 1, scale=12)
  brate.cc <- pyears(Surv(sdat.c$cmonth, sdat.c$dead) ~ 1, scale=12)
  
  # Check exposure prevalence in base data and incomplete data
  table(sdat$treat) %>% prop.table()
  table(sdat.c$treat) %>% prop.table()
}

################################################################################
# Estimate associations with outcome and censoring
# One model for the hazard of the outcome event
# One model for the hazard of censoring (this is simply the reverse of the 
# outcome, aka 1-Y)

# Note: the Plasmode package actually selects the second independent variable in 
# the formula for the exposure, rather than the first if you are having the
# package run the model. If you run it yourself it gets it right. 
################################################################################

nsim <- 1000

# sample(1:100000, 1) # 85285
set.seed(85285)

sdat.c <- readRDS("Y:/programs/dan/analyticdat/baseccdat.rds")

os1 <- coxph(Surv(sdat.c$cmonth, sdat.c$dead) ~ treat + genderf + 
               reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
               smokey + surgery + age, 
             data=sdat.c, x=T)

# Censoring hazard
oc1 <- coxph(Surv(sdat.c$cmonth, sdat.c$ndead) ~ treat + genderf +
               reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
               smokey + surgery + age,
             data=sdat.c, x=T)

# treatment probability
# tp1 <- glm(treat ~ genderf +
#              reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#              smokey + 
#              surgery + age, 
#           family=binomial(), 
#           control=glm.control(trace=TRUE),
#           data=sdat.c)
# 
# osf1 <- Surv(sdat.c$cmonth, sdat.c$dead) ~ treat + genderf + 
#   reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#   smokey + surgery + age
# ocf1 <- Surv(sdat.c$cmonth, sdat.c$ndead) ~ treat + genderf +
#   reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#   smokey + surgery + age
# tpf1 <- treat ~ genderf +
#   reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
#   smokey + surgery + age



# Set number of iterations for mega sim
megan <- 100000

# Note on the PlasodeSur package:
   # It seems that - if simulating exposure in addition to outcome - the outcome
   # simulation is not dependent on the exposure simulation. This is contrary to
   # what is stated in the documentation. Possibly this error occurred wile
   # trying to save computation time, as the censoring and outcome survival
   # functions are assigned only once, rather than reassigned for each exposure
   # simulation.

plas.sim <- function(eor, ns, sz){
  PlasmodeSur(
    objectOut = os1,
    objectCen = oc1,
    # objectExp = tp1,
    # formulaOut = osf1,
    # formulaCen = ocf1,
    # formulaExp = tpf1,
    idVar = sdat.c$patientid,
    #idVar = "patientid",
    #data=sdat.c,
    effectOR = eor,
    nsim=ns,
    size=sz
    )
}

# Run simulation #1: effect OR = 1.0
sor10 <- plas.sim(1, nsim, nrow(sdat.c)) # Primary simulation
sor10d <- plas.sim(1, nsim, 2*nrow(sdat.c)) # Double sample size
m.sor10 <- plas.sim(1, 1, megan) # single very large simulation

# Run simulation #2: effect OR = 0.9
sor09 <- plas.sim(0.9, nsim, nrow(sdat.c)) # Primary simulation
sor09d <- plas.sim(0.9, nsim, 2*nrow(sdat.c)) # Double sample size
m.sor09 <- plas.sim(0.9, 1, megan) # single very large simulation

# Run simulation #2: effect OR = 0.9
sor05 <- plas.sim(0.5, nsim, nrow(sdat.c)) # Primary simulation
sor05d <- plas.sim(0.5, nsim, 2*nrow(sdat.c)) # Double sample size
m.sor05 <- plas.sim(0.5, 1, megan) # single very large simulation


# Save these in case we want to rerun the ampute section (next)
sim.dat <- list()
sim.dat[["sor10"]] <- sor10
sim.dat[["sor09"]] <- sor09
sim.dat[["sor05"]] <- sor05
saveRDS(sim.dat, "Y:/programs/dan/plasdat1000.Rdata")

sim.dat.d <- list()
sim.dat.d[["sor10"]] <- sor10d
sim.dat.d[["sor09"]] <- sor09d
sim.dat.d[["sor05"]] <- sor05d
saveRDS(sim.dat.d, "Y:/programs/dan/plasdat1000d.Rdata")

tsim.dat <- list()
tsim.dat[["sor10"]] <- m.sor10
tsim.dat[["sor09"]] <- m.sor09
tsim.dat[["sor05"]] <- m.sor05
saveRDS(tsim.dat, "Y:/programs/dan/plasdatm.Rdata")


### Generate missingness in simulated data #####################################

create.missing <- function(sdat.m, s123, nsim=1000){
  # Create base missing data pattern
  mpat <- tribble(
    ~event,     ~time,     ~treat,  ~genderf, ~reth_black, 
    ~reth_hisp, ~reth_oth, ~age,    ~smokey, 
    ~surgery,  ~ecog12, ~ecog24,
    1, 1, 1, 1, 1,
    1, 1, 1, 1,
    1, 0, 0,
    
    1, 1, 1, 1, 0,
    0, 0, 1, 0,
    1, 0, 0,
    
    1, 1, 1, 1, 0,
    0, 0, 1, 0,
    1, 1, 1,
  )
  
  # Create MAR weights
  # We are using log age for the weights so that age doesn't dominate
  # dichotomous vars.
  mar.w <- tribble(
    ~event,     ~time,     ~treat,  ~genderf, ~reth_black, 
    ~reth_hisp, ~reth_oth, ~age,    ~smokey, 
   ~surgery,    ~ecog12, ~ecog24,
    # Weights for pattern 1
    1, 1, 1, 1, 1,
    1, 1,.1, 1,
    1, 0, 0,
    
    # Weights for pattern 2
    1, 1, 1, 1, 0,
    0, 0,.1, 0,
    1, 0, 0,
    
    # Weights for pattern 3
    1, 1, 1, 1, 0,
    0, 0,.1, 0,
    1, 2, 4
  )
  
  # create MNAR weights, missingness.
  mnar.w <- tribble(
    ~event,     ~time,     ~treat,  ~genderf, ~reth_black, 
    ~reth_hisp, ~reth_oth, ~age,    ~smokey, 
    ~surgery,   ~ecog12, ~ecog24,
    
    # Same weights for all patterns
    1, 1, 1, 1, 1,
    1, 1,.1, 1, 
    1, 2, 4,
    
    1, 1, 1, 1, 0,
    0, 0,.1, 0,
    1, 2, 4,
    
    1, 1, 1, 1, 0,
    0, 0,.1, 0,
    1, 2, 4
  )
  
  # Set pattern frequency
  mfreq <- c(0.8, 0.1, 0.1)
  
  # Set loop values
  mprop.sc <- c(0.1, 0.5, 0.7)
  mmech.sc <- c("MCAR", "MAR", "MNAR")
  mor.lab <- c("or10", "or09", "or05")
  mprop.lab <- c("m10", "m50", "m70")
  
  # Initialize list to contain all simulations
  f.list <- list()
  #f.list[["maindat"]] <- sdat.m %>% select(-b.ecogvalue)
  
  iter <- 0
  
  # Begin loop of missing data treatment effect scenarios
  for(o in 1:length(s123)){
    cur_sim <- s123[[o]]
    # Begin loop of missing data proportion scenarios
    for(i in 1:length(mprop.sc)){
      mprop <- mprop.sc[i]
      
      # Begin loop of missing data mechanism scenarios (3)
      for(m in 1:length(mmech.sc)){
        
        # Set missing data mechanism 
        mmech <- mmech.sc[m]
        mpat <- mpat
        
        # set weights according to mechanism
        if(mmech == "MCAR"){
          freq.w <- NULL
          mmech.w <- NULL
        } else if(mmech == "MAR"){
          freq.w <- mfreq
          mmech.w <- mar.w
        } else{
          freq.w <- mfreq
          mmech.w <- mnar.w
        }
        
        ### Begin loop of nsim ###
        # Note: if this is taking up a lot of memory, can trim redundant variables
        # and re-merge later.
        s.list <- list() # initialize empty storage list
        for(s in 1:nsim){
          sim.index <- s
          # Pull simulation i data
          wdat <- cur_sim[["Sim_Data"]] %>% 
            select(paste("ID", sim.index, sep=""),
                   paste("EVENT", sim.index, sep=""),
                   paste("TIME", sim.index, sep="")
            )
          # paste("EXPOSURE", sim.index, sep=""))
          names(wdat) <- c("patientid", "event", "time")
          
          # Attach full data to simulation i
          id <- wdat$patientid
          wdat <- dplyr::left_join(wdat, sdat.m, by="patientid") %>% select(-patientid) %>%
            mutate(event = ifelse(event == T, 1, 0))
          
          # Add simulation i data to list
          if(i == 1 & m == 1){
            f.list[[mor.lab[o]]][["simdat"]][[paste("sim", sim.index, sep = "")]] <- wdat #%>% 
            #select(-mage)
          }
          
          # Run ampute on full data to simulate missingness in ECOG
          # type set to tail means that extreme weights (positive or negative)
          # indicate a greater probability of missingness.
          wdat.m <- mice::ampute(wdat, prop=mprop, freq=mfreq, patterns=mpat, 
                                 mech=mmech, weights = mmech.w, type="RIGHT")
          
          # Set the missing data pattern deployed (numbers indicate pattern index
          # from mpat. 4 indicates no missing)
          s.list[[paste("sim", sim.index, sep = "")]][["pattern"]] <- 
            ifelse(is.na(wdat.m$amp$ecog12),
                   ifelse(is.na(wdat.m$amp$reth_black), 2, 1), 
                   ifelse(is.na(wdat.m$amp$reth_black), 3, 4))
          s.list[[paste("sim", sim.index, sep = "")]][["iter"]] <- sim.index
          
          iter <- iter + 1
          
          # Print progress
          if(iter %% 50 == 0){
            print(paste("Completed iteration", iter, "of", 
                        length(s123)*length(mprop.sc)*length(mmech.sc)*nsim))
            effs <- survival::coxph(Surv(wdat$time, wdat$event) ~ treat + genderf +
                                      reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                                      smokey  + surgery + age, data=wdat)
            print(paste("SIM OR:", mor.lab[o], "EFF:", exp(effs$coefficients[1])))
          }
          
        } # end nsim loop
        
        f.list[[mor.lab[o]]][[mmech]][[mprop.lab[i]]] <- s.list
      } # end missing data iteration loop
    } # end missing data proportion loop
  } # end treatment effect loop
  
  # Store pattern reference
  mpat.na <- mpat %>% rbind(1)
  mpat.na[mpat.na == 0] <- NA
  f.list[["patterns"]] <- mpat.na # %>% select(-mage)
  
  return(f.list)
}

# Remove unneeded variables from primary data
sdat.c <- readRDS("Y:/programs/dan/analyticdat/baseccdat.rds")
sdat.m <- sdat.c %>% select(-dead, -cmonth, -ndead)

nsim <- 1000

slist <- readRDS("Y:/programs/dan/plasdat1000.Rdata")
missing.sim <- create.missing(sdat.m, slist)
save(missing.sim, file = paste0("Y:/programs/dan/analyticdat/simdat", nsim, ".Rdata"))

slistd <- readRDS("Y:/programs/dan/plasdat1000d.Rdata")
missing.sim.d <- create.missing(sdat.m, slistd)
save(missing.sim.d, file = paste0("Y:/programs/dan/analyticdat/simdat", nsim, "d.Rdata"))

nsim2 <- 1
slist1 <- readRDS("Y:/programs/dan/plasdatm.Rdata")
missing.sim.1 <- create.missing(sdat.m, slist1, nsim2)
save(missing.sim.1, file = paste0("Y:/programs/dan/analyticdat/simdat", nsim2, ".Rdata"))
