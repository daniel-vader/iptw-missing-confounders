# Analysis-related functions:
    # Fit models to simulated data with missingness
    # Fit models to large data set without missingness (used as "truth")
    # Delist model fit statistics and summarize missingness and truth effects in
        # one data frame.
    # Graph distribution of effect estimates by model, missing data proportion,
        # and missing data mechanism.

# Author: Daniel Vader

################################################################################
# Run analysis on different missing data simulations using
# complete case, multiple imputation, or propensity score calibration approaches.
################################################################################

# l = nested data list
# approach = missing data approach (nomiss, cc, mi, psc)
# 
md.analyze <- function(l, approach="cc", misclass=F, diff.miss=F, timer=F, grid=F){
  
  # Store missing patterns reference
  pats <- l$patterns
  
  # Remove missing pattern reference from list
  l$patterns <- NULL
  
  # If running with differential misclassification, generate ecog category
     # probabilities (given death)
  #if(diff.miss){
    eprobs <- ecog.prob(grid)
  #}
  
  # treatment effect
  lapply(seq_along(l), function(li){
    l1 <- l[[li]]
    mp <- names(l)[li]
    # If no missing data (no confounder control) #########################
    if(approach == "nomiss"){
      
      scene <- lapply(l1[[1]], FUN=function(dc){
        # Misclassify outcome and outcome times
        if(misclass){
          mcd <- misclass.outcome(dc$event, dc$time, 
                                  dc$ecog12, dc$ecog24,
                                  eprobs[[mp]], diff.miss)
          dc$event <- mcd[[1]]
          dc$time <- mcd[[2]]
        }
        #print(table(is.na(dc$time)))
        
        
        # Fit model
        m <- survival::coxph(Surv(dc$time, dc$event) ~ treat, data=dc)
        
        coef <- m$coefficients[1]
        se <- sqrt(m$var) %>% as.numeric()
        
        # Calculate 95% CI
        lci <- coef - 1.96*se
        uci <- coef + 1.96*se
        
        # Output treatment coefficient, se, and 95% CI
        fl <- c(coef, se, lci, uci)
        return(fl)
      })
      
      scene.df <- do.call(rbind.data.frame, scene)
      names(scene.df) <- c("coef", "se", "lci", "uci")
      return(scene.df)
      
    # If no missing data (with confounder control) #######################
    } else if(approach == "nomissc"){
      scene <- lapply(l1[[1]], FUN=function(dc){
        # Misclassify outcome and outcome times
        if(misclass){
          mcd <- misclass.outcome(dc$event, dc$time, 
                                  dc$ecog12, dc$ecog24,
                                  eprobs[[mp]], diff.miss)
          dc$event <- mcd[[1]]
          dc$time <- mcd[[2]]
        }
        
        # Fit treatment propensity model
        m.ps <- glm(treat ~ genderf +
                      reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                      smokey + surgery + age, 
                    family=binomial(link="logit"),
                    data=dc)
        dc$ps <- predict.glm(m.ps, type="response")
        
        # Set numerator for stabilization
        num <- mean(dc$treat)
        
        # Calculate IPTWs
        dc$iptw <- ifelse(dc$treat == 1,
                          num/dc$ps, # weights for treated
                          (1-num)/(1-dc$ps)) # weights for untreated
        
        # Fit propensity score model
        m <- survival::coxph(Surv(dc$time, dc$event) ~ treat, data=dc, weights = dc$iptw)
        
        coef <- m$coefficients[1]
        se <- sqrt(m$var) %>% as.numeric()
        
        # Calculate 95% CI
        lci <- coef - 1.96*se
        uci <- coef + 1.96*se
        
        # Output treatment coefficient, se, and 95% CI
        fl <- c(coef, se, lci, uci)
        return(fl)
      })
      
      scene.df <- do.call(rbind.data.frame, scene)
      names(scene.df) <- c("coef", "se", "lci", "uci")
      return(scene.df)
      
      
      # Else run the missing data methods ######################################
    } else {
    
      #mechanism
      lapply(l1[-1], FUN = function(l2) {
        
        # proportion
        pblapply(l2, FUN = function(l3){
          # sim
          scene <- lapply(l3, FUN = function(d){
            # Grab data for current sim
            dc <- l1[[1]][[d[["iter"]]]]
            
            # Misclassify outcome and outcome times
            if(misclass){
              mcd <- misclass.outcome(dc$event, dc$time, 
                                      dc$ecog12, dc$ecog24,
                                      eprobs[[mp]], diff.miss)
              dc$event <- mcd[[1]]
              dc$time <- mcd[[2]]
            }
            
            # Remove data by missing data pattern
            dc$cpat <- d[["pattern"]]
            
            dc$rindex <- as.numeric(rownames(dc))
            dc2 <- NULL
            for(z in 1:nrow(pats)){
              cc <- dc %>% filter(cpat == z) %>%
                sweep(MARGIN=2, c(as.numeric(pats[.$cpat[z],]), 1, 1),'*')
              dc2 <- rbind(dc2, cc)
            }
            dc <- dc2 %>% arrange(rindex) %>% select(-cpat, -rindex) 
            
            if(timer){st <- Sys.time()}
            
            # If complete case analysis ##########################################
            if(approach == "cc"){
              
              # drop all rows with missing data
              dc.f <- dc %>% na.omit()
              
              # Fit treatment propensity model
              m.ps <- glm(treat ~   genderf +
                            reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                            smokey + surgery + age, 
                          family=binomial(link="logit"),
                          data=dc.f)
              dc.f$ps <- predict.glm(m.ps, type="response")
              
              # Set numerator for stabilization
              num <- mean(dc.f$treat)
              
              # Calculate IPTWs
              dc.f$iptw <- ifelse(dc.f$treat == 1,
                                num/dc.f$ps, # weights for treated
                                (1-num)/(1-dc.f$ps)) # weights for untreated
              
              # Fit propensity score model
              old.warn <- options(warn=2)
              bcatch <- try(
                m <- survival::coxph(Surv(dc.f$time, dc.f$event) ~ treat, 
                                     data=dc.f, weights = dc.f$iptw))
              if(class(bcatch) != "try-error"){
                coef <- m$coefficients[1]
                se <- sqrt(m$var) %>% as.numeric()
              }
              
              
              # If multiple imputation ###########################################
            } else if(approach == "mi") {
              # Attach the subject-specific Nelson-Aaalen estimator of the 
              # cumulative hazard rate for use in imputation models
              # White & Royston (2009) doi: 10.1002/sim.3618 
              dc$htna <- nelsonaalen(as.data.frame(dc), timevar=time, statusvar=event) 
              dc$ecog <- ifelse(dc$ecog12 == 0 & dc$ecog24 == 0, 1,
                                ifelse(dc$ecog12 == 1, 2, 3)) %>% 
                factor(levels = c(1,2,3), labels=c("PS < 1", "1 <= PS < 2", "2 <= PS"))
              dc$reth <- ifelse(dc$reth_black == 0 & dc$reth_hisp == 0 & dc$reth_oth == 0, 1,
                                ifelse(dc$reth_black == 1, 2,
                                       ifelse(dc$reth_hisp == 1, 3, 4))) %>%
                factor(levels = c(1,2,3,4), labels=c("White", "Black", "Hispanic", "Other"))
              
              # remove indicator ecog vars, time  (mi will use all available vars)
              dc2 <- dc %>% select(-ecog12, -ecog24, -time, -reth_black,
                                   -reth_hisp, -reth_oth) 
              
              dc2 <- dc2 %>% mutate(genderf=as.factor(genderf), 
                                    smokey=as.factor(smokey),
                                    surgery=as.factor(surgery))
                 
              
              # Imputation methods
              meth <- dc2[1,] %>% 
                mutate(event="", treat="", genderf="logreg", age="", smokey="logreg",
                       surgery="logreg", htna="", ecog="polr",
                       reth="polyreg") #%>% as.character()
              
              # Create imputed data sets, use predictive mean matching to estimate
              # ECOG
              impd <- mice.par(dc2, m=10, maxit=5, nnodes=5, method=meth, print=F)
              
              # Fit treatment propensity model
              m.ps <- with(impd, 
                           glm(treat ~  genderf +
                                 reth + 
                                 ecog + smokey + surgery + age, 
                               
                          family=binomial(link="logit")))
              
              
              # Setup vectors for storage
              coef.v <- NULL
              var.v <- NULL
              
              # Set numerator for stabilization
              num <- mean(impd$data$treat)
              
              for(ms in 1:length(m.ps$analyses)){
                pst <- predict.glm(m.ps$analyses[[ms]], type="response")
                
                # Calculate IPTWs
                iptw <- ifelse(impd$data$treat == 1,
                                    num/pst, # weights for treated
                                    (1-num)/(1-pst)) # weights for untreated 
                
                # Fit model
                m <- survival::coxph(Surv(dc$time, dc$event) ~ dc$treat,
                                weights = iptw)
                
                # Pull coefficient and variance
                coef.v <- c(coef.v, m$coefficients[1])
                var.v <- c(var.v, as.numeric(m$var))
              }
              
              # Pool
              md.pooled <- pool.scalar(coef.v, var.v, n=nrow(dc))
              
              # pull pooled treatment coefficient and se
              coef <- md.pooled$qbar
              se <- sqrt(md.pooled$ubar)
              
              # create dummy bcatch for later if statement
              bcatch <- "n"
  
                      
              # If propensity score calibration ##################################
            } else if(approach == "psc"){
              
              # Calculate error prone propensity score (model with all variables 
                 # containing complete data)
              m.ep <- glm(treat ~ genderf + surgery + age, 
                          family=binomial(link="logit"),
                          data=dc)
              dc$ps.ep <- predict.glm(m.ep, type="response")# ep probabilities
              dc$l.ps.ep <- log(dc$ps.ep/(1-dc$ps.ep)) # get ep log odds
              
              # Catch sims with rank deficent fit
              if(length(m.ep$coefficients) > m.ep$rank){
                print("Error prone PS: Rank deficient fit")
                print(summary(m.ep))
              }
              
              # Calculate gold standard ps (model with all variables and rows containing complete data)
              dc.f <- na.omit(dc)
              m.gs <- glm(treat ~  genderf + 
                            reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                            smokey + surgery + age, 
                          family=binomial(link="logit"),
                          data=dc.f)
              dc.f$ps.gs <- predict.glm(m.gs, type="response") # gs probabilities
              dc.f$l.ps.gs <- log(dc.f$ps.gs/(1-dc.f$ps.gs)) # gs log odds
              
              # Catch sims with rank deficent fit
              if(length(m.gs$coefficients) > m.gs$rank){
                print("Gold standard PS: Rank deficient fit")
                print(summary(m.gs))
              }
              
              
              
              # Regress the gold standard ps on the error prone ps
              # Some reference papers have this as a linear model but log linear
              # makes more sense given all values are between 0 and 1.
              # TODO: take log odds of ps.gs and ps.ep and then run a linear model. 
                  # Then transform back to probability scale.
              m.epgs <- lm(l.ps.gs ~ treat + l.ps.ep, data=dc.f)
              
              # Catch sims with rank deficent fit
              if(length(m.epgs$coefficients) > m.epgs$rank){
                print("Error prone GS: Rank deficient fit")
                print(summary(m.epgs))
              }
              
              
              # Assign gold standard ps to subjects with complete data and gs regressed
              # ep ps to subjects with incomplete data.
              epgs.odds <- exp(predict(m.epgs, dc))
              dc$epgs.ps <- epgs.odds/(1+epgs.odds)
              
              dc$gs.ps <- predict.glm(m.gs, dc, type="response")
              
              dc$ps <- ifelse(is.na(dc$gs.ps), dc$epgs.ps, dc$gs.ps)
              
              # Estimate stabilized IPTWs (ATE)
              num <- mean(dc$treat) # P(treat = 1)
              if(min(dc$ps) < 0){
                dc.min <- dc[which.min(dc$ps),]
                print(paste("m.ps=",dc.min$ps, "m.epgs=", dc.min$epgs.ps,
                            "m.gs=", dc.min$gs.ps, "iter"=d[["iter"]]))
              }
              dc$iptw <- ifelse(dc$treat == 1,
                                num/dc$ps, # weights for treated
                                (1-num)/(1-dc$ps)) # weights for untreated
              
              # Fit model
              old.warn <- options(warn=2)
              bcatch <- try(
                  m <- survival::coxph(Surv(dc$time, dc$event) ~ dc$treat,
                                       data=dc, weights = dc$iptw))
              if(class(bcatch) != "try-error"){
                  coef <- m$coefficients
                  se <- sqrt(m$var)
              }
            } 
            
            # Drop simulation if Cox model failed
            if(class(bcatch) != "try-error"){
              # Calculate 95% CI
              lci <- coef - 1.96*se
              uci <- coef + 1.96*se
            
              # Output treatment coefficient, se, and 95% CI
              fl <- c(coef, se, lci, uci)
            } else{
              print(paste("Failed on data simulation: ", d[["iter"]]))
              fl <- c(NA, NA, NA, NA)
            }
            if(!timer){
              return(fl)
            } else{
              et <- Sys.time()
              elapsed <- et-st
              return(elapsed)
            }
          }) # end scene
          if(!timer){
            scene.df <- do.call(rbind.data.frame, scene)
            names(scene.df) <- c("coef", "se", "lci", "uci")
            return(scene.df)
          }else{
            return(scene)
          }
        }) # end proportion
      }) # End mechanism
    } # End if else no missing
  }) # End treatment effect
}

################################################################################
# Misclassify outcomes
# Takes dichotomous outcome vector o, continuous time vector t
# Logical diff determines whether misclassification is differential
# b and c indicates variable categories that missingness is differential upon
################################################################################
misclass.outcome <- function(o, t, b, c,  prps, diff=F,
                             scale.b = 0.9, # sens in b is 0.9 * sens in group a
                             scale.c = 0.8 # sens in c is 0.8 * sens in group a
                             ){
  sens <- .902
  spec <- .966
  n.rand <- runif(length(o), 0, 1)
  if(!diff){
    n.o <- ifelse(o == 1, 
                  # If truly died, check against sensitivity
                  ifelse(n.rand > sens, 0, 1),
                  # If truly alive, check against specificity
                  ifelse(n.rand > spec, 1, 0))
  }else{
    pa <- prps[1]
    pb <- prps[2]
    pc <- prps[3]
    
    sa <- (sens/(pa+pb*scale.b+pc*scale.c))
    sb <- scale.b*sa
    sc <- scale.c*sa
    
    # If this scales a number greater than or equal to 1, scale back to a more
       # realistic probability
    sa <- if(sa >= 1){.999
      }else(sa)
    sb <- if(sb >= 1){.999
      }else(sb)
    sc <- if(sc >= 1){.999
      }else(sc)
    
    n.o <- ifelse(o == 1, 
                  # If truly died, check against sensitivity
                  ifelse(b == 1, 
                         # If group b, use group b sensitivity
                         ifelse(n.rand > sb, 0, 1),
                         # If not group a but if group c, use group c sensitivity
                         ifelse(c == 1, ifelse(n.rand > sc, 0, 1),
                                # Else must be group a, use group a sensitivity
                                ifelse(n.rand > sa, 0, 1))),
                  # If truly alive, check against specificity
                  ifelse(n.rand > spec, 1, 0))
  }
  
  # Alter times if misclassified
  n.rand2 <- runif(length(o), 0, 1)
  n.t <- ifelse(n.o < o, # If false negative
                ifelse(t <= 4, # If patient died at or before 3 months
                       ifelse(t == 1, t, # If died at month 1, no change
                              ifelse(n.rand2 <= .5, t, t-1)), # 50% roll back 1 month
                       ifelse(n.rand2 <= .30, t,
                              ifelse(n.rand2 <= .70, t-1,
                                     ifelse(n.rand2 <= .9, t-2, t-3)))), 
                t) # else keep same
  
  return(list(n.o, n.t))
}

ecog.prob <- function(grid){
  if(grid){
    md <- readRDS("../analyticdat/baseccdat.rds")
    bd <- readRDS("../plasdatm.Rdata")
  }else{
    md <- readRDS("/programs/dan/analyticdat/baseccdat.rds")
    bd <- readRDS("/programs/dan/plasdatm.Rdata")
  }
  
  bd.names <- names(bd) %>% sub('.', '', .)
  plist <- list()
  for(i in 1:length(bd)){
    cd <- bd[[i]]
    fd <- tibble(patientid = cd$Sim_Data[, "ID1"], 
                 event = cd$Sim_Data[, "EVENT1"], 
                 time = cd$Sim_Data[, "TIME1"],
                 #treat = sim$Sim_Data[, "EXPOSURE1"]
    ) %>%
      left_join(md, by="patientid") %>%
      filter(event == 1)
    
    p.ecog12 <- sum(fd$ecog12)/nrow(fd)
    p.ecog24 <- sum(fd$ecog24)/nrow(fd)
    p.ecog01 <- 1 - p.ecog12 - p.ecog24
    plist[[bd.names[i]]] <- (c(p.ecog12, p.ecog24, p.ecog01))
  }
  return(plist)
}

################################################################################
# Simulate "truth" from large data set
################################################################################
tsim <- function(sim, main, boosted = F){
  
  # Merge sim data with main data ############################################
  dc <- tibble(patientid = sim$Sim_Data[, "ID1"], 
               event = sim$Sim_Data[, "EVENT1"], 
               time = sim$Sim_Data[, "TIME1"],
               #treat = sim$Sim_Data[, "EXPOSURE1"]
               ) %>%
    left_join(main, by="patientid")
  
  # fit ps model #############################################################
  if(!boosted){ # PS estimated with logistic regression
    # Calculate ps
    m.ps <- glm(treat ~ genderf +
                  reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                  smokey + surgery + age, 
                family=binomial(link="logit"),
                data=dc)
    dc$ps <- predict.glm(m.ps, type="response")
    
  } else{ # PS estimated with boosted CART
    dc2 <- dc %>% 
      select(-patientid, -event, -time, -dead, -ndead, -cmonth) %>%
      as.data.frame()
    
    # Find ATE propensity scores using twang (boosted CART)
    pscore.cart <- ps(treat ~ genderf + reth_black + reth_hisp + reth_oth + ecog12 + 
                        ecog24 + smokey + surgery + age,
                      n.trees = 10000,
                      data=dc2)
    
    dc$ps <- pscore.cart$ps$ks.mean.ATE
    m.ps <- NA
  }
  
  # Calculate stabilized weights
  num <- mean(dc$treat)
  dc$iptw <- ifelse(dc$treat == 1,
                    num/dc$ps, # weights for treated
                    (1-num)/(1-dc$ps)) # weights for untreated
  
  # Fit model
  m2 <- survival::coxph(Surv(time, event) ~ treat, data=dc, weights = dc$iptw)
  
  coef <- m2$coefficients[1]
  se <- sqrt(m2$var) %>% as.numeric()
  
  # Calculate 95% CI
  lci <- coef - 1.96*se
  uci <- coef + 1.96*se
  
  # Output treatment coefficient, se, and 95% CI
  fl <- tribble(~coef, ~se, ~lci, ~uci,
                coef, se, lci, uci)
  out <- list()
  out[["results"]] <- fl
  out[["ps_model"]] <- m.ps
  out[["cox_model"]] <- m2
  out[["iptw"]] <- dc$iptw
  return(out)
}

# Old version of tsim 
tsim.old <- function(sim, main, n, coxfit=F){
  fl <- tibble(coef = numeric(),
               se = numeric(),
               lci = numeric(),
               uci = numeric(),
               ps = character())
  
  for(i in 1:n){
    # Merge sim data with main data ############################################
    dc <- tibble(patientid = sim$Sim_Data[, paste0("ID", i)], 
                event = sim$Sim_Data[, paste0("EVENT", i)], 
                time = sim$Sim_Data[, paste0("TIME", i)]) %>%
      left_join(main, by="patientid")
    
    # fit cox model ############################################################
    if(coxfit){
      m <- survival::coxph(Surv(time, event) ~ treat +  genderf + 
                             reth_black + reth_hisp + reth_oth + 
                             ecog12 + ecog24 + smokey + surgery +
                             age, data=dc)
      coef <- m$coef["treat"] %>% as.numeric()
      se <- sqrt(m$var[1,1])
      
      # Calculate 95% CI
      lci <- coef - 1.96*se
      uci <- coef + 1.96*se
      
      # Output treatment coefficient, se, and 95% CI
      fl <- fl %>% add_row(coef=coef, se=se, lci=lci, uci=uci, ps="no")
    }

    # fit ps model #############################################################
    
    # Calculate ps
    m.gs <- glm(treat ~   genderf +
                  reth_black + reth_hisp + reth_oth + ecog12 + ecog24 + 
                  smokey + surgery + age, 
                family=binomial(link="logit"),
                data=dc)
    dc$ps <- predict.glm(m.gs, type="response")
    
    num <- mean(dc$treat)
    
    dc$iptw <- ifelse(dc$treat == 1,
                      num/dc$ps, # weights for treated
                      (1-num)/(1-dc$ps)) # weights for untreated
    
    # Fit model
    m2 <- survival::coxph(Surv(time, event) ~ treat, data=dc, weights = dc$iptw)
    
    coef <- m2$coefficients[1]
    se <- sqrt(m2$var) %>% as.numeric()
    
    # Calculate 95% CI
    lci <- coef - 1.96*se
    uci <- coef + 1.96*se
    
    # Output treatment coefficient, se, and 95% CI
    fl <- fl %>% add_row(coef=coef, se=se, lci=lci, uci=uci, ps="yes")
    
    # Trakcer
    print(paste("Iteration", i, "of", n))
  }
  return(fl)
}

################################################################################
# Change lists of effects into dataframe (tibble) format. Calculate mean "truth"
# effects.
################################################################################

delist_eff <- function(dat, tp, psv="no", miss=T, mset){
  simhr <- c("Effect = ln(1)", "Effect = ln(0.9)", "Effect = ln(0.5)")
  d <- tibble(mech=character(), prop=numeric(), est=numeric(), est.hr=numeric(),
              est.true=numeric(), est.hr.true=numeric(), cover=numeric())
  
  tp10 <- tp$sor10$results %>% pull(coef) 
  tp09 <- tp$sor09$results %>% pull(coef) 
  tp05 <- tp$sor05$results %>% pull(coef) 
  
  tp2 <- c(tp10, tp09, tp05)
  
  if(miss){
    mechs <- names(dat[[1]])
    #props <- c(.1,.3,.5, .7)
    props <- c(.1, .5, .7)
    
    d <- tibble(mech=character(), prop=numeric(), est=numeric(), est.hr=numeric(),
                est.true=numeric(), est.hr.true=numeric(), bias.hr=numeric(),
                cover=numeric())
    for(i in 1:length(dat)){
      for(j in 1:length(dat[[1]])){
        for(k in 1:length(dat[[1]][[1]])){
          td <- dat[[i]][[j]][[k]] %>%
            mutate(mech=mechs[j], prop=props[k], est.true=tp2[i], bias=coef-tp2[i],
                   est.hr.true=exp(est.true), est.hr=exp(coef), sim.hr=simhr[i],
                   bias.hr = est.hr - est.hr.true,
                   cover=ifelse(est.true < lci | est.true > uci, 0, 1))
          d <- rbind(d, td)
        }
      }
    }
    
    return(d)
  } else {
    for(i in 1:length(dat)){
      d1 <- dat[i] %>% as.data.frame()
      colnames(d1) <- c("coef", "se", "lci", "uci")
      td <- d1 %>% 
        mutate(mech="MCAR", prop=0, est.true=tp2[i], bias=coef-tp2[i], 
                est.hr.true=exp(est.true), est.hr=exp(coef), sim.hr=simhr[i],
                bias.hr = est.hr - est.hr.true,
                cover=ifelse(est.true < lci | est.true > uci, 0, 1))
      d <- rbind(d, td)
    }
    return(d)
  }
}

delist_dur <- function(dat){
  method <- c("CC", "MI", "PSC")
  mechs <- names(dat[[1]][[1]])
  props <- c(.1, .5, .7)
  simhr <- c("Effect = ln(1)", "Effect = ln(0.9)", "Effect = ln(0.5)")
  d <- tibble(approach=character(), mech=character(), prop=numeric(), sim.hr=character(), dur=numeric())

  for(i in 1:length(dat)){
    for(j in 1:length(dat[[1]])){
      for(k in 1:length(dat[[1]][[1]])){
        for(l in 1:length(dat[[1]][[1]][[1]])){
          td <- dat[[i]][[j]][[k]][[l]][[1]] %>% as.double(units="secs") %>%
            round(digits=0)
          d <- d %>% add_row(approach=method[i], mech=mechs[k], prop=props[l], sim.hr=simhr[j],
                       dur=td)
        }
      }
    }
  }
  
  return(d)
}

# Delist and aggregate all approaches at once.
delist_all <- function(cc, psc, mi, cd, cdc, tr){
  
  effs <- list()
  eff.cc.d <- delist_eff(cc, tr, ps="yes") %>% 
    mutate(approach="CC")
  eff.psc.d <- psc %>% delist_eff(tr, ps="yes") %>% 
    mutate(approach="PSC")
  eff.mi.d <- mi %>% delist_eff(tr, ps="yes") %>%
    mutate(approach="MI")
  
  eff.cd.d <- delist_eff(cd, tr, ps="yes", miss=F, mset="Univ") %>% 
    mutate(approach="Univ")
  eff.cdc.d <- delist_eff(cdc, tr, ps="yes", miss=F, mset="Multiv") %>% 
    mutate(approach="Multiv")
  
  effs[["eff.all"]] <- rbind(eff.cc.d, eff.psc.d, eff.mi.d) %>% 
    mutate(mech = factor(mech, levels=c("MCAR", "MAR", "MNAR")))
  
  effs[["eff.all.cd"]] <- rbind(eff.cc.d, eff.psc.d, eff.mi.d, eff.cd.d, eff.cdc.d) %>% 
    mutate(mech = factor(mech, levels=c("MCAR", "MAR", "MNAR")))
  
  return(effs)
}

################################################################################
# Plot effects and bias at simulated HR (thr) by missingness mechanism 
# and proportion missing.
################################################################################
plot.effects <- function(d){
  ggplot(d, aes(x=mech, y=est.hr, fill=factor(prop))) +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", color="black", size=1, 
                 shape=15, position=position_dodge(.75)) +
    facet_grid(sim.hr ~ approach) +
    geom_hline(aes(yintercept = est.hr.true)) +
    scale_fill_viridis(discrete=T) +
    ylab("Hazard Ratio") +
    labs(fill="Proportion \n Missing") +
    theme(axis.title.x = element_blank())
}

plot.bias <- function(d, lim=T, alt=F){
  if(alt){
    d <- d %>% 
      mutate(prop = paste0(prop*100, "%"),
             sim.hr = paste0("True HR = ",round(est.hr.true, 2)),
             approach = ifelse(approach == "CC", "Complete case",
                               ifelse(approach == "MI", "Multiple imputation",
                                      "Propensity score\ncalibration"))) %>%
      group_by(sim.hr, prop, mech, approach) %>% 
      summarize(y025=quantile(bias.hr, 0.025, na.rm=T),
                y25=quantile(bias.hr, 0.25, na.rm=T),
                y50=quantile(bias.hr, 0.50, na.rm=T),
                y75=quantile(bias.hr, 0.75, na.rm=T),
                y975=quantile(bias.hr, 0.975, na.rm=T),
                ymean=mean(bias.hr, na.rm=T))
    pl <- ggplot(d, aes(x=prop)) +
      geom_boxplot(aes(ymin=y025, lower=y25, middle=y50, upper=y75, ymax=y975, 
                       fill=factor(approach)),
                   stat="identity",
                   outlier.shape=NA) +
      geom_point(aes(y=ymean, shape="Mean", group=factor(approach)), 
                 position=position_dodge(.9), size=2) +
      geom_point(aes(y=y50, shape="Median", group=factor(approach)), 
                 position=position_dodge(.9), alpha=0, size=3) +
      # stat_summary(fun=mean, aes(shape="Mean"), geom="point", color="black", size=1, 
      #              #shape=15, 
      #              position=position_dodge(.75), show.legend = T) +
      # stat_summary(fun = median, aes(shape="Median"), geom = "point", color="black", size=1, 
      #   #shape=95, 
        # position=position_dodge(.75), alpha=0, show.legend = T) +
      scale_shape_manual("Measures of\ncentral tendency", values=c(20, 95)) +
      facet_grid(sim.hr ~ mech) +
      geom_hline(aes(yintercept = 0)) +
      scale_fill_viridis(discrete=T, begin=0.3, end=1)  +
      #scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
      ylab("Bias") +
      xlab("% missing data") +
      labs(fill="Missing data\napproach") +
      guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) +
      theme_bw()
    
    if(lim){
      pl <- pl + 
        coord_cartesian(ylim=c(-0.45,0.3)) +
        scale_y_continuous(breaks = c(-0.45, -.30, -.15, 0, 0.15, 0.3, 0.45))
      
    } else {
      pl <- pl + coord_cartesian(ylim=c(-1,4))
    }
  } else {
    pl <- ggplot(d, aes(x=mech, y=bias, fill=factor(prop))) +
      geom_boxplot() +
      stat_summary(fun=mean, geom="point", color="black", size=1, 
                   shape=15, position=position_dodge(.75)) +
      facet_grid(sim.hr ~ approach) +
      geom_hline(aes(yintercept = 0)) +
      scale_fill_viridis(discrete=T, begin=0.1, end=1) +
      #scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
      ylab("Absolute bias") +
      labs(fill="Proportion \n Missing") +
      theme_bw() + 
      theme(axis.title.x = element_blank())
    if(lim){
      pl <- pl + ylim(-1.1,0.6)
    }
  }
  
  return(pl)
}

# Tabulate bias quantification statistics on the HR scale
tab.bias <- function(effs){
  coverage.errors <- effs %>% 
    # mutate(sim.hr = ifelse(sim.hr == "Simulated HR = 0.5", "hr05",
    #                        ifelse(sim.hr == "Simulated HR = 0.9", "hr09",
    #                               ifelse(sim.hr == "Simulated HR = 1", "hr10", "NA")))
    #        ) %>%
    mutate(sim.hr = round(est.hr.true, 2), 
           prop = paste0(prop*100, "%")) %>%
    group_by(sim.hr, prop, mech, approach) %>% 
    summarize(bias=mean(bias.hr, na.rm=T), mse=mean(bias.hr^2, na.rm=T),
              coverage95 = mean(cover, na.rm=T)) %>%
    tidyr::pivot_wider(names_from=c(mech), 
                       values_from=c(bias, mse, coverage95)) %>%
    mutate(sortb = case_when(approach == "Univ" ~ 1,
                             approach == "Multiv" ~ 2,
                             TRUE ~ 3)) %>%
    arrange(sim.hr, prop, sortb, approach) %>%
    relocate(sim.hr, prop, approach, 
             bias_MCAR, mse_MCAR, coverage95_MCAR,
             bias_MAR, mse_MAR, coverage95_MAR,
             bias_MNAR, mse_MNAR, coverage95_MNAR
             
             # bias_MCAR_hr05, mse_MCAR_hr05, coverage95_MCAR_hr05,
             # bias_MAR_hr05, mse_MAR_hr05, coverage95_MAR_hr05,
             # bias_MNAR_hr05, mse_MNAR_hr05, coverage95_MNAR_hr05,
             # 
             # bias_MCAR_hr09, mse_MCAR_hr09, coverage95_MCAR_hr09,
             # bias_MAR_hr09, mse_MAR_hr09, coverage95_MAR_hr09,
             # bias_MNAR_hr09, mse_MNAR_hr09, coverage95_MNAR_hr09,
             # 
             # bias_MCAR_hr10, mse_MCAR_hr10, coverage95_MCAR_hr10,
             # bias_MAR_hr10, mse_MAR_hr10, coverage95_MAR_hr10,
             # bias_MNAR_hr10, mse_MNAR_hr10, coverage95_MNAR_hr10
             )
}

tab.eff <- function(effs){
  effect_summary <- effs %>% 
    mutate(sim.hr = paste0("log(", round(est.hr.true, 2), ")"),
           prop = paste0(prop*100, "%")) %>%
    group_by(sim.hr, prop, mech, approach) %>% 
    summarize(#eff=mean(coef, na.rm=T),
              empirical_error=sqrt(var(coef, na.rm=T)), 
              model_error=mean(se, na.rm=T)) %>%
    tidyr::pivot_wider(names_from=mech, 
                       values_from=c(#eff, 
                                     empirical_error, model_error)) %>%
    mutate(sortb = case_when(approach == "Univ" ~ 1,
                             approach == "Multiv" ~ 2,
                             TRUE ~ 3)) %>%
    arrange(sim.hr, prop, sortb, approach) %>%
    relocate(sim.hr, prop, approach, 
             #eff_MCAR, 
             empirical_error_MCAR, model_error_MCAR,
             #eff_MAR, 
             empirical_error_MAR, model_error_MAR,
             #eff_MNAR, 
             empirical_error_MNAR, model_error_MNAR)
  return(effect_summary)
}

# Plot 95% CI Coverage probabilities and MSE
plot.mse95 <- function(dat, lim1=T, lim2=T){
  eb95 <- dat$eff.all %>% 
    mutate(prop = paste0(prop*100, "%"),
           sim.hr = paste0("True HR = ",round(est.hr.true, 2)),
           approach = ifelse(approach == "CC", "Complete case",
                             ifelse(approach == "MI", "Multiple imputation",
                                    "Propensity score\ncalibration"))) %>% 
    group_by(sim.hr, prop, mech, approach) %>% 
    summarize(cover95 = mean(cover, na.rm=T),
              MSE = mean(bias.hr^2, na.rm=T))
  
  pl1 <- ggplot(eb95, aes(x=factor(prop), y=cover95, fill=factor(approach))) +
    geom_col(position="dodge", width = 0.8) +
    geom_hline(aes(yintercept = .95)) +
    facet_grid(sim.hr ~ factor(mech)) +
    scale_fill_viridis(discrete=T, begin=.3, end=1) +
    #scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
    ylab("95% CI coverage probability") +
    xlab("% Missing data") +
    labs(fill="Missing data\nApproach") +
    theme_bw()
  
  pl2 <- ggplot(eb95, aes(x=factor(prop), y=MSE, fill=factor(approach))) +
    geom_col(position="dodge", width=0.8) +
    facet_grid(sim.hr ~ factor(mech)) +
    #ylim(0.7,1) +
    scale_fill_viridis(discrete=T, begin=.3, end=1) +
    #scale_y_break(c(0, .02, .04, .08, .16)) +
    ylab("Mean squared error") +
    xlab("% Missing data") +
    labs(fill="Missing data\nApproach") +
    theme_bw()
  
  if(lim1){
    pl1 <- pl1 + coord_cartesian(ylim=c(0.45, 1.0)) +
      scale_y_continuous(breaks = c(0.5, .6, .7, .8, .9, 1))
  }
  if(lim2){
    pl2 <- pl2 + coord_cartesian(ylim=c(0, .05))
  }
  
  return(list(cover95 = pl1, mse = pl2))
}

# Let's load RDATA files with different names
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}