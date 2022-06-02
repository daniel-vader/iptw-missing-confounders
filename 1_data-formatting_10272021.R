# Title: Format Flatiron UCC Data
# Notes: Load primary data file(s), format, and output analytic data sets. Adapted
   # from Joanna Harton's data formatting script. Final data set excludes
   # subjects who did not receive therapy.
   # Data sets created include: cleaned_flatiron_ucc_[system date]
   #                            deident_flatiron_ucc_[system date]
   #                            deident_comorb_flatiron_ucc_[system date]
# Author: Daniel Vader
# Date created: 07/15/2020
################################################################################

# Load packages
  # Note: if lubridate is having trouble loading Rcpp.dll uninstall the Rcpp
  # package and install the older version that matches your current R install.
library(dplyr)
library(lubridate)
library(readr)
library(comorbidity)


# CHANGE THE FOLLOWING THREE VARIABLE VALUES WHEN PROCESSING NEW DATA
# Set base path for data files.
  # Not setting WD to avoid accidentally outputting
  # files to the wrong place.
path_base <- "/data/updated data 20211013/edm_bladder_232021/"

# Set date of source data in MMYYYY format. May append this to output file name.
source_date <- "102021"

# Path to output processed data set
path_out <- "/programs/dan/analyticdat/"

# Treatment classifications path
path_trtclass <- "/programs/dan/other/"


################################################################################
### Import data
################################################################################

### Main data ###
demo <- readr::read_csv(
  paste(path_base, "demographics.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

death <- readr::read_csv(
  paste(path_base, "enhanced_mortality_v2.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

biomarker <- readr::read_csv(
  paste(path_base, "enhanced_advurothelialbiomarkers.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

therapy <- readr::read_csv(
  paste(path_base, "lineoftherapy.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

enhanced <- readr::read_csv(
  paste(path_base, "enhanced_advurothelial.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

diagnosis <- readr::read_csv(
  paste(path_base, "diagnosis.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

ecog <- readr::read_csv(
  paste(path_base, "ecog.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

therapy <- readr::read_csv(paste(path_base, "lineoftherapy.csv", sep="")) %>%
  dplyr::rename_all(.funs = tolower)

### Load limited visit and drug episode data
  # (used to determine date of last contact). ###

# Visit has the dates for lab tests, treatment, and vitals
visit <- readr::read_csv(paste(path_base, "visit.csv", sep=""), 
                         col_types = cols_only(PatientID = col_character(), 
                                               VisitDate = col_date())) %>%
  dplyr::rename_all(.funs = tolower) %>%
  dplyr::rename(activitydate = visitdate)

# Drug episode has additional drug administration information. 
  # Since no one in this cohort (not enhanced) has abstracted data, this may not
  # actually be important.
drug <- readr::read_csv(paste(path_base, "drugepisode.csv", sep=""), 
                        col_types = cols_only(PatientID = col_character(), 
                                              EpisodeDate = col_date())) %>%
  dplyr::rename_all(.funs = tolower) %>%
  dplyr::rename(activitydate = episodedate)

################################################################################
# Find the date of last contact for each subject (censoring date if not dead)
################################################################################

# Set the censoring date based on the most recent date from visit and drug
last.date <- rbind(visit, drug) %>% 
  dplyr::group_by(patientid) %>% 
  dplyr::summarize(censordate = max(activitydate))

# Date information from death added later


################################################################################
### Demographics, censoring date, enhanced cohort data. 
### Merge, set death status, set last observation date based on censoring date
### and death date.
################################################################################


# join last.date with demographics and death data
dat0 <- dplyr::left_join(demo, last.date, by="patientid") %>%
  left_join(death, by="patientid") %>% 
  left_join(enhanced, by="patientid")

# Dateofdeath as date
# since death dates are only Month/Year set all dates to 15th of the month as
# per flatiron recommendations
dat0$dateofdeath <- ifelse(!is.na(dat0$dateofdeath) & nchar(dat0$dateofdeath > 4),
                            paste(dat0$dateofdeath, "-15", sep=""), NA) %>%  
  as.Date()

# Create binary indicator for whether subject has observed death
dat0$dead <- ifelse(is.na(dat0$dateofdeath), 0, 1)

# Create date of last observation variable. Considers death and censoring date.
  # Death date is always considered the date of last contact if available.
dat0$datelastcontact <- 
  dplyr::if_else(dat0$dead == 1, dat0$dateofdeath, dat0$censordate)


################################################################################
### Set up variables for first line and second line therapy.
### Only first-line therapy is currently really used in later code 
### (comorbidities). Consider updating if want to fully incorporate second line
### therapy.
################################################################################

# Merge first line and second line therapy data in wide format and exclude
  # any subjects who received no therapy
firsttherapy <- therapy %>% dplyr::filter(linenumber == 1) %>% 
  dplyr::select(patientid, linename, startdate, enddate) %>%
  dplyr::rename(therapy1=linename, startdate1=startdate, enddate1=enddate)

secondtherapy <- therapy %>% dplyr::filter(linenumber == 2) %>%
  dplyr::select(patientid, linename, startdate, enddate) %>%
  dplyr::rename(therapy2=linename, startdate2=startdate, enddate2=enddate)

dat <- 
  dplyr::full_join(dat0, firsttherapy, by="patientid") %>%
  dplyr::full_join(secondtherapy, by = "patientid") %>% 
  filter(!is.na(therapy1))

# In some cases date of last observation will occur prior to the start of therapy.
  # This is because death is only month-specific, and is assigned the 15th day of
  # the month when death occurs. So if a patient dies the same
  # month they start therapy, death may appear to be prior to therapy start.
# temp <- dat %>% select(patientid, startdate1, enddate1, datelastcontact, dateofdeath) %>% 
          # mutate(diff=dateofdeath-datelastcontact) %>% filter(datelastcontact < startdate1)

# Load treatment classifications
trtclass <- readr::read_csv(paste(path_trtclass, 
                                  "treatment_class102621.csv", sep="")) %>%
  mutate(class = tolower(class))

trt.classify <- function(v){
  im <- trtclass %>% dplyr::filter(class == "immuno") %>% dplyr::pull(trt)
  ca <- trtclass %>% dplyr::filter(class == "carbo") %>% dplyr::pull(trt)
  bo <- trtclass %>% dplyr::filter(class == "both") %>% dplyr::pull(trt)
  ne <- trtclass %>% dplyr::filter(class == "neither") %>% dplyr::pull(trt)
  ni <- trtclass %>% dplyr::filter(class == "notimmuno") %>% dplyr::pull(trt)
  c <- ifelse(is.na(v), NA,
              ifelse(v %in% im, "immuno", 
                     ifelse(v %in% ca, "carbo",
                            ifelse(v %in% bo, "both",
                                   ifelse(v %in% ne, "neither", 
                                          ifelse(v %in% ni, "notimmuno", 
                                                 "unknown"))))))
  return(c)
}

dat$trtclass1 <- trt.classify(dat$therapy1)
dat$trtclass2 <- trt.classify(dat$therapy2)

# Export new treatments for classification
# write.csv(table(dat$therapy1[dat$trtclass1 == "unknown"]), paste0(path_out, "2021-10-26_trt-class-table.csv"))

# Output all unclassified therapies. Can be used to help update
   # treatment classificaiton CSV when using new data.
# table(c(dat$therapy1[dat$trtclass1 == "unknown"], 
#         dat$therapy2[dat$trtclass2 == "unknown"])) %>% as.data.frame() %>%
#   readr::write_csv(paste(path_out, "noclass_firstsecond.csv"))


################################################################################
### Classify race/ethnicity. Re-label smoking status.
################################################################################
dat$race.ethnicity <- 
  dplyr::if_else(dat$ethnicity == "Hispanic or Latino" & !is.na(dat$ethnicity), 
                 "Hispanic or Latino",
  dplyr::if_else(is.na(dat$race), "Unknown",
  dplyr::if_else(dat$race == "Hispanic or Latino", "Hispanic or Latino",
  dplyr::if_else(dat$race == "Black or African American", "Black",
  dplyr::if_else(dat$race %in% c("Asian", "Other Race"), "Other", dat$race)))))

# Re-label smoking status variable.
dat$smokingstatus <- plyr::mapvalues(dat$smokingstatus, 
                                     from=c("Unknown/not documented"), 
                                     to=c("Unknown"))


################################################################################
### ECOG: Assign patients baseline ECOG values
################################################################################

ecog2 <- ecog %>% 
  filter(ecog$patientid %in% dat$patientid) %>%
  dplyr::left_join(dat[,c("patientid", "startdate1")], by="patientid")

# How many days prior to the start of treatment was ECOG measured
ecog2$startdiff <- ecog2$startdate1 - ecog2$ecogdate

# Only keep ECOG measurements that occurred at most 7 days after or at most 61
  # days prior to the start of treatment.
  # Then, select the ECOG measurement for each patient that is closest to the 
  # start date. If there is a tie, ECOG = mean(ECOG).
ecog3 <- ecog2 %>% filter(startdiff >= -7, startdiff < 62) %>% 
  group_by(patientid) %>% slice_min(abs(startdiff), with_ties=T) %>%
  dplyr::summarize(b.ecogvalue = mean(ecogvalue), b.ecogdate=min(ecogdate), 
                   b.startdiff=min(startdiff)) 

# Do the same thing but for measurements that occurred 62 days or more prior to
  # treatment, keeping observations closest to 62 days.
ecog4 <- ecog2 %>% filter(startdiff >= 62) %>% 
  group_by(patientid) %>% slice_min(abs(startdiff), with_ties=T) %>%
  dplyr::summarize(o.ecogvalue = mean(ecogvalue), o.ecogdate=min(ecogdate), 
                   o.startdiff=min(startdiff))

# Join ecog data with main data
dat2 <- full_join(dat, ecog3, by="patientid") %>% full_join(ecog4, by = "patientid")

###############################################################################
### De-identify data by replacing dates with relative dates ###
###############################################################################

dat3 <- dat2

# Set baseline year, month, and date to use to create relative times.
dat3$baseline_year <- as.numeric(format(dat3$startdate1, "%Y"))
dat3$baseline_month <- floor_date(dat3$startdate1, unit="months")
dat3$baseline_date <- dat3$startdate1

# Set age
dat3$age <- dat3$baseline_year - dat3$birthyear

# Set relative months. Date of death is actually only month specific so this
# makes more sense.
dat3$startmonth1 <- 
  interval(dat3$baseline_month, floor_date(dat3$startdate1, unit="months")) %/% months(1)
dat3$endmonth1 <- 
  interval(dat3$baseline_month, floor_date(dat3$enddate1, unit="months")) %/% months(1)
dat3$deathmonth <- 
  interval(dat3$baseline_month, floor_date(dat3$dateofdeath, unit="months")) %/% months(1)
dat3$censormonth <- 
  interval(dat3$baseline_month, floor_date(dat3$censordate, unit="months")) %/% months(1)
dat3$lastcontactmonth <- 
  interval(dat3$baseline_month, floor_date(dat3$datelastcontact, unit="months")) %/% months(1)

dat3$diagnosismonth <- 
  interval(dat3$baseline_month, floor_date(dat3$diagnosisdate, unit="months")) %/% months(1)

# Set relative dates
dat3$startdate1 <- as.numeric(dat3$startdate1 - dat3$baseline_date)
dat3$enddate1 <- as.numeric(dat3$enddate1 - dat3$baseline_date)

dat3$startdate2 <- as.numeric(dat3$startdate2 - dat3$baseline_date)
dat3$enddate2 <- as.numeric(dat3$enddate2 - dat3$baseline_date)
dat3$dateofdeath <- as.numeric(dat3$dateofdeath - dat3$baseline_date)
dat3$datelastcontact <- as.numeric(dat3$datelastcontact - dat3$baseline_date)

dat3$diagnosisdate <- as.numeric(dat3$diagnosisdate - dat3$baseline_date) 

dat3$advanceddiagnosisdate <- as.numeric(dat3$advanceddiagnosisdate - dat3$baseline_date)
dat3$b.ecogdate <- as.numeric(dat3$b.ecogdate - dat3$baseline_date)
dat3$o.ecogdate <- as.numeric(dat3$o.ecogdate - dat3$baseline_date)

dat3$surgerydate <- as.numeric(dat3$surgerydate - dat3$baseline_date)

dat4 <- dat3 %>% select(-c(baseline_date, b.startdiff, o.startdiff, baseline_month))


################################################################################
### Add comorbidity info. Code based on ICD codes and calculate Elixhauser score.
################################################################################

# Add start dates for first and second line therapy to full diagnosis data.
   # Drop all subjects who are not receiving therapy
all.diagnosis <- diagnosis %>% dplyr::filter(patientid %in% dat2$patientid)
all.diagnosis <- dat2 %>% dplyr::select(patientid, startdate1) %>% 
  dplyr::left_join(all.diagnosis, by="patientid") %>% 
  filter(!is.na(diagnosisdate))

# Set difference between comorbidity diagnosis and treatment. 
   # If diagnosis date is missing set to NA
all.diagnosis$timediff <- 
  ifelse(is.na(all.diagnosis$diagnosisdate), NA,
                        as.numeric(all.diagnosis$startdate1 - 
                                     all.diagnosis$diagnosisdate))

# Filter for diagnoses with 1 year of lookback
diag1yr <- all.diagnosis %>% dplyr::filter(timediff <= 365, timediff >=0) %>% 
  right_join(dat2[,c("patientid")], by="patientid")

# Fill in a dummy code system for patients with no observed comorbidities
diag1yr$diagnosiscodesystem <- ifelse(is.na(diag1yr$diagnosiscodesystem), 
                                      "ICD-10-CM", diag1yr$diagnosiscodesystem)


# Classify when diagnoses occurred
# diagcat.classifier <- function(v){
#   if (all(is.na(v))) {
#     "Missing start or diagnosis date"
#   } else if (min(v, na.rm=T) >= 0) {
#     "Only diagnoses before/at therapy start"
#   } else if (max(v, na.rm=T) < 0 ) {
#     "Only diagnoses after therapy start"
#   } else {
#     "Diagnoses before and after"
#   }
# }
# all.diagnosis <- all.diagnosis %>% 
#   dplyr::group_by(patientid) %>% 
#   dplyr::summarize(diagnosis.category = diagcat.classifier(timediff)) %>%
#   dplyr::right_join(all.diagnosis, by="patientid")


# Diagnostic code must but in all upper case with no punctuation to work with
   # comorbidity package.
diag1yr$code <- gsub("[[:punct:]]", "", diag1yr$diagnosiscode) %>%
  toupper()


# Assign comorbidity status  and elixhauser score based on ICD-10 and ICD-9 
   # diagnostic code. Use heirarchy of comorbidities
elix10 <- diag1yr %>% dplyr::filter(diagnosiscodesystem == "ICD-10-CM") %>%
  comorbidity(id="patientid", code="code", score="elixhauser", assign0=T,
              icd="icd10") 
elix9 <- diag1yr %>% dplyr::filter(diagnosiscodesystem=="ICD-9-CM") %>%
  comorbidity(id="patientid", code="code", score="elixhauser", assign0=T,
              icd="icd9")

# When a subject has both ICD9 and ICD10 coded conditions, combine their data
   # Set condition indicator to 1 if it shows up in either ICD9 or ICD10 data
   # for patient. After combining, re-calculate scores.
cb910 <- function(v){
  x <- if(sum(v)>0){
    1
  } else{0}
  return(x)
}

elixall <- rbind(elix10, elix9) %>% dplyr::group_by(patientid) %>%
  dplyr::summarize(
    chf = cb910(chf),
    carit = cb910(carit),
    valv = cb910(valv),
    pcd = cb910(pcd),
    pvd = cb910(pvd),
    hypunc = cb910(hypunc),
    hypc = cb910(hypc),
    para = cb910(para),
    ond = cb910(ond),
    cpd = cb910(cpd),
    diabunc = cb910(diabunc),
    diabc = cb910(diabc),
    hypothy = cb910(hypothy),
    rf = cb910(rf),
    ld = cb910(ld),
    pud = cb910(pud),
    aids = cb910(aids),
    lymph = cb910(lymph),
    metacanc = cb910(metacanc),
    solidtum = cb910(solidtum),
    rheumd = cb910(rheumd),
    coag = cb910(coag),
    obes = cb910(obes),
    wloss = cb910(wloss),
    fed = cb910(fed),
    blane = cb910(blane),
    dane = cb910(dane),
    alcohol = cb910(alcohol),
    drug = cb910(drug),
    psycho = cb910(psycho),
    depre = cb910(depre)
  ) %>% mutate(assign0 = T) # Using hierarchy of comorbidity so just set this
                              # vector to TRUE for now for use with the VW algorithm. 

dat5 <- full_join(dat4, elixall, by="patientid", all=TRUE)

elixna <- function(d){
  
}

# Re-calculate score with combined variables 
   # (using code from comorbidity package source)
   # For missing data simulation project, may want to turn this into its
   # own function. Using van Walraven et al. algorithm (2009) 
elixall$wscore_vw <- 
  with(elixall, 
          chf * 7 + 
          carit * 5 + 
          valv * (-1) + 
          pcd * 4 + 
          pvd * 2 + 
          ifelse(hypunc == 1 | hypc == 1, 1, 0) * 0 + 
          para * 7 + 
          ond * 6 + 
          cpd * 3 +
          diabunc * ifelse(diabc == 1 & assign0, 0, 0) + 
          diabc * 0 + 
          hypothy * 0 + 
          rf * 5 + 
          ld * 11 + 
          pud * 0 + 
          aids * 0 + 
          lymph * 9 + 
          metacanc * 12 + 
          solidtum * ifelse(metacanc == 1 & assign0, 0, 4) + 
          rheumd * 0 + 
          coag * 3 + 
          obes * (-4) + 
          wloss * 6 + 
          fed * 5 + 
          blane * (-2) + 
          dane * (-2) +
          alcohol * 0 + 
          drug * (-7) +
          psycho * 0 + 
          depre * (-3)
        )

# Join with full de-identified data 
dat5 <- full_join(dat4, elixall, by="patientid", all=TRUE)

# Subjects with no comorbidities on record get a score of 0

################################################################################
# Write to file
################################################################################

write.csv(dat5, 
          file = paste(path_out, "deident_comorb_flatiron_ucc_",  Sys.Date(),#source_date,
                       ".csv", sep=""), row.names = FALSE)
