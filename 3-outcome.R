## ---------------------------
##
## Program: 3. Outcomes and follow-up
##
## Purpose: Flag the occurrence and date of first occurrence of the outcome of interest since cohort entry.
## Here, the outcome is all-cause mortality, where date of death (dod) is the earliest among:
## Emis dod, CPRD dod, HES dod, and ONS dod. 
## 
## Define follow-up for the two following analyses types: 
## - intention-to-treat analysis (ITT): patients are analyzed according to their original 
## exposure group. Follow-up ends at the earliest among: death ('dod'), occurrence of the event 
## of interest ('event_date'), departure from the CPRD ('regend'), end of linkage to HES or ONS, 
## last available data ('lcd'), or study follow-up end ('study_follow_up_end')
## - as-treated: patients are analyzed according to their actual exposure group, and are thus
## censored when they switch or discontinue their treatment. Follow-up ends at the earliest among
## the above plus treatment switch ('switch_date') or treatment discontinuation ('disc_date'). 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes: For this project, both CPRD and HES are used to identify outcomes. The date of first occurrence
## of the outcome is first identified in CPRD and HES separately, and then the earliest between both
## is used to define first outcome occurrence date. 
##
## ---------------------------

# analysis: main, flex_grace_period, or 90_day_grace_period

analysis <- ''

#### LOAD PACKAGES ####

library(lubridate)
library(readxl)
library(dplyr)
library(magrittr)
library(haven)
library(parallel)
library(data.table)

#### SPECIFY STUDY DESIGN ####

study_start = ymd(20190101)
study_end = ymd(20221231)
study_follow_up_end = ymd(20240331)

#### DEFINE PATHS ####

if (analysis == 'main' | analysis == '') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
} else if (analysis == 'flex_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/flex_grace_period" 
} else if (analysis == '90_day_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period" 
} 

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
path_linkage_1 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt1'
path_linkage_2 <- 'Z:/EPI/Protocol 24_004042/Data linkage/Results/Aurum_linked/Final_pt2'

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_censored.rds', sep = '/'))
switched_to <- readRDS(file = paste(path_main, 'switched_to.rds', sep = '/'))

setwd(path_cohort)

outcome_desc <- "outcome_desc.txt"
writeLines("Outcome description:", outcome_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = outcome_desc, append= TRUE)

#### A. DEFINE ALL-CAUSE MORTALITY/DATE OF DEATH ####

# merge with ONS dod
col_classes <- c("character", "character", "character", "character", "character", "character", "character", "character")
death_ons1 <- fread(paste(path_linkage_1, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons2 <- fread(paste(path_linkage_2, 'death_patient_24_004042_DM.txt', sep = '/'), colClasses = col_classes)
death_ons <- bind_rows(death_ons1, death_ons2)

rm(death_ons1, death_ons2)

length(which(cohort$id %in% death_ons$patid))
cat(paste('Patients with ONS information:', length(which(cohort$id %in% death_ons$patid)), '\n'), file = outcome_desc, append= TRUE)

death_ons %<>% 
  select(patid, dod) %>% 
  rename(ons_dod = dod)

cohort <- merge(cohort, death_ons, by.x = 'id', by.y = 'patid', all.x = TRUE)

cohort$ons_dod <- dmy(cohort$ons_dod)

dod_match <- cohort %>% 
  select(id, emis_dod, cprd_dod, ons_dod) %>% 
  mutate(same_dod = if_else(!is.na(emis_dod) & emis_dod == cprd_dod &
                            !is.na(cprd_dod) & cprd_dod == ons_dod &
                            !is.na(ons_dod), 1, 0),
         no_ons_dod = if_else( (!is.na(emis_dod) | !(is.na(cprd_dod))) & is.na(ons_dod), 1, 0)
         )

length(which(dod_match$same_dod == 1))  
cat(paste('Matching CPRD/EMIS/ONS dod:', length(which(dod_match$same_dod == 1)) , '\n'), file = outcome_desc, append= TRUE)

length(which(dod_match$no_ons_dod == 1))  
cat(paste('Missing ONS dod with EMIS or CPRD dod:', length(which(dod_match$no_ons_dod == 1)) , '\n'), file = outcome_desc, append= TRUE)


# set dod as earliest among CORD, and ONS
cohort$dod <- pmin(
  cohort$emis_dod,
  cohort$cprd_dod,
  cohort$ons_dod,
  na.rm = TRUE
)

# remove patients with issues in dod encoding
length(which(cohort$dod<cohort$birthdate))

cohort <- cohort %>% 
  mutate(dod_issue = if_else(!is.na(dod), # if patient died
                             if_else(dod<birthdate, 1, 0), # ensure date of death is after date of birth
                             0))

cat('Number of patients with issues in date of death encoding:', length(which(cohort$dod_issue==1)))
cat(paste('Number of patients with issues in date of death encoding:', length(which(cohort$dod_issue==1)), '\n'), file = outcome_desc, append = TRUE)

length(which(cohort$dod_issue==0))

cohort %<>% filter (dod_issue == 0) %>% select (-dod_issue)
length(unique(cohort$id))

summary(cohort$dod)

cat('Number of patients who experienced event on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)))
cat(paste('Number of patients who experienced event on day of cohort entry (removed):', length(which(cohort$dod == cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who experienced event prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)))
cat(paste('Number of patients who experienced event prior to cohort entry (removed):', length(which(cohort$dod < cohort$entry_date)), '\n'), file = outcome_desc, append = TRUE)

length(which(cohort$entry_date<cohort$itt_event_date | is.na(cohort$itt_event_date)))

cohort <- cohort %>% 
  filter (entry_date<dod | is.na(dod))

cat('Number of patients who died during study period:', length(which(!is.na(cohort$dod))))
cat(paste('Number of patients who died during study period:', length(which(!is.na(cohort$dod))), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who died during study period on SNRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'snri')))
cat(paste('Number of patients who died during study period on SNRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'snri')), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients who died during study period on SSRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'ssri')))
cat(paste('Number of patients who died during study period on SSRIs:', length(which(!is.na(cohort$dod) & cohort$trt == 'ssri')), '\n'), file = outcome_desc, append = TRUE)


#### DEFINING FOLLOW-UP ####

## ITT
# define ITT exit excluding event timing
cohort$itt_exit_date <- pmin(
  cohort$lcd, 
  cohort$regend,
  cohort$dod,
  study_follow_up_end,
  na.rm = TRUE
)

# define ITT event
# no event if no date of death or if died after ITT exit
cohort %<>%
  mutate (itt_event = if_else(is.na(dod) | (!is.na(dod) & dod > itt_exit_date), 0, 1),
          itt_event_date = if_else(itt_event == 1, dod, NA))

cohort$itt_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$itt_exit_date)), 'days')

## AT
# define AT exit excluding event timing
cohort$at_exit_date <- pmin(
  cohort$itt_exit_date,
  cohort$disc_date,
  cohort$switch_date,
  na.rm = TRUE
)

# define AT events
cohort <- cohort %>% 
  mutate (at_event = if_else(itt_event == 1 & itt_event_date <= at_exit_date, 1, 0),
          at_event_date = if_else(itt_event_date <= at_exit_date, itt_event_date, NA))

cohort$at_follow_up <- time_length(interval(as.Date(cohort$entry_date), as.Date(cohort$at_exit_date)), 'days')

# define discontinuation and switch flags
# disc: patients discontinued if they exited the cohort because their trt disc occurred first
cohort <- cohort %>% 
  mutate(
    disc = if_else(
      !is.na(disc_date) & at_exit_date == disc_date & # if there is a disc date and it led to cohort exit
        ((!is.na(itt_event_date) & at_exit_date != itt_event_date) | is.na(itt_event_date)), # and if there was an event then disc date was not on the day of the event
      1, 
      0
    )
  )

# switch: patients switched if they exited the cohort because their switch occurred first
cohort <- cohort %>% 
  mutate(
    switch = if_else(
      !is.na(switch_date) & at_exit_date == switch_date & 
        ((!is.na(itt_event_date) & at_exit_date != itt_event_date) | is.na(itt_event_date)), 
      1, 
      0
    )
  )


# examine follow up times
summary(cohort$itt_follow_up)
summary(cohort$at_follow_up)

test <- cohort %>% 
  filter(itt_follow_up <= 0) %>% 
  mutate(regend_before_exit = if_else(regend <= itt_exit_date, 1, 0)) %>% 
  select(id, entry_date, regend, itt_exit_date, at_exit_date, dod, itt_follow_up, at_follow_up, regend_before_exit)

# some patients have a negative ITT follow-up
# because their regimen ended prior to their cohort exit
length(which(cohort$itt_follow_up<0))
length(which(cohort$itt_follow_up==0)) 

cat('Number of patients with 0 or negative ITT follow-up (removed):', length(test$id))
cat(paste('Number of patients with 0 or negative ITT follow-up (removed):', length(test$id), '\n'), file = outcome_desc, append = TRUE)

cohort %<>% filter (itt_follow_up > 0)

summary(cohort$itt_follow_up)
summary(cohort$at_follow_up)

length(which(cohort$at_follow_up<=0)) 

# some patients have 0 AT follow-up because switched on the same day they entered
cat('Number of patients with 0 or negative AT follow-up (removed):', length(which(cohort$at_follow_up<=0)) )
cat(paste('Number of patients with 0 or negative AT follow-up (removed):', length(which(cohort$at_follow_up<=0)) , '\n'), file = outcome_desc, append = TRUE)

cohort %<>% filter (at_follow_up > 0)

# ensure no patients who discontinued have a follow-up less than the grace period or an event
test <- cohort %>% 
  filter(disc == 1)

summary(test$at_follow_up)
summary(test$at_event)

test <- cohort %>% # should be empty
  filter(disc == 1 & at_event == 1) %>% 
  select(id, entry_date, at_event_date, disc_date)

#### DEFINE OVERALL CENSOR DATE FROM TRT DISCONTINUATION OR SWITCH ####

cohort %<>%
  mutate (censor = if_else (switch == 1 | disc == 1, 1, 0),
          censor_date = if_else(censor == 1, pmin(switch_date, disc_date, na.rm = TRUE), NA))

table(cohort$censor)
length(which(is.na(cohort$censor_date)))

cat('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor))
cat(paste('Number of patients censored due to trt switch or discontinuation:', sum(cohort$censor), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients censored due to trt switch:', sum(cohort$switch))
cat(paste('Number of patients censored due to trt switch:', sum(cohort$switch), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients censored due to trt discontinuation:', sum(cohort$disc))
cat(paste('Number of patients censored due to trt discontinuation:', sum(cohort$disc), '\n'), file = outcome_desc, append = TRUE)

cat('Number of patients not censored:', sum(cohort$censor == 0))
cat(paste('Number of patients not censored:', sum(cohort$censor == 0), '\n'), file = outcome_desc, append = TRUE)

# describe trt switches
switched_to %<>% select(id, trt_seq)
cohort <- merge(cohort, switched_to, by='id', all.x = TRUE)
table(cohort[cohort$switch ==1, 'trt_seq'])
cat(paste('Switch SNRI to SSRI:', table(cohort[cohort$switch ==1, 'trt_seq']), '\n')[[1]], file = outcome_desc, append = TRUE)
cat(paste('Switch SSRI to SNRI:', table(cohort[cohort$switch ==1, 'trt_seq']), '\n')[[2]], file = outcome_desc, append = TRUE)

saveRDS(cohort, file = paste(path_cohort, 'antidepressant_cohort_outcome.rds', sep='/'))
