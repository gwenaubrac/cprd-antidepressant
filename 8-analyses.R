## ---------------------------
##
## Program: 8. Incidence rates, hazard ratios, and survival curves
##
## Purpose: 
## 1. Calculate incidence rates, rate ratios, and rate differences with and without weights for the outcome by exposure group.
## 2. Calculate cox proportional hazard ratios with and without rates for the outcome by exposure group. 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

# analysis: flex_grace_period, 90_day_grace_period
# male, female, young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed

analysis <- ''

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(ggplot2)
library(writexl)
library(lubridate)

#### DEFINE PATHS ####

path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 

if (analysis == 'main' | analysis == '') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/main" 
} else if (analysis == 'male') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/male" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/male" 
} else if (analysis == 'female') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/female" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/female" 
} else if (analysis == 'young') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/young" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/young" 
} else if (analysis == 'old') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/old" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/old" 
} else if (analysis == '2019') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2019" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2019" 
} else if (analysis == '2020') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2020" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2020" 
} else if (analysis == '2021') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2021" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2021" 
} else if (analysis == '2022') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/2022" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/2022" 
} else if (analysis == 'flex_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/flex_grace_period" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/sensitivity/flex_grace_period" 
} else if (analysis == '90_day_grace_period') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/sensitivity/90_day_grace_period" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/sensitivity/90_day_grace_period" 
} else if (analysis == 'depressed') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/depressed" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/depressed" 
} else if (analysis == 'not_depressed') {
  path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/subgroup/not_depressed" 
  path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/subgroup/not_depressed" 
} 


path_comorb_cprd <- "Z:/EPI/Protocol 24_004042/Gwen/data/comorbidities/Aurum codes comorbidities"
comorb_cprd_files <- list.files(path_comorb_cprd, pattern = '.xlsx', all.files = TRUE, full.names = TRUE)

setwd(path_results)

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))
cohort_long <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_ipcw.rds', sep = '/'))
covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_main, 'dec_comorb.rds', sep = '/'))

if (analysis == 'depressed' | analysis == 'not_depressed') {
  comorbidities <- comorbidities[!comorbidities %in% c('depression')]
  base_comorb <- base_comorb[!base_comorb %in% c('depression_base')]
  dec_comorb <- dec_comorb[!dec_comorb %in% c('depression_d1', 'depression_d2', 'depression_d3',
                                                  'depression_d4', 'depression_d5', 'depression_d6', 
                                                  'depression_d7', 'depression_d8', 'depression_d9')]
}

times_dec <- readRDS(paste(path_cohort, 'times_dec.rds', sep = '/'))
times_dec <- times_dec[1:9] # remove last decile (100%) for splitting later on
times_dec

table(cohort$trt)

#### ITT ANALYSIS: INCIDENCE RATES ####

## Compute weighted person-time and events to calculate weighted IR
## for unstabilized and stabilized IPTW using ITT_cohort_analytic df
## note: IR calculated per YEAR

# convert event to a numeric for calculations
cohort_analytic_itt <- cohort
cohort_analytic_itt$itt_event <- as.numeric(cohort_analytic_itt$itt_event)

cohort_analytic_itt %<>%
  mutate (
    person_time = as.numeric(itt_exit_date - entry_date),
    iptw_person_time = person_time * iptw,
    siptw_person_time = person_time * siptw,
    iptw_event = itt_event * iptw,
    siptw_event = itt_event * siptw,
  )

saveRDS(cohort_analytic_itt, file = paste(path_cohort, 'cohort_analytic_itt.rds', sep = '/'))

## Summary of incidence over entire study period

# quick overview table of events
table(cohort$itt_event, cohort$trt_dummy)
table(cohort_analytic_itt$itt_event, cohort_analytic_itt$trt_dummy)

# no weights
itt_ir <- cohort_analytic_itt %>%
  summarise (
    n = length(unique(id)),
    n_ref = n_distinct(id[trt_dummy == 0]), # number of patients in reference trt group
    n_comp= n_distinct(id[trt_dummy == 1]), # number of patients in comparator trt group
    n_events = sum(itt_event), # total events
    n_events_ref = n_distinct(id[itt_event == 1 & trt_dummy == 0]), # number of events in ref trt group
    n_events_comp = n_distinct(id[itt_event == 1 & trt_dummy == 1]), # number of events in comparator trt group
    total_follow_up_yrs = sum(person_time)/365, # total person-time (years)
    IR = n_events / total_follow_up_yrs, # overall incidence rate (per year)
    IR_ref = n_events_ref/sum(person_time[trt_dummy == 0]) * 365, # incidence rate in reference group (per year)
    IR_comp = n_events_comp/sum(person_time[trt_dummy == 1]) * 365, # incidence rate in comparator group (per year)
    IRR = IR_comp / IR_ref, # incidence rate ratio
    IRD = IR_comp - IR_ref # incidence rate difference
  )

itt_ir$model <- 'itt_ir'

# iptw and siptw weights
itt_ir_iptw <- cohort_analytic_itt %>%
  summarise (
    n_events = sum(iptw_event),
    n_events_ref = sum(iptw_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_event),
    total_follow_up_stab_yrs = sum(siptw_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

itt_ir_iptw$model <- 'itt_ir_iptw'

#### ITT ANALYSIS: HAZARD RATIOS ####

## no weights
cox_itt <- coxph(
  Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  cluster = id,
  robust = TRUE
)

# check the proportinal hazard assumption visually and with test
ph <- cox.zph(cox_itt, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) # check proportional hazard assumption, all p-values should be large (>0.05)
ph

png(filename = 'cox_itt_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt)
saveRDS(cox_itt, file = paste(path_results, 'cox_itt.rds', sep = '/'))

## IPTW weights

# unstab
cox_itt_iptw <- coxph(
  Surv(as.numeric(itt_exit_date, entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = iptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_itt_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_itt_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt_iptw)
saveRDS(cox_itt_iptw, file = paste(path_results, 'cox_itt_iptw.rds', sep = '/'))

# stab
cox_itt_siptw <- coxph(
  Surv(as.numeric(itt_exit_date, entry_date), itt_event) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = siptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_itt_siptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_itt_siptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_itt_siptw)
saveRDS(cox_itt_siptw, file = paste(path_results, 'cox_itt_siptw.rds', sep = '/'))
rm(cox_itt, cox_itt_iptw, cox_itt_siptw)

#### AT ANALYSIS: ANALYTIC COHORT DATAFRAME ####

cohort_analytic_at <- cohort

## Create analytic dataframe with long cohort split by event or at exit time

# time: time to event or cohort at cohort exit
# counting_time: time in days between cohort entry and time at exit
# event_at_time: whether patient experienced the event at counting time 

cohort_analytic_at$time <- pmin(cohort_analytic_at$at_exit_date, cohort_analytic_at$at_event_date, na.rm = TRUE)
cohort_analytic_at$event_at_time <- ifelse(is.na(cohort_analytic_at$at_event_date), 0, ifelse(cohort_analytic_at$time == cohort_analytic_at$at_event_date, 1, 0))
cohort_analytic_at$event_counting_time <- as.numeric(difftime(cohort_analytic_at$at_event_date, cohort_analytic_at$entry_date, units = 'days'))
cohort_analytic_at$counting_time <- as.numeric(difftime(cohort_analytic_at$time, cohort_analytic_at$entry_date, units = "days"))
cohort_analytic_at <- cohort_analytic_at[order(cohort_analytic_at$counting_time),]

cohort_analytic_at$Tstart <- 0
cohort_analytic_at <- survSplit(cohort_analytic_at, cut=times_dec, end="counting_time",
                                start="Tstart", event="event_at_time")

names(cohort_analytic_at)[names(cohort_analytic_at) == 'counting_time'] <- 'Tstop'
cohort_analytic_at$event_at_tstop <- if_else(is.na(cohort_analytic_at$event_counting_time), 0, # indicator for having event for given time interval
                                             if_else(cohort_analytic_at$Tstop == cohort_analytic_at$event_counting_time, 1, 0)) 

# retrieve covariate values at deciles
for (i in 1:length(comorbidities)) {
  comorb_name <- comorbidities[i]
  cohort_analytic_at$comorb <- if_else(
    cohort_analytic_at$Tstart == 0,
    cohort_analytic_at[[paste(comorb_name, 'base', sep = '_')]],
    if_else(
      cohort_analytic_at$Tstart == times_dec[[1]],
      cohort_analytic_at[[paste(comorb_name, 'd1', sep = '_')]],
      if_else(
        cohort_analytic_at$Tstart == times_dec[[2]],
        cohort_analytic_at[[paste(comorb_name, 'd2', sep = '_')]],
        if_else(
          cohort_analytic_at$Tstart == times_dec[[3]],
          cohort_analytic_at[[paste(comorb_name, 'd3', sep = '_')]],
          if_else(
            cohort_analytic_at$Tstart == times_dec[[4]],
            cohort_analytic_at[[paste(comorb_name, 'd4', sep = '_')]],
            if_else(
              cohort_analytic_at$Tstart == times_dec[[5]],
              cohort_analytic_at[[paste(comorb_name, 'd5', sep = '_')]],
              if_else(
                cohort_analytic_at$Tstart == times_dec[[6]],
                cohort_analytic_at[[paste(comorb_name, 'd6', sep = '_')]],
                if_else(
                  cohort_analytic_at$Tstart == times_dec[[7]],
                  cohort_analytic_at[[paste(comorb_name, 'd7', sep = '_')]],
                  if_else(
                    cohort_analytic_at$Tstart == times_dec[[8]],
                    cohort_analytic_at[[paste(comorb_name, 'd8', sep = '_')]],
                    cohort_analytic_at[[paste(comorb_name, 'd9', sep = '_')]])
                )
              )
            )
          )
        )
      )
    )
  )
  names(cohort_analytic_at)[names(cohort_analytic_at) == 'comorb'] <- comorb_name
}

table(cohort_analytic_at$event_at_tstop, cohort_analytic_at$trt)

# add IPCW corresponding to time interval
cohort_analytic_at$decTstart <- if_else(
  cohort_analytic_at$Tstart >= 0 & cohort_analytic_at$Tstart < times_dec[1],
  0,
  if_else(
    cohort_analytic_at$Tstart >= times_dec[1] & cohort_analytic_at$Tstart < times_dec[2],
    times_dec[1],
    if_else(
      cohort_analytic_at$Tstart >= times_dec[2] & cohort_analytic_at$Tstart < times_dec[3],
      times_dec[2],
      if_else(
        cohort_analytic_at$Tstart >= times_dec[3] & cohort_analytic_at$Tstart < times_dec[4],
        times_dec[3],
        if_else(
          cohort_analytic_at$Tstart >= times_dec[4] & cohort_analytic_at$Tstart < times_dec[5],
          times_dec[4],
          if_else(
            cohort_analytic_at$Tstart >= times_dec[5] & cohort_analytic_at$Tstart < times_dec[6],
            times_dec[5],
            if_else(
              cohort_analytic_at$Tstart >= times_dec[6] & cohort_analytic_at$Tstart < times_dec[7],
              times_dec[6],
              if_else(
                cohort_analytic_at$Tstart >= times_dec[7] & cohort_analytic_at$Tstart < times_dec[8],
                times_dec[7],
                if_else(
                  cohort_analytic_at$Tstart >= times_dec[8] & cohort_analytic_at$Tstart < times_dec[9],
                  times_dec[8], 
                  times_dec[9]
                )
              )
            )
          )
        )
      )
    )
  )
)

ipcw_weights <- cohort_long %>% dplyr::select(id, Tstart, ipcw_lag, sipcw_lag, ipcw_nonlag, sipcw_nonlag, ipcw_pool, sipcw_pool, ipcw_mod, sipcw_mod)

cohort_analytic_at <- cohort_analytic_at %>%
  dplyr::left_join(ipcw_weights, by = c('id', 'decTstart' = 'Tstart'), relationship = 'many-to-one')

saveRDS(cohort_analytic_at, file = paste(path_results, 'cohort_analytic_at.rds', sep = '/'))

#### AT ANALYSIS: INCIDENCE RATES ####

## Compute weighted person-time and events to calculate weighted IR
## for unstabilized and stabilized IPTW, IPCW, and combination weights
## for the AT analytic dataframe.

# convert event to a numeric for calculations
cohort_analytic_at$at_event <- as.numeric(cohort_analytic_at$at_event)

cohort_analytic_at %<>%
  mutate (
    person_time = (Tstop - Tstart),
    
    # IPTW and sIPTW
    iptw_person_time = person_time * iptw,
    iptw_event = event_at_tstop * iptw,
    
    siptw_person_time = person_time * siptw,
    siptw_event = event_at_tstop * siptw,
    
    
    # IPCW and sIPCW - lagged
    ipcw_lag_person_time = person_time * ipcw_lag,
    ipcw_lag_event = event_at_tstop * ipcw_lag,
    
    sipcw_lag_person_time = person_time * sipcw_lag,
    sipcw_lag_event = event_at_tstop * sipcw_lag,
    
    iptw_ipcw_lag_person_time = person_time * iptw * ipcw_lag,
    iptw_ipcw_lag_event = event_at_tstop * iptw * ipcw_lag,
    
    siptw_sipcw_lag_person_time = person_time * siptw * sipcw_lag,
    siptw_sipcw_lag_event = event_at_tstop * siptw * sipcw_lag,
    
    
    # ipcw_nonlag and sipcw_nonlag - non-lagged
    ipcw_nonlag_person_time = person_time * ipcw_nonlag,
    ipcw_nonlag_event = event_at_tstop * ipcw_nonlag,
    
    sipcw_nonlag_person_time = person_time * sipcw_nonlag,
    sipcw_nonlag_event = event_at_tstop * sipcw_nonlag,
    
    iptw_ipcw_nonlag_person_time = person_time * iptw * ipcw_nonlag,
    iptw_ipcw_nonlag_event = event_at_tstop * iptw * ipcw_nonlag,
    
    siptw_sipcw_nonlag_person_time = person_time * siptw * sipcw_nonlag,
    siptw_sipcw_nonlag_event = event_at_tstop * siptw * sipcw_nonlag,
    
    
    # ipcw_pool and sipcw_pool
    ipcw_pool_person_time = person_time * ipcw_pool,
    ipcw_pool_event = event_at_tstop * ipcw_pool,
    
    sipcw_pool_person_time = person_time * sipcw_pool,
    sipcw_pool_event = event_at_tstop * sipcw_pool,
    
    iptw_ipcw_pool_person_time = person_time * iptw * ipcw_pool,
    iptw_ipcw_pool_event = event_at_tstop * iptw * ipcw_pool,
    
    siptw_sipcw_pool_person_time = person_time * siptw * sipcw_pool,
    siptw_sipcw_pool_event = event_at_tstop * siptw * sipcw_pool,
    
    # ipcw_mod and sipcw_mod - non-lagged modified d1
    ipcw_mod_person_time = person_time * ipcw_mod,
    ipcw_mod_event = event_at_tstop * ipcw_mod,
    
    sipcw_mod_person_time = person_time * sipcw_mod,
    sipcw_mod_event = event_at_tstop * sipcw_mod,
    
    iptw_ipcw_mod_person_time = person_time * iptw * ipcw_mod,
    iptw_ipcw_mod_event = event_at_tstop * iptw * ipcw_mod,
    
    siptw_sipcw_mod_person_time = person_time * siptw * sipcw_mod,
    siptw_sipcw_mod_event = event_at_tstop * siptw * sipcw_mod
  )

## Summary of incidence over entire study period

# no weights
at_ir <- cohort_analytic_at %>%
  summarise (
    n = length(unique(id)),
    n_ref = n_distinct(id[trt_dummy == 0]), # number of patients in reference trt group
    n_comp= n_distinct(id[trt_dummy == 1]), # number of patients in comparator trt group
    n_events = sum(event_at_tstop), # total events
    n_events_ref = n_distinct(id[event_at_tstop == 1 & trt_dummy == 0]), # number of events in ref trt group
    n_events_comp = n_distinct(id[event_at_tstop == 1 & trt_dummy == 1]), # number of events in comparator trt group
    total_follow_up_yrs = sum(person_time)/365, # total person-time (years)
    IR = n_events / total_follow_up_yrs, # overall incidence rate (per year)
    IR_ref = n_events_ref/sum(person_time[trt_dummy == 0]) * 365, # incidence rate in reference group (per year)
    IR_comp = n_events_comp/sum(person_time[trt_dummy == 1]) * 365, # incidence rate in comparator group (per year)
    IRR = IR_comp / IR_ref, # incidence rate ratio
    IRD = IR_comp - IR_ref # incidence rate difference
  )

at_ir$model <- 'at_ir'

# iptw and siptw weights
at_ir_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_event),
    n_events_ref = sum(iptw_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_event),
    total_follow_up_stab_yrs = sum(siptw_person_time) / 365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw$model <- 'at_ir_iptw'

# lagged ipcw weights
at_ir_ipcw_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_lag_event),
    n_events_ref = sum(ipcw_lag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_lag_event),
    total_follow_up_stab_yrs = sum(sipcw_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_lag_event[trt_dummy == 0])/sum(sipcw_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_lag_event[trt_dummy == 1])/sum(sipcw_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_lag$model <- 'at_ir_ipcw_lag'

at_ir_iptw_ipcw_lag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_lag_event),
    n_events_ref = sum(iptw_ipcw_lag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_lag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_lag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_lag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_lag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_lag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_lag_event[trt_dummy == 0])/sum(siptw_sipcw_lag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_lag_event[trt_dummy == 1])/sum(siptw_sipcw_lag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_lag$model <- 'at_ir_iptw_ipcw_lag'

# non-lagged ipcw weights
at_ir_ipcw_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_nonlag_event),
    n_events_ref = sum(ipcw_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_nonlag_event),
    total_follow_up_stab_yrs = sum(sipcw_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_nonlag_event[trt_dummy == 0])/sum(sipcw_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_nonlag_event[trt_dummy == 1])/sum(sipcw_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_nonlag$model <- 'at_ir_ipcw_nonlag'

at_ir_iptw_ipcw_nonlag <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_nonlag_event),
    n_events_ref = sum(iptw_ipcw_nonlag_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_nonlag_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_nonlag_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_nonlag_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_nonlag_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_nonlag_event[trt_dummy == 0])/sum(siptw_sipcw_nonlag_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_nonlag_event[trt_dummy == 1])/sum(siptw_sipcw_nonlag_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_nonlag$model <- 'at_ir_iptw_ipcw_nonlag'


# pooled ipcw weights
at_ir_ipcw_pool <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_pool_event),
    n_events_ref = sum(ipcw_pool_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_pool_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_pool_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_pool_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_pool_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_pool_event),
    total_follow_up_stab_yrs = sum(sipcw_pool_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_pool_event[trt_dummy == 0])/sum(sipcw_pool_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pool_event[trt_dummy == 1])/sum(sipcw_pool_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_pool$model <- 'at_ir_ipcw_pool'

at_ir_iptw_ipcw_pool <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_pool_event),
    n_events_ref = sum(iptw_ipcw_pool_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_pool_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_pool_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_pool_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_pool_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_pool_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_pool_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_pool_event[trt_dummy == 0])/sum(siptw_sipcw_pool_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pool_event[trt_dummy == 1])/sum(siptw_sipcw_pool_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_pool$model <- 'at_ir_iptw_ipcw_pool'


# modified non-lagged ipcw weights
at_ir_ipcw_mod <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_mod_event),
    n_events_ref = sum(ipcw_mod_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_mod_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_mod_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_mod_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_mod_event),
    total_follow_up_stab_yrs = sum(sipcw_mod_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_mod_event[trt_dummy == 0])/sum(sipcw_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_mod_event[trt_dummy == 1])/sum(sipcw_mod_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_mod$model <- 'at_ir_ipcw_mod'

at_ir_iptw_ipcw_mod <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_mod_event),
    n_events_ref = sum(iptw_ipcw_mod_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_mod_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_mod_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_mod_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_mod_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_mod_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_mod_event[trt_dummy == 0])/sum(siptw_sipcw_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_mod_event[trt_dummy == 1])/sum(siptw_sipcw_mod_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw_ipcw_mod$model <- 'at_ir_iptw_ipcw_mod'

#### AT ANALYSIS: HAZARD RATIOS ####

## NO WEIGHTS
cox_at <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at)
saveRDS(cox_at, file = paste(path_results, 'cox_at.rds', sep = '/'))

## IPTW weights
# unstab
cox_at_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw)
saveRDS(cox_at_iptw, file = paste(path_results, 'cox_at_iptw.rds', sep = '/'))

# stab
cox_at_siptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw)
saveRDS(cox_at_siptw, file = paste(path_results, 'cox_at_siptw.rds', sep = '/'))

## LAGGED WEIGHTS
# unstab
cox_at_ipcw_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_lag)
saveRDS(cox_at_ipcw_lag, file = paste(path_results, 'cox_at_ipcw_lag.rds', sep = '/'))

# stab
cox_at_sipcw_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_lag)
saveRDS(cox_at_sipcw_lag, file = paste(path_results, 'cox_at_sipcw_lag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_lag)
saveRDS(cox_at_iptw_ipcw_lag, file = paste(path_results, 'cox_at_iptw_ipcw_lag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_lag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_lag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_lag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_lag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_lag)
saveRDS(cox_at_siptw_sipcw_lag, file = paste(path_results, 'cox_at_siptw_sipcw_lag.rds', sep = '/'))

## NON-LAGGED WEIGHTS
# unstab
cox_at_icpw_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_icpw_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_icpw_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_icpw_nonlag)
saveRDS(cox_at_icpw_nonlag, file = paste(path_results, 'cox_at_ipcw_nonlag.rds', sep = '/'))

# stab
cox_at_sipcw_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_nonlag.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_nonlag)
saveRDS(cox_at_sipcw_nonlag, file = paste(path_results, 'cox_at_sipcw_nonlag.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_nonlag)
saveRDS(cox_at_iptw_ipcw_nonlag, file = paste(path_results, 'cox_at_iptw_ipcw_nonlag.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_nonlag <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_nonlag,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_nonlag, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_nonlag_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_nonlag)
saveRDS(cox_at_siptw_sipcw_nonlag, file = paste(path_results, 'cox_at_siptw_sipcw_nonlag.rds', sep = '/'))

## POOLED WEIGHTS 
# unstab
cox_at_ipcw_pool <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_pool,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_pool, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_pool_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_pool)
saveRDS(cox_at_ipcw_pool, file = paste(path_results, 'cox_at_ipcw_pool.rds', sep = '/'))

# stab
cox_at_sipcw_pool <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_pool,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_pool, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_pool_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_pool)
saveRDS(cox_at_sipcw_pool, file = paste(path_results, 'cox_at_sipcw_pool.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_pool <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_pool,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_pool, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_pool_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_pool)
saveRDS(cox_at_iptw_ipcw_pool, file = paste(path_results, 'cox_at_iptw_ipcw_pool.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_pool <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_pool,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_pool, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_pool_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_pool)
saveRDS(cox_at_siptw_sipcw_pool, file = paste(path_results, 'cox_at_siptw_sipcw_pool.rds', sep = '/'))


## MODIFIED NON-LAGGED WEIGHTS
# unstab
cox_at_ipcw_mod <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_ipcw_mod, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_ipcw_mod_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_ipcw_mod)
saveRDS(cox_at_ipcw_mod, file = paste(path_results, 'cox_at_ipcw_mod.rds', sep = '/'))

# stab
cox_at_sipcw_mod <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_sipcw_mod, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_sipcw_mod_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_sipcw_mod)
saveRDS(cox_at_sipcw_mod, file = paste(path_results, 'cox_at_sipcw_mod.rds', sep = '/'))

# iptw
cox_at_iptw_ipcw_mod <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_iptw_ipcw_mod, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_iptw_ipcw_mod.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_iptw_ipcw_mod)
saveRDS(cox_at_iptw_ipcw_mod, file = paste(path_results, 'cox_at_iptw_ipcw_mod.rds', sep = '/'))

# siptw
cox_at_siptw_sipcw_mod <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_at_siptw_sipcw_mod, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_at_siptw_sipcw_mod_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_at_siptw_sipcw_mod)
saveRDS(cox_at_siptw_sipcw_mod, file = paste(path_results, 'cox_at_siptw_sipcw_mod.rds', sep = '/'))

#### DATAFRAME WITH INCIDENCE RATES FOR EACH MODEL ####

incidence_rates <- bind_rows(
  itt_ir,
  itt_ir_iptw,
  at_ir,
  at_ir_iptw,
  at_ir_ipcw_lag,
  at_ir_iptw_ipcw_lag,
  at_ir_ipcw_nonlag,
  at_ir_iptw_ipcw_nonlag,
  at_ir_ipcw_pool,
  at_ir_iptw_ipcw_pool,
  at_ir_ipcw_mod,
  at_ir_iptw_ipcw_mod
)

incidence_rates %<>% relocate(model)

write_xlsx(incidence_rates, paste(path_results, 'incidence_rates.xlsx', sep ='/'))

