## ---------------------------
##
## Program: 8. Incidence rates, hazard ratios, and survival curves
##
## Purpose: 
## 1. Calculate incidence rates, rate ratios, and rate differences with and without weights for the outcome by exposure group.
## 2. Calculate cox proportional hazard ratios with and without rates for the outcome by exposure group and over time. 
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
times_dec <- times_dec[1:9] # remove last decile (100%)
times_dec

#### ITT ANALYSIS: ANALTYIC COHORT DATAFRAME ####

cohort_analytic_itt <- cohort

## Create analytic dataframe with long cohort split into intervals

# time: time to event or itt cohort exit
# counting_time: time in days between cohort entry and cohort exit
# event_at_time: whether patient experienced event at counting time

cohort_analytic_itt$event_at_time <- ifelse(is.na(cohort_analytic_itt$itt_event_date), 0, ifelse(cohort_analytic_itt$itt_exit_date == cohort_analytic_itt$itt_event_date, 1, 0))
cohort_analytic_itt$counting_time <- as.numeric(difftime(cohort_analytic_itt$itt_exit_date, cohort_analytic_itt$entry_date, units = "days"))
cohort_analytic_itt$event_counting_time <- as.numeric(difftime(cohort_analytic_itt$itt_event_date, cohort_analytic_itt$entry_date, units = 'days'))
cohort_analytic_itt <- cohort_analytic_itt[order(cohort_analytic_itt$counting_time),]

# note: will run into issues if someone's follow up ends on same day they enter (counting time of 0)
length(which(cohort$entry_date == cohort$switch_date))
length(which(cohort$entry_date == cohort$itt_exit_date))
length(which(cohort_analytic_itt$counting_time == 0))

cohort_analytic_itt %<>%
  filter(counting_time>0)

# split into deciles of censoring
cohort_analytic_itt$Tstart <- 0
cohort_analytic_itt <- survSplit(cohort_analytic_itt, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'event_at_time')
names(cohort_analytic_itt)[names(cohort_analytic_itt) == 'counting_time'] <- 'Tstop' 
cohort_analytic_itt$event_at_tstop <- if_else(is.na(cohort_analytic_itt$event_counting_time), 0, # indicator for having event for given time interval
                                              if_else(cohort_analytic_itt$Tstop == cohort_analytic_itt$event_counting_time, 1, 0)) 

for (i in 1:length(comorbidities)) {
  comorb_name <- comorbidities[i]
  cohort_analytic_itt$comorb <- if_else(
    cohort_analytic_itt$Tstart == 0,
    cohort_analytic_itt[[paste(comorb_name, 'base', sep = '_')]],
    if_else(
      cohort_analytic_itt$Tstart == times_dec[[1]],
      cohort_analytic_itt[[paste(comorb_name, 'd1', sep = '_')]],
      if_else(
        cohort_analytic_itt$Tstart == times_dec[[2]],
        cohort_analytic_itt[[paste(comorb_name, 'd2', sep = '_')]],
        if_else(
          cohort_analytic_itt$Tstart == times_dec[[3]],
          cohort_analytic_itt[[paste(comorb_name, 'd3', sep = '_')]],
          if_else(
            cohort_analytic_itt$Tstart == times_dec[[4]],
            cohort_analytic_itt[[paste(comorb_name, 'd4', sep = '_')]],
            if_else(
              cohort_analytic_itt$Tstart == times_dec[[5]],
              cohort_analytic_itt[[paste(comorb_name, 'd5', sep = '_')]],
              if_else(
                cohort_analytic_itt$Tstart == times_dec[[6]],
                cohort_analytic_itt[[paste(comorb_name, 'd6', sep = '_')]],
                if_else(
                  cohort_analytic_itt$Tstart == times_dec[[7]],
                  cohort_analytic_itt[[paste(comorb_name, 'd7', sep = '_')]],
                  if_else(
                    cohort_analytic_itt$Tstart == times_dec[[8]],
                    cohort_analytic_itt[[paste(comorb_name, 'd8', sep = '_')]],
                    cohort_analytic_itt[[paste(comorb_name, 'd9', sep = '_')]])
                )
              )
            )
          )
        )
      )
    )
  )
  names(cohort_analytic_itt)[names(cohort_analytic_itt) == 'comorb'] <- comorb_name
}

saveRDS(cohort_analytic_itt, file = paste(path_results, 'cohort_analytic_itt.rds', sep = '/'))


#### ITT ANALYSIS: INCIDENCE RATES ####

## Compute weighted person-time and events to calculate weighted IR
## for unstabilized and stabilized IPTW using ITT_cohort_analytic df
## note: IR calculated per YEAR

# convert event to a numeric for calculations
cohort_analytic_itt$itt_event <- as.numeric(cohort_analytic_itt$itt_event)

cohort_analytic_itt %<>%
  mutate (
    person_time = (Tstop - Tstart),
    iptw_person_time = person_time * iptw,
    siptw_person_time = person_time * siptw,
    iptw_event = event_at_tstop * iptw,
    siptw_event = event_at_tstop * siptw,
  )

itt_person_time <- cohort_analytic_itt %>% # get person time per id for visualization later
  group_by(id) %>% 
  mutate(itt_pt = sum(person_time),
         itt_pt_iptw = sum(iptw_person_time),
         itt_pt_siptw = sum(siptw_person_time)) %>% 
  select(id, itt_pt, itt_pt_iptw, itt_pt_siptw) %>% 
  slice(1)

saveRDS(itt_person_time, file = paste(path_cohort, 'itt_person_time.rds', sep = '/'))

itt_events <- cohort_analytic_itt %>% # get events per year for later
  group_by(year) %>% 
  mutate(total_events = sum(event_at_tstop),
         total_events_siptw = sum(siptw_event)) %>% 
  select(year, total_events, total_events_siptw) %>% 
  slice(1)

saveRDS(itt_events, file = paste(path_cohort, 'itt_events.rds', sep = '/'))
rm(itt_events, itt_person_time)

## Summary of incidence over entire study period

# quick overview table of events
table(cohort$itt_event, cohort$trt_dummy)
table(cohort_analytic_itt$event_at_tstop, cohort_analytic_itt$trt_dummy)

# no weights
itt_ir <- cohort_analytic_itt %>%
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
    total_follow_up_stab_yrs = sum(siptw_person_time), 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

itt_ir_iptw$model <- 'itt_ir_iptw'

#### ITT ANALYSIS: HAZARD RATIOS ####

## no weights
cox_no_weight <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_itt,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_no_weight, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) # check proportional hazard assumption, all p-values should be large (>0.05)
ph

png(filename = 'cox_no_weight_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_no_weight)
saveRDS(cox_no_weight, file = paste(path_results, 'cox_no_weight.rds', sep = '/'))

## IPTW weights

# unstab
cox_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = iptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_iptw)
saveRDS(cox_iptw, file = paste(path_results, 'cox_iptw.rds', sep = '/'))

# stab
cox_siptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_itt,
  weights = siptw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_siptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_siptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_siptw)
saveRDS(cox_siptw, file = paste(path_results, 'cox_siptw.rds', sep = '/'))
rm(cox_no_weight, cox_iptw, cox_siptw)

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

# note: will run into issues if someone's follow up ends on same day they enter (counting time of 0)
cohort_analytic_at %<>%
  filter(counting_time>0)

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

ipcw_weights <- cohort_long %>% dplyr::select(id, Tstart, ipcw, sipcw, ipcw_nl, sipcw_nl, ipcw_pooled, sipcw_pooled, ipcw_nl_mod, sipcw_nl_mod)

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
    ipcw_person_time = person_time * ipcw,
    ipcw_event = event_at_tstop * ipcw,
    
    sipcw_person_time = person_time * sipcw,
    sipcw_event = event_at_tstop * sipcw,
    
    iptw_ipcw_person_time = person_time * iptw * ipcw,
    iptw_ipcw_event = event_at_tstop * iptw * ipcw,
    
    siptw_sipcw_person_time = person_time * siptw * sipcw,
    siptw_sipcw_event = event_at_tstop * siptw * sipcw,
    
    
    # IPCW_nl and sIPCW_nl - non-lagged
    ipcw_nl_person_time = person_time * ipcw_nl,
    ipcw_nl_event = event_at_tstop * ipcw_nl,
    
    sipcw_nl_person_time = person_time * sipcw_nl,
    sipcw_nl_event = event_at_tstop * sipcw_nl,
    
    iptw_ipcw_nl_person_time = person_time * iptw * ipcw_nl,
    iptw_ipcw_nl_event = event_at_tstop * iptw * ipcw_nl,
    
    siptw_sipcw_nl_person_time = person_time * siptw * sipcw_nl,
    siptw_sipcw_nl_event = event_at_tstop * siptw * sipcw_nl,
    
    
    # IPCW_pooled and sIPCW_pooled
    ipcw_pooled_person_time = person_time * ipcw_pooled,
    ipcw_pooled_event = event_at_tstop * ipcw_pooled,
    
    sipcw_pooled_person_time = person_time * sipcw_pooled,
    sipcw_pooled_event = event_at_tstop * sipcw_pooled,
    
    iptw_ipcw_pooled_person_time = person_time * iptw * ipcw_pooled,
    iptw_ipcw_pooled_event = event_at_tstop * iptw * ipcw_pooled,
    
    siptw_sipcw_pooled_person_time = person_time * siptw * sipcw_pooled,
    siptw_sipcw_pooled_event = event_at_tstop * siptw * sipcw_pooled,
    
    # IPCW_nl_mod and sIPCW_nl_mod - non-lagged modified q1
    ipcw_nl_mod_person_time = person_time * ipcw_nl_mod,
    ipcw_nl_mod_event = event_at_tstop * ipcw_nl_mod,
    
    sipcw_nl_mod_person_time = person_time * sipcw_nl_mod,
    sipcw_nl_mod_event = event_at_tstop * sipcw_nl_mod,
    
    iptw_ipcw_nl_mod_person_time = person_time * iptw * ipcw_nl_mod,
    iptw_ipcw_nl_mod_event = event_at_tstop * iptw * ipcw_nl_mod,
    
    siptw_sipcw_nl_mod_person_time = person_time * siptw * sipcw_nl_mod,
    siptw_sipcw_nl_mod_event = event_at_tstop * siptw * sipcw_nl_mod
  )

at_person_time <- cohort_analytic_at %>% # get person time per id for later
  group_by(id) %>%
  mutate(
    at_pt = sum(person_time),
    at_pt_siptw_sipcw = sum(siptw_sipcw_person_time),
    at_pt_siptw_sipcw_nl = sum(siptw_sipcw_nl_person_time),
    at_pt_siptw_sipcw_pooled = sum(siptw_sipcw_pooled_person_time)
  ) %>%
  select(
    id,
    at_pt,
    at_pt_siptw_sipcw,
    at_pt_siptw_sipcw_nl,
    at_pt_siptw_sipcw_pooled
  ) %>%
  slice(1)

saveRDS(at_person_time, file = paste(path_cohort, 'at_person_time.rds', sep = '/'))

at_events <- cohort_analytic_at %>% # get events per year for later
  group_by(year) %>%
  mutate(
    total_events = sum(event_at_tstop),
    total_events_siptw_sipcw = sum(siptw_sipcw_event),
    total_events_siptw_sipcw_nl = sum(siptw_sipcw_nl_event),
    total_events_siptw_sipcw_pooled = sum(siptw_sipcw_pooled_event),
  ) %>%
  select(
    year,
    total_events,
    total_events_siptw_sipcw,
    total_events_siptw_sipcw_nl,
    total_events_siptw_sipcw_pooled
  ) %>% 
slice(1)

saveRDS(at_events, file = paste(path_cohort, 'at_events.rds', sep = '/'))
rm(at_person_time, at_events)

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
    total_follow_up_stab_yrs = sum(siptw_person_time), 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_iptw$model <- 'at_ir_iptw'

# lagged ipcw weights
at_ir_ipcw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_event),
    n_events_ref = sum(ipcw_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_event),
    total_follow_up_stab_yrs = sum(sipcw_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_event[trt_dummy == 0])/sum(sipcw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_event[trt_dummy == 1])/sum(sipcw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw$model <- 'at_ir_ipcw'

at_ir_ipcw_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_event),
    n_events_ref = sum(iptw_ipcw_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_event[trt_dummy == 0])/sum(siptw_sipcw_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_event[trt_dummy == 1])/sum(siptw_sipcw_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_iptw$model <- 'at_ir_ipcw_iptw'

# non-lagged ipcw weights
at_ir_ipcw_nl <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_nl_event),
    n_events_ref = sum(ipcw_nl_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_nl_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_nl_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_nl_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_nl_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_nl_event),
    total_follow_up_stab_yrs = sum(sipcw_nl_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_nl_event[trt_dummy == 0])/sum(sipcw_nl_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_nl_event[trt_dummy == 1])/sum(sipcw_nl_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_nl$model <- 'at_ir_ipcw_nl'

at_ir_ipcw_nl_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_nl_event),
    n_events_ref = sum(iptw_ipcw_nl_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_nl_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_nl_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_nl_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_nl_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_nl_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_nl_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_nl_event[trt_dummy == 0])/sum(siptw_sipcw_nl_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_nl_event[trt_dummy == 1])/sum(siptw_sipcw_nl_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_nl_iptw$model <- 'at_ir_ipcw_nl_iptw'


# pooled ipcw weights
at_ir_ipcw_pooled <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_pooled_event),
    n_events_ref = sum(ipcw_pooled_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_pooled_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_pooled_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_pooled_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_pooled_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_pooled_event),
    total_follow_up_stab_yrs = sum(sipcw_pooled_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_pooled_event[trt_dummy == 0])/sum(sipcw_pooled_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_pooled_event[trt_dummy == 1])/sum(sipcw_pooled_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_pooled$model <- 'at_ir_ipcw_pooled'

at_ir_ipcw_pooled_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_pooled_event),
    n_events_ref = sum(iptw_ipcw_pooled_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_pooled_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_pooled_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_pooled_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_pooled_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_pooled_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_pooled_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 0])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 1])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_pooled_iptw$model <- 'at_ir_ipcw_pooled_iptw'


# modified non-lagged ipcw weights
at_ir_ipcw_nl_mod <- cohort_analytic_at %>%
  summarise (
    n_events = sum(ipcw_nl_mod_event),
    n_events_ref = sum(ipcw_nl_mod_event[trt_dummy == 0]),
    n_events_comp = sum(ipcw_nl_mod_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(ipcw_nl_mod_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(sipcw_nl_mod_event),
    total_follow_up_stab_yrs = sum(sipcw_nl_mod_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(sipcw_nl_mod_event[trt_dummy == 0])/sum(sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(sipcw_nl_mod_event[trt_dummy == 1])/sum(sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_nl_mod$model <- 'at_ir_ipcw_nl_mod'

at_ir_ipcw_nl_mod_iptw <- cohort_analytic_at %>%
  summarise (
    n_events = sum(iptw_ipcw_nl_mod_event),
    n_events_ref = sum(iptw_ipcw_nl_mod_event[trt_dummy == 0]),
    n_events_comp = sum(iptw_ipcw_nl_mod_event[trt_dummy == 1]),
    total_follow_up_yrs = sum(iptw_ipcw_nl_mod_person_time)/365,
    IR = n_events / total_follow_up_yrs, 
    IR_ref = n_events_ref/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp = n_events_comp/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
    IRR = IR_comp / IR_ref,
    IRD = IR_comp - IR_ref,
    n_events_stab = sum(siptw_sipcw_nl_mod_event),
    total_follow_up_stab_yrs = sum(siptw_sipcw_nl_mod_person_time)/365, 
    IR_stab = n_events_stab / total_follow_up_stab_yrs, 
    IR_ref_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 0])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
    IR_comp_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 1])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
    IRR_stab = IR_comp_stab / IR_ref_stab,
    IRD_stab = IR_comp_stab - IR_ref_stab
  )

at_ir_ipcw_nl_mod_iptw$model <- 'at_ir_ipcw_nl_mod_iptw'

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
cox_lagged <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_lagged, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_lagged_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_lagged)
saveRDS(cox_lagged, file = paste(path_results, 'cox_lagged.rds', sep = '/'))

# stab
cox_stab_lagged <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_stab_lagged, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_stab_lagged_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_stab_lagged)
saveRDS(cox_stab_lagged, file = paste(path_results, 'cox_stab_lagged.rds', sep = '/'))

# iptw
cox_lagged_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_lagged_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_lagged_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_lagged_iptw)
saveRDS(cox_lagged_iptw, file = paste(path_results, 'cox_lagged_iptw.rds', sep = '/'))

# siptw
cox_stab_lagged_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_stab_lagged_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_stab_lagged_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_stab_lagged_iptw)
saveRDS(cox_stab_lagged_iptw, file = paste(path_results, 'cox_stab_lagged_iptw.rds', sep = '/'))

## NON-LAGGED WEIGHTS
# unstab
cox_nl <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_nl,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl)
saveRDS(cox_nl, file = paste(path_results, 'cox_nl.rds', sep = '/'))

# stab
cox_nl_stab <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_nl,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_stab, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_stab_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_stab)
saveRDS(cox_nl_stab, file = paste(path_results, 'cox_nl_stab.rds', sep = '/'))

# iptw
cox_nl_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_nl,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_iptw)
saveRDS(cox_nl_iptw, file = paste(path_results, 'cox_nl_iptw.rds', sep = '/'))

# siptw
cox_nl_stab_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_nl,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_stab_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_stab_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_stab_iptw)
saveRDS(cox_nl_stab_iptw, file = paste(path_results, 'cox_nl_stab_iptw.rds', sep = '/'))

## POOLED WEIGHTS 
# unstab
cox_pooled <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_pooled,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_pooled, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_pooled_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_pooled)
saveRDS(cox_pooled, file = paste(path_results, 'cox_pooled.rds', sep = '/'))

# stab
cox_pooled_stab <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_pooled,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_pooled_stab, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_pooled_stab_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_pooled_stab)
saveRDS(cox_pooled_stab, file = paste(path_results, 'cox_pooled_stab.rds', sep = '/'))

# iptw
cox_pooled_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_pooled,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_pooled_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_pooled_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_pooled_iptw)
saveRDS(cox_pooled_iptw, file = paste(path_results, 'cox_pooled_iptw.rds', sep = '/'))

# siptw
cox_pooled_stab_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_pooled,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_pooled_stab_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_pooled_stab_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_pooled_stab_iptw)
saveRDS(cox_pooled_stab_iptw, file = paste(path_results, 'cox_pooled_stab_iptw.rds', sep = '/'))


## MODIFIED NON-LAGGED WEIGHTS
# unstab
cox_nl_mod <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = ipcw_nl_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_mod, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_mod_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_mod)
saveRDS(cox_nl_mod, file = paste(path_results, 'cox_nl_mod.rds', sep = '/'))

# stab
cox_nl_mod_stab <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = sipcw_nl_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_mod_stab, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_mod_stab_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_mod_stab)
saveRDS(cox_nl_mod_stab, file = paste(path_results, 'cox_nl_mod_stab.rds', sep = '/'))

# iptw
cox_nl_mod_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = iptw * ipcw_nl_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_mod_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_mod_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_mod_iptw)
saveRDS(cox_nl_mod_iptw, file = paste(path_results, 'cox_nl_mod_iptw.rds', sep = '/'))

# siptw
cox_nl_mod_stab_iptw <- coxph(
  Surv(Tstart, Tstop, event_at_tstop) ~ trt_dummy,
  data = cohort_analytic_at,
  weights = siptw * sipcw_nl_mod,
  cluster = id,
  robust = TRUE
)

ph <- cox.zph(cox_nl_mod_stab_iptw, transform="km", terms=TRUE, singledf=FALSE, global=TRUE) 
ph

png(filename = 'cox_nl_mod_stab_iptw_ph.png', width = 3000, height = 3000, res = 500)
plot(ph)
dev.off()

summary(cox_nl_mod_stab_iptw)
saveRDS(cox_nl_mod_stab_iptw, file = paste(path_results, 'cox_nl_mod_stab_iptw.rds', sep = '/'))

rm(cox_at, cox_at_iptw, cox_at_siptw)
rm(cox_lagged, cox_stab_lagged, cox_lagged_iptw, cox_stab_lagged_iptw)
rm(cox_nl, cox_nl_iptw, cox_nl_stab, cox_nl_stab_iptw)
rm(cox_pooled, cox_pooled_stab, cox_pooled_iptw, cox_pooled_stab_iptw)

#### DATAFRAME WITH INCIDENCE RATES FOR EACH MODEL ####

incidence_rates <- bind_rows(
  itt_ir,
  itt_ir_iptw,
  at_ir,
  at_ir_iptw,
  at_ir_ipcw,
  at_ir_ipcw_iptw,
  at_ir_ipcw_nl,
  at_ir_ipcw_nl_iptw,
  at_ir_ipcw_pooled,
  at_ir_ipcw_pooled_iptw,
  at_ir_ipcw_nl_mod,
  at_ir_ipcw_nl_mod_iptw
)

incidence_rates %<>% relocate(model)

write_xlsx(incidence_rates, paste(path_results, 'incidence_rates.xlsx', sep ='/'))

