## ---------------------------
##
## Program: 9. Bootstrapping 
##
## Purpose: Bootstrap cohort to obtain percentile confidence intervals for IR, IRR, and IRD.
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-09
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

library(parallel)
library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(writexl)
library(lubridate)
library(dplyr)
library(magrittr)
library(fastDummies)
library(ggplot2)
library(tidyr)
library(splines)
library(boot)

set.seed(123) # for reproducibility

#### DEFINE PATHS ####

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


path_cprdA <- "Z:/EPI/Protocol 24_004042/dataA"
path_cprdB <- "Z:/EPI/Protocol 24_004042/dataB"
path_cprdC <- "Z:/EPI/Protocol 24_004042/dataC (no followup)" 
path_exposure <- "Z:/EPI/Protocol 24_004042/Gwen/data/exposures/antidepressants" 
path_comorb_cprd <- "Z:/EPI/Protocol 24_004042/Gwen/data/comorbidities/Aurum codes comorbidities"

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))

covariates <- readRDS(file = paste(path_cohort, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_cohort, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_cohort, 'base_comorb.rds', sep = '/'))
quart_comorb <- readRDS(file = paste(path_cohort, 'quart_comorb.rds', sep = '/'))

comorbidities <- comorbidities[!comorbidities %in% c('hypokalemia', 'hypocalcemia', 'hypomagnesemia')]
base_comorb <- base_comorb[!base_comorb %in% c('hypokalemia_base', 'hypocalcemia_base', 'hypomagnesemia_base')]
quart_comorb <- quart_comorb[!quart_comorb %in% c('hypokalemia_Q1', 'hypocalcemia_Q1', 'hypomagnesemia_Q1', 
                                                  'hypokalemia_Q2', 'hypocalcemia_Q2', 'hypomagnesemia_Q2',
                                                  'hypokalemia_Q3', 'hypocalcemia_Q3', 'hypomagnesemia_Q3')]

covariates

#### BOOTSTRAP FUNCTION ####

# test run once not as a function
# to extract column names

bs <- function(data, indices) {
  
  d <- data[indices,]
  d %<>% mutate(boot_id = row_number())
  
  #### CALCUALTE IPTW ####
  
  base_variables <- c(covariates, base_comorb)
  
  if (analysis == 'male' | analysis == 'female') {
    base_variables <- base_variables[!base_variables %in% c('sex')]
  } else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
    base_variables <- base_variables[!base_variables %in% c('year', 'month_year')]
  } else if (analysis == 'depressed' | analysis == 'not_depressed') {
    base_variables <- base_variables[!base_variables %in% c('depression_base')]
  }
  
  d <- d %>%
    filter (!is.na(sex))
  
  model <- reformulate(base_variables, 'trt_dummy')
  
  iptw_fit <- glm(model, 
                  family = binomial(), 
                  data = d)
  
  prop_score <- if_else(d$trt_dummy == 0, 
                        1 - predict(iptw_fit, type = 'response'),
                        predict(iptw_fit, type = 'response'))
  
  d$ps <- prop_score
  d$iptw <- 1/prop_score
  
  numer_fit <- glm(trt_dummy~1, family = binomial(), data = d)
  pn_trt <- predict(numer_fit, type = 'response')
  
  denom_fit <- glm(model,
                   family = binomial(),
                   data = d)
  
  
  pd_trt <- predict(denom_fit, type = 'response')

  d$siptw <- if_else(d$trt_dummy == 0, 
                          ((1-pn_trt) / (1-pd_trt)),
                          pn_trt/pd_trt)
  
  #### CALCULATE IPCW ####
  
  d_long <- d

  if (analysis == 'male' | analysis == 'female') {
    variables <- variables[!variables %in% c('sex')]
  } else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
    variables <- variables[!variables %in% c('year', 'month_year')]
  } else if (analysis == 'depressed' | analysis == 'not_depressed') {
    variables <- variables[!variables %in% c('depression')]
  }
  
  base_variables <- c(covariates, base_comorb)

  if (analysis == 'male' | analysis == 'female') {
    base_variables <- base_variables[!base_variables %in% c('sex')]
  } else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
    base_variables <- base_variables[!base_variables %in% c('year', 'month_year')]
  } else if (analysis == 'depressed' | analysis == 'not_depressed') {
    base_variables <- base_variables[!base_variables %in% c('depression_base')]
  }
  
  d_long$time <- d_long$at_exit_date
  d_long$uncensored_at_time <- if_else(is.na(d_long$censor_date), 1, 
                                            if_else(d_long$time == d_long$censor_date, 0, 1)) 
  d_long$counting_time <- as.numeric(difftime(d_long$time, d_long$entry_date, units = 'days'))
  d_long$censor_counting_time <- as.numeric(difftime(d_long$censor_date, d_long$entry_date, units = 'days'))
  d_long <- d_long[order(d_long$counting_time),]
 
  times_quart <- readRDS(paste(path_cohort, 'times_quart.rds', sep = '/'))
  
  d_long %<>%
    filter(counting_time>0)
  
  d_long$Tstart <- 0
  d_long <- survSplit(d_long, cut = times_quart, end = 'counting_time', start = 'Tstart', event = 'uncensored_at_time')
  names(d_long)[names(d_long) == 'counting_time'] <- 'Tstop' 
  d_long$uncensored_at_tstop <- if_else(is.na(d_long$censor_counting_time), 1, 
                                             if_else(d_long$Tstop == d_long$censor_counting_time, 0, 1)) 
  for (i in 1:length(comorbidities)) {
    comorb_name <- comorbidities[i]
    d_long$comorb <- if_else(
      d_long$Tstart == 0,
      d_long[[paste(comorb_name, 'base', sep = '_')]],
      if_else(
        d_long$Tstart == times_quart[[1]],
        d_long[[paste(comorb_name, 'Q1', sep = '_')]],
        if_else(
          d_long$Tstart == times_quart[[2]],
          d_long[[paste(comorb_name, 'Q2', sep = '_')]],
          d_long[[paste(comorb_name, 'Q3', sep = '_')]]
        )
      )
    )
    names(d_long)[names(d_long) == 'comorb'] <- comorb_name
  }
  
  # lagged model
  
  d_long_ref <- d_long %>% 
    filter(trt_dummy == 0)
  
  subset_Q1 <- subset(d_long_ref, Tstart == 0)
  subset_Q2 <- subset(d_long_ref, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_ref, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_ref, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')
  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  
  predict_Q2 <- predict(uncensored_in_Q1, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q2, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q3, subset_Q4, type = 'response')
  
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens'] <- 1
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4
  
  d_long_ref <- d_long_ref %>%
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 
  
  d_long_ref$ipcw <- 1/d_long_ref$p_uncens_cum
  d_long_ref$sipcw <- d_long_ref$p_uncens_base/d_long_ref$p_uncens_cum
  
  d_long$ipcw <- 0
  d_long$sipcw <- 0
  d_long[d_long$trt_dummy == 0, 'ipcw'] <- d_long_ref$ipcw
  d_long[d_long$trt_dummy == 0, 'sipcw'] <- d_long_ref$sipcw
  
  d_long_comp <- d_long %>% 
    filter(trt_dummy == 1)
  
  subset_Q1 <- subset(d_long_comp, Tstart == 0)
  subset_Q2 <- subset(d_long_comp, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_comp, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_comp, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')

  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  
  predict_Q2 <- predict(uncensored_in_Q1, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q2, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q3, subset_Q4, type = 'response')
  
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens'] <- 1
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4
  
  d_long_comp <- d_long_comp %>% 
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 
  
  d_long_comp$ipcw <- 1/d_long_comp$p_uncens_cum
  d_long_comp$sipcw <- d_long_comp$p_uncens_base/d_long_comp$p_uncens_cum
  
  d_long[d_long$trt_dummy == 1, 'ipcw'] <- d_long_comp$ipcw
  d_long[d_long$trt_dummy == 1, 'sipcw'] <- d_long_comp$sipcw
  
  # non-lagged model
  
  d_long_ref <- d_long %>% 
    filter(trt_dummy == 0)
  
  subset_Q1 <- subset(d_long_ref, Tstart == 0)
  subset_Q2 <- subset(d_long_ref, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_ref, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_ref, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')
  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  uncensored_in_Q4 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q4, family = binomial)
  
  predict_Q1 <- predict(uncensored_in_Q1, subset_Q1, type = 'response')
  predict_Q2 <- predict(uncensored_in_Q2, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q3, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q4, subset_Q4, type = 'response')
  
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens'] <- predict_Q1
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4

  d_long_ref <- d_long_ref %>%
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 

  d_long_ref$ipcw <- 1/d_long_ref$p_uncens_cum
  d_long_ref$sipcw <- d_long_ref$p_uncens_base/d_long_ref$p_uncens_cum

  d_long$ipcw_nl <- 0
  d_long$sipcw_nl <- 0
  d_long[d_long$trt_dummy == 0, 'ipcw_nl'] <- d_long_ref$ipcw
  d_long[d_long$trt_dummy == 0, 'sipcw_nl'] <- d_long_ref$sipcw
  
  d_long_comp <- d_long %>% 
    filter(trt_dummy == 1)
  
  subset_Q1 <- subset(d_long_comp, Tstart == 0)
  subset_Q2 <- subset(d_long_comp, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_comp, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_comp, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')
  
  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  uncensored_in_Q4 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q4, family = binomial)
  
  predict_Q1 <- predict(uncensored_in_Q1, subset_Q1, type = 'response')
  predict_Q2 <- predict(uncensored_in_Q2, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q3, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q4, subset_Q4, type = 'response')
  
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens'] <- predict_Q1
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4

  d_long_comp <- d_long_comp %>% 
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 

  d_long_comp$ipcw <- 1/d_long_comp$p_uncens_cum
  d_long_comp$sipcw <- d_long_comp$p_uncens_base/d_long_comp$p_uncens_cum

  d_long[d_long$trt_dummy == 1, 'ipcw_nl'] <- d_long_comp$ipcw
  d_long[d_long$trt_dummy == 1, 'sipcw_nl'] <- d_long_comp$sipcw
  
  # pooled model 
  
  d_long_ref <- d_long %>% 
    filter(trt_dummy == 0) %>% 
    group_by(boot_id) %>% 
    mutate (quart = as.factor(row_number()))

  baseline_uncensored <- glm(reformulate('1', 'uncensored_at_tstop'), data = d_long_ref, family = binomial())
  d_long_ref$p_uncens_base <- predict(baseline_uncensored, type = 'response')

  variables_pooled <- c(variables, 'quart')
  
  uncensored_pooled_ref <- glm(reformulate(variables_pooled, 'uncensored_at_tstop'), data = d_long_ref, family = binomial) 
  
  d_long_ref$p_uncens_pooled <- predict(uncensored_pooled_ref, type = 'response')
  
  d_long_ref <- d_long_ref %>%  
    group_by(boot_id) %>% 
    mutate (p_uncens_cum_pooled = cumprod(p_uncens_pooled)) 
  
  d_long_ref$ipcw <- 1/d_long_ref$p_uncens_cum_pooled
  d_long_ref$sipcw <- d_long_ref$p_uncens_base/d_long_ref$p_uncens_cum_pooled
  
  d_long$ipcw_pooled <- 0
  d_long$sipcw_pooled <- 0
  
  d_long[d_long$trt_dummy == 0, 'ipcw_pooled'] <- d_long_ref$ipcw
  d_long[d_long$trt_dummy == 0, 'sipcw_pooled'] <- d_long_ref$sipcw
  
  d_long_comp <- d_long %>% 
    filter(trt_dummy == 1) %>% 
    group_by(boot_id) %>% 
    mutate (quart = as.factor(row_number()))
  
  baseline_uncensored <- glm(reformulate('1', 'uncensored_at_tstop'), data = d_long_comp, family = binomial())
  d_long_comp$p_uncens_base <- predict(baseline_uncensored, type = 'response')
  
  variables_pooled <- c(variables, 'quart')
  
  uncensored_pooled_comp <- glm(reformulate(variables_pooled, 'uncensored_at_tstop'), data = d_long_comp, family = binomial)
  d_long_comp$p_uncens_pooled <- predict(uncensored_pooled_comp, type = 'response')

  d_long_comp %<>% 
    group_by(boot_id) %>% 
    mutate (p_uncens_cum_pooled = cumprod(p_uncens_pooled)) 
  
  d_long_comp$ipcw <- 1/d_long_comp$p_uncens_cum_pooled
  d_long_comp$sipcw <- d_long_comp$p_uncens_base/d_long_comp$p_uncens_cum_pooled

  d_long[d_long$trt_dummy == 1, 'ipcw_pooled'] <- d_long_comp$ipcw
  d_long[d_long$trt_dummy == 1, 'sipcw_pooled'] <- d_long_comp$sipcw
  
  # modified non-lagged weights
  
  d_long_ref <- cohort_long %>% 
    filter(trt_dummy == 0)
  
  subset_Q1 <- subset(d_long_ref, Tstart == 0)
  subset_Q2 <- subset(d_long_ref, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_ref, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_ref, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')
  
  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  uncensored_in_Q4 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q4, family = binomial)
  
  predict_Q1 <- predict(uncensored_in_Q1, subset_Q1, type = 'response')
  predict_Q2 <- predict(uncensored_in_Q2, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q3, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q4, subset_Q4, type = 'response')
  
  d_long_ref[d_long_ref$Tstart==0, 'p_uncens'] <- predict_Q1
  d_long_ref[d_long_ref$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_ref[d_long_ref$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_ref[d_long_ref$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4
  
  d_long_ref <- d_long_ref %>%
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 
  
  d_long_ref <- d_long_ref %>%
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = if_else(row_number() == 1, 1, p_uncens_cum))
  
  d_long_ref$ipcw <- 1/d_long_ref$p_uncens_cum
  d_long_ref$sipcw <- d_long_ref$p_uncens_base/d_long_ref$p_uncens_cum

  cohort_long$ipcw_nl_mod <- 0
  cohort_long$sipcw_nl_mod <- 0
  cohort_long[cohort_long$trt_dummy == 0, 'ipcw_nl_mod'] <- d_long_ref$ipcw
  cohort_long[cohort_long$trt_dummy == 0, 'sipcw_nl_mod'] <- d_long_ref$sipcw
  
  d_long_comp <- cohort_long %>% 
    filter(trt_dummy == 1)
  
  subset_Q1 <- subset(d_long_comp, Tstart == 0)
  subset_Q2 <- subset(d_long_comp, Tstart == times_quart[[1]])
  subset_Q3 <- subset(d_long_comp, Tstart == times_quart[[2]])
  subset_Q4 <- subset(d_long_comp, Tstart == times_quart[[3]])
  
  baseline_uncensored_q1 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q1, family = binomial())
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored, type = 'response')
  
  baseline_uncensored_Q2 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q2, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_Q2, type = 'response')
  
  baseline_uncensored_Q3 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q3, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_Q3, type = 'response')
  
  baseline_uncensored_Q4 <- glm(reformulate('1', 'uncensored_at_tstop'), data = subset_Q4, family = binomial())
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_Q4, type = 'response')
  
  
  uncensored_in_Q1 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q1, family = binomial)
  uncensored_in_Q2 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q2, family = binomial)
  uncensored_in_Q3 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q3, family = binomial)
  uncensored_in_Q4 <- glm(reformulate(variables, 'uncensored_at_tstop'), data = subset_Q4, family = binomial)
  
  predict_Q1 <- predict(uncensored_in_Q1, subset_Q1, type = 'response')
  predict_Q2 <- predict(uncensored_in_Q2, subset_Q2, type = 'response')
  predict_Q3 <- predict(uncensored_in_Q3, subset_Q3, type = 'response')
  predict_Q4 <- predict(uncensored_in_Q4, subset_Q4, type = 'response')
  
  d_long_comp[d_long_comp$Tstart==0, 'p_uncens'] <- predict_Q1
  d_long_comp[d_long_comp$Tstart==times_quart[[1]], 'p_uncens'] <- predict_Q2
  d_long_comp[d_long_comp$Tstart==times_quart[[2]], 'p_uncens'] <- predict_Q3
  d_long_comp[d_long_comp$Tstart==times_quart[[3]], 'p_uncens'] <- predict_Q4
  
  d_long_comp <- d_long_comp %>% 
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = cumprod(p_uncens)) 
  
  d_long_comp <- d_long_comp %>%
    group_by(boot_id) %>% 
    mutate (p_uncens_cum = if_else(row_number() == 1, 1, p_uncens_cum))
  
  d_long_comp$ipcw <- 1/d_long_comp$p_uncens_cum
  d_long_comp$sipcw <- d_long_comp$p_uncens_base/d_long_comp$p_uncens_cum
  
  d_long[d_long$trt_dummy == 1, 'ipcw_nl_mod'] <- d_long_comp$ipcw
  d_long[d_long$trt_dummy == 1, 'sipcw_nl_mod'] <- d_long_comp$sipcw

  
  ### 1. ITT ANALYSES ###
  
  ## A. PREPARING DATA ##
  
  itt_d <- d
  itt_d$event_at_time <- if_else(is.na(itt_d$itt_event_date), 0, if_else(itt_d$itt_exit_date == itt_d$itt_event_date, 1, 0))
  itt_d$counting_time <- as.numeric(difftime(itt_d$itt_exit_date, itt_d$entry_date, units = "days"))
  itt_d$event_counting_time <- as.numeric(difftime(itt_d$itt_event_date, itt_d$entry_date, units = 'days'))
  itt_d$censor_counting_time <- as.numeric(difftime(itt_d$at_exit_date, itt_d$entry_date, units = 'days'))
  itt_d <- itt_d[order(itt_d$counting_time),]
  
  itt_d %<>%
    filter(counting_time>0)
  
  itt_d$Tstart <- 0
  itt_d <- survSplit(itt_d, cut = times_quart, end = 'counting_time', start = 'Tstart', event = 'event_at_time')
  names(itt_d)[names(itt_d) == 'counting_time'] <- 'Tstop' 
  itt_d$event_at_tstop <- if_else(is.na(itt_d$event_counting_time), 0, 
                                                if_else(itt_d$Tstop == itt_d$event_counting_time, 1, 0)) 
  
  for (i in 1:length(comorbidities)) {
    comorb_name <- comorbidities[i]
    itt_d$comorb <- if_else(
      itt_d$Tstart == 0,
      itt_d[[paste(comorb_name, 'base', sep = '_')]],
      if_else(
        itt_d$Tstart == times_quart[[1]],
        itt_d[[paste(comorb_name, 'Q1', sep = '_')]],
        if_else(
          itt_d$Tstart == times_quart[[2]],
          itt_d[[paste(comorb_name, 'Q2', sep = '_')]],
          itt_d[[paste(comorb_name, 'Q3', sep = '_')]]
        )
      )
    )
    names(itt_d)[names(itt_d) == 'comorb'] <- comorb_name
  }
  
  ## B. INCIDENCE RATES ##
  
  itt_d$itt_event <- as.numeric(itt_d$itt_event)
  
  itt_d %<>%
    mutate (
      person_time = (Tstop - Tstart),
      iptw_person_time = person_time * iptw,
      siptw_person_time = person_time * siptw,
      iptw_event = event_at_tstop * iptw,
      siptw_event = event_at_tstop * siptw,
    )
  
  ## over entire study
  
  # no weights
  itt_ir <- itt_d %>%
    summarise (
      IR = sum(event_at_tstop) / (sum(person_time)/365), 
      IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
      IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365, 
      IRR = IR_comp / IR_ref, 
      IRD = IR_comp - IR_ref,
      IR_stab = 0, 
      IR_ref_stab = 0,
      IR_comp_stab = 0,
      IRR_stab = 0,
      IRD_stab = 0
    )
  
  itt_ir <- data.frame(t(itt_ir))
  colnames(itt_ir) <- 'itt.overall'
  itt_ir$row_names <- rownames(itt_ir)
  itt_ir <- itt_ir %>% relocate(row_names)
  
  # iptw and siptw weights
  itt_ir_iptw <- itt_d %>%
    summarise (
      IR = sum(iptw_event) / sum(iptw_person_time)/365, 
      IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_event) / sum(siptw_person_time), 
      IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  itt_ir_iptw <- data.frame(t(itt_ir_iptw))
  colnames(itt_ir_iptw) <- 'itt.overall.iptw'
  itt_ir_iptw$row_names <- rownames(itt_ir_iptw)
  itt_ir_iptw <- itt_ir_iptw %>% relocate(row_names)
  
  ir_results <- left_join(itt_ir, itt_ir_iptw, by = 'row_names')
  
  ## over time
  
  # no weights
  itt_ir_time <- itt_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(event_at_tstop) / (sum(person_time)/365), 
      IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
      IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365, 
      IRR = IR_comp / IR_ref, 
      IRD = IR_comp - IR_ref,
      IR_stab = 0, 
      IR_ref_stab = 0,
      IR_comp_stab = 0,
      IRR_stab = 0,
      IRD_stab = 0
    )
  
  itt_ir_time <- data.frame(t(itt_ir_time))
  colnames(itt_ir_time) <- itt_ir_time[1,]
  itt_ir_time <- itt_ir_time[-1,]
  colnames(itt_ir_time) <- paste('itt', colnames(itt_ir_time), sep = '.')
  itt_ir_time$row_names <- rownames(itt_ir_time)
  itt_ir_time <- itt_ir_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, itt_ir_time, by = 'row_names')
  
  # iptw and siptw weights
  itt_ir_iptw_time <- itt_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_event) / sum(iptw_person_time)/365, 
      IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_event) / sum(siptw_person_time), 
      IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  itt_ir_iptw_time <- data.frame(t(itt_ir_iptw_time))
  colnames(itt_ir_iptw_time) <- itt_ir_iptw_time[1,]
  itt_ir_iptw_time <- itt_ir_iptw_time[-1,]
  colnames(itt_ir_iptw_time) <- paste('itt', colnames(itt_ir_iptw_time), sep = '.')
  itt_ir_iptw_time$row_names <- rownames(itt_ir_iptw_time)
  itt_ir_iptw_time <- itt_ir_iptw_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, itt_ir_iptw_time, by = 'row_names', suffix = c('', '.iptw'))
  
  ### 2. AT ANALYSES ###
  
  ## A. PREPARING DATA ##
  
  at_d <- d
  
  at_d$time <- pmin(at_d$at_exit_date, at_d$at_event_date, na.rm = TRUE)
  at_d$event_at_time <- ifelse(is.na(at_d$at_event_date), 0, ifelse(at_d$time == at_d$at_event_date, 1, 0))
  at_d$event_counting_time <- as.numeric(difftime(at_d$at_event_date, at_d$entry_date, units = 'days'))
  at_d$censor_counting_time <- as.numeric(difftime(at_d$at_exit_date, at_d$entry_date, units = 'days'))
  at_d$counting_time <- as.numeric(difftime(at_d$time, at_d$entry_date, units = "days"))
  at_d <- at_d[order(at_d$counting_time),]
  
  at_d %<>%
    filter(counting_time>0)
  
  at_d$Tstart <- 0
  at_d <- survSplit(at_d, cut=times_quart, end="counting_time",
                                  start="Tstart", event="event_at_time")
  
  names(at_d)[names(at_d) == 'counting_time'] <- 'Tstop'
  at_d$event_at_tstop <- if_else(is.na(at_d$event_counting_time), 0, 
                                               if_else(at_d$Tstop == at_d$event_counting_time, 1, 0)) 
  
  # retrieve covariate values at quartiles
  for (i in 1:length(comorbidities)) {
    comorb_name <- comorbidities[i]
    at_d$comorb <- if_else(
      at_d$Tstart == 0,
      at_d[[paste(comorb_name, 'base', sep = '_')]],
      if_else(
        at_d$Tstart == times_quart[[1]],
        at_d[[paste(comorb_name, 'Q1', sep = '_')]],
        if_else(
          at_d$Tstart == times_quart[[2]],
          at_d[[paste(comorb_name, 'Q2', sep = '_')]],
          at_d[[paste(comorb_name, 'Q3', sep = '_')]]
        )
      )
    )
    names(at_d)[names(at_d) == 'comorb'] <- comorb_name
  }
  
  # add IPCW corresponding to time interval
  at_d$quartTstart <- ifelse(at_d$Tstart >= 0 & at_d$Tstart < times_quart[1], 0,
                                           ifelse(at_d$Tstart >= times_quart[1] & at_d$Tstart < times_quart[2], times_quart[1],
                                                  ifelse(at_d$Tstart >= times_quart[2] & at_d$Tstart < times_quart[3], times_quart[2], times_quart[3])))
  
  ipcw_weights <- d_long %>% select(boot_id, Tstart, ipcw, sipcw, ipcw_nl, sipcw_nl, ipcw_pooled, sipcw_pooled, ipcw_nl_mod, sipcw_nl_mod)
  
  at_d <- at_d %>%
    dplyr::left_join(ipcw_weights, by = c('boot_id', 'quartTstart' = 'Tstart'), relationship = 'many-to-one')
  
  ## B. INCIDENCE RATES ##
  
  at_d$at_event <- as.numeric(at_d$at_event)
  
  at_d %<>%
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
      
      # IPCW_nl and sIPCW_nl - non-lagged
      ipcw_nl_mod_person_time = person_time * ipcw_nl_mod,
      ipcw_nl_mod_event = event_at_tstop * ipcw_nl_mod,
      
      sipcw_nl_mod_person_time = person_time * sipcw_nl_mod,
      sipcw_nl_mod_event = event_at_tstop * sipcw_nl_mod,
      
      iptw_ipcw_nl_mod_person_time = person_time * iptw * ipcw_nl_mod,
      iptw_ipcw_nl_mod_event = event_at_tstop * iptw * ipcw_nl_mod,
      
      siptw_sipcw_nl_mod_person_time = person_time * siptw * sipcw_nl_mod,
      siptw_sipcw_nl_mod_event = event_at_tstop * siptw * sipcw_nl_mod
    )
  
  
  ## over entire study
  
  # no weights
  at_ir <- at_d %>%
    summarise (
      IR = sum(event_at_tstop) / (sum(person_time)/365), 
      IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
      IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365, 
      IRR = IR_comp / IR_ref, 
      IRD = IR_comp - IR_ref,
      IR_stab = 0, 
      IR_ref_stab = 0,
      IR_comp_stab = 0,
      IRR_stab = 0,
      IRD_stab = 0
    )
  
  at_ir <- data.frame(t(at_ir))
  colnames(at_ir) <- 'at.overall'
  at_ir$row_names <- rownames(at_ir)
  at_ir <- at_ir %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir, by = 'row_names')
  
  # iptw and siptw weights
  at_ir_iptw <- at_d %>%
    summarise (
      IR = sum(iptw_event) / sum(iptw_person_time)/365, 
      IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_event) / sum(siptw_person_time), 
      IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_iptw <- data.frame(t(at_ir_iptw))
  colnames(at_ir_iptw) <- 'at.overall.iptw'
  at_ir_iptw$row_names <- rownames(at_ir_iptw)
  at_ir_iptw <- at_ir_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_iptw, by = 'row_names')
  
  # lagged weights
  at_ir_ipcw <- at_d %>%
    summarise (
      IR = sum(ipcw_event) / sum(ipcw_person_time)/365, 
      IR_ref = sum(ipcw_event[trt_dummy == 0])/sum(ipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_event[trt_dummy == 1])/sum(ipcw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_event) / sum(sipcw_person_time), 
      IR_ref_stab = sum(sipcw_event[trt_dummy == 0])/sum(sipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_event[trt_dummy == 1])/sum(sipcw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw <- data.frame(t(at_ir_ipcw))
  colnames(at_ir_ipcw) <- 'at.overall.ipcw'
  at_ir_ipcw$row_names <- rownames(at_ir_ipcw)
  at_ir_ipcw <- at_ir_ipcw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw, by = 'row_names')
  
  at_ir_ipcw_iptw <- at_d %>%
    summarise (
      IR = sum(iptw_ipcw_event) / sum(iptw_ipcw_person_time)/365, 
      IR_ref = sum(iptw_ipcw_event[trt_dummy == 0])/sum(iptw_ipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_event[trt_dummy == 1])/sum(iptw_ipcw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_event) / sum(siptw_sipcw_person_time), 
      IR_ref_stab = sum(siptw_sipcw_event[trt_dummy == 0])/sum(siptw_sipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_event[trt_dummy == 1])/sum(siptw_sipcw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_iptw <- data.frame(t(at_ir_ipcw_iptw))
  colnames(at_ir_ipcw_iptw) <- 'at.overall.ipcw.iptw'
  at_ir_ipcw_iptw$row_names <- rownames(at_ir_ipcw_iptw)
  at_ir_ipcw_iptw <- at_ir_ipcw_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_iptw, by = 'row_names')
  
  # non-lagged weights
  at_ir_ipcw_nl <- at_d %>%
    summarise (
      IR = sum(ipcw_nl_event) / sum(ipcw_nl_person_time)/365, 
      IR_ref = sum(ipcw_nl_event[trt_dummy == 0])/sum(ipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_nl_event[trt_dummy == 1])/sum(ipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_nl_event) / sum(sipcw_nl_person_time), 
      IR_ref_stab = sum(sipcw_nl_event[trt_dummy == 0])/sum(sipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_nl_event[trt_dummy == 1])/sum(sipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl <- data.frame(t(at_ir_ipcw_nl))
  colnames(at_ir_ipcw_nl) <- 'at.overall.ipcw_nl'
  at_ir_ipcw_nl$row_names <- rownames(at_ir_ipcw_nl)
  at_ir_ipcw_nl <- at_ir_ipcw_nl %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl, by = 'row_names')
  
  at_ir_ipcw_nl_iptw <- at_d %>%
    summarise (
      IR = sum(iptw_ipcw_nl_event) / sum(iptw_ipcw_nl_person_time)/365, 
      IR_ref = sum(iptw_ipcw_nl_event[trt_dummy == 0])/sum(iptw_ipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_nl_event[trt_dummy == 1])/sum(iptw_ipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_nl_event) / sum(siptw_sipcw_nl_person_time), 
      IR_ref_stab = sum(siptw_sipcw_nl_event[trt_dummy == 0])/sum(siptw_sipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_nl_event[trt_dummy == 1])/sum(siptw_sipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_iptw <- data.frame(t(at_ir_ipcw_nl_iptw))
  colnames(at_ir_ipcw_nl_iptw) <- 'at.overall.ipcw_nl.iptw'
  at_ir_ipcw_nl_iptw$row_names <- rownames(at_ir_ipcw_nl_iptw)
  at_ir_ipcw_nl_iptw <- at_ir_ipcw_nl_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_iptw, by = 'row_names')
  
  # pooled weights
  at_ir_ipcw_pooled <- at_d %>%
    summarise (
      IR = sum(ipcw_pooled_event) / sum(ipcw_pooled_person_time)/365, 
      IR_ref = sum(ipcw_pooled_event[trt_dummy == 0])/sum(ipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_pooled_event[trt_dummy == 1])/sum(ipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_pooled_event) / sum(sipcw_pooled_person_time), 
      IR_ref_stab = sum(sipcw_pooled_event[trt_dummy == 0])/sum(sipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_pooled_event[trt_dummy == 1])/sum(sipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_pooled <- data.frame(t(at_ir_ipcw_pooled))
  colnames(at_ir_ipcw_pooled) <- 'at.overall.ipcw_pooled'
  at_ir_ipcw_pooled$row_names <- rownames(at_ir_ipcw_pooled)
  at_ir_ipcw_pooled <- at_ir_ipcw_pooled %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_pooled, by = 'row_names')
  
  at_ir_ipcw_pooled_iptw <- at_d %>%
    summarise (
      IR = sum(iptw_ipcw_pooled_event) / sum(iptw_ipcw_pooled_person_time)/365, 
      IR_ref = sum(iptw_ipcw_pooled_event[trt_dummy == 0])/sum(iptw_ipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_pooled_event[trt_dummy == 1])/sum(iptw_ipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_pooled_event) / sum(siptw_sipcw_pooled_person_time), 
      IR_ref_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 0])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 1])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_pooled_iptw <- data.frame(t(at_ir_ipcw_pooled_iptw))
  colnames(at_ir_ipcw_pooled_iptw) <- 'at.overall.ipcw_pooled.iptw'
  at_ir_ipcw_pooled_iptw$row_names <- rownames(at_ir_ipcw_pooled_iptw)
  at_ir_ipcw_pooled_iptw <- at_ir_ipcw_pooled_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_pooled_iptw, by = 'row_names')
  
  # modified non-lagged weigths
  
  # non-lagged weights
  at_ir_ipcw_nl_mod <- at_d %>%
    summarise (
      IR = sum(ipcw_nl_mod_event) / sum(ipcw_nl_mod_person_time)/365, 
      IR_ref = sum(ipcw_nl_mod_event[trt_dummy == 0])/sum(ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_nl_mod_event[trt_dummy == 1])/sum(ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_nl_mod_event) / sum(sipcw_nl_mod_person_time), 
      IR_ref_stab = sum(sipcw_nl_mod_event[trt_dummy == 0])/sum(sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_nl_mod_event[trt_dummy == 1])/sum(sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_mod <- data.frame(t(at_ir_ipcw_nl_mod))
  colnames(at_ir_ipcw_nl_mod) <- 'at.overall.ipcw_nl_mod'
  at_ir_ipcw_nl_mod$row_names <- rownames(at_ir_ipcw_nl_mod)
  at_ir_ipcw_nl_mod <- at_ir_ipcw_nl_mod %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_mod, by = 'row_names')
  
  at_ir_ipcw_nl_mod_iptw <- at_d %>%
    summarise (
      IR = sum(iptw_ipcw_nl_mod_event) / sum(iptw_ipcw_nl_mod_person_time)/365, 
      IR_ref = sum(iptw_ipcw_nl_mod_event[trt_dummy == 0])/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_nl_mod_event[trt_dummy == 1])/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_nl_mod_event) / sum(siptw_sipcw_nl_mod_person_time), 
      IR_ref_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 0])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 1])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_mod_iptw <- data.frame(t(at_ir_ipcw_nl_mod_iptw))
  colnames(at_ir_ipcw_nl_mod_iptw) <- 'at.overall.ipcw_nl_mod.iptw'
  at_ir_ipcw_nl_mod_iptw$row_names <- rownames(at_ir_ipcw_nl_mod_iptw)
  at_ir_ipcw_nl_mod_iptw <- at_ir_ipcw_nl_mod_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_mod_iptw, by = 'row_names')
  
  ## over time
  
  # no weights
  at_ir_time <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(event_at_tstop) / (sum(person_time)/365), 
      IR_ref = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 0])/sum(person_time[trt_dummy == 0]) * 365,
      IR_comp = n_distinct(boot_id[event_at_tstop == 1 & trt_dummy == 1])/sum(person_time[trt_dummy == 1]) * 365, 
      IRR = IR_comp / IR_ref, 
      IRD = IR_comp - IR_ref,
      IR_stab = 0, 
      IR_ref_stab = 0,
      IR_comp_stab = 0,
      IRR_stab = 0,
      IRD_stab = 0
    )
  
  at_ir_time <- data.frame(t(at_ir_time))
  colnames(at_ir_time) <- at_ir_time[1,]
  at_ir_time <- at_ir_time[-1,]
  colnames(at_ir_time) <- paste('at', colnames(at_ir_time), sep = '.')
  at_ir_time$row_names <- rownames(at_ir_time)
  at_ir_time <- at_ir_time %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_time, by = 'row_names')
  
  # iptw and siptw weights
  at_ir_iptw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_event) / sum(iptw_person_time)/365, 
      IR_ref = sum(iptw_event[trt_dummy == 0])/sum(iptw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_event[trt_dummy == 1])/sum(iptw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_event) / sum(siptw_person_time), 
      IR_ref_stab = sum(siptw_event[trt_dummy == 0])/sum(siptw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_event[trt_dummy == 1])/sum(siptw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_iptw <- data.frame(t(at_ir_iptw))
  colnames(at_ir_iptw) <- at_ir_iptw[1,]
  at_ir_iptw <- at_ir_iptw[-1,]
  colnames(at_ir_iptw) <- paste('at', colnames(at_ir_iptw), sep = '.')
  at_ir_iptw$row_names <- rownames(at_ir_iptw)
  at_ir_iptw <- at_ir_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_iptw, by = 'row_names', suffix = c('', '.iptw'))
  
  # lagged weights
  at_ir_ipcw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(ipcw_event) / sum(ipcw_person_time)/365, 
      IR_ref = sum(ipcw_event[trt_dummy == 0])/sum(ipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_event[trt_dummy == 1])/sum(ipcw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_event) / sum(sipcw_person_time), 
      IR_ref_stab = sum(sipcw_event[trt_dummy == 0])/sum(sipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_event[trt_dummy == 1])/sum(sipcw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw <- data.frame(t(at_ir_ipcw))
  colnames(at_ir_ipcw) <- at_ir_ipcw[1,]
  at_ir_ipcw <- at_ir_ipcw[-1,]
  colnames(at_ir_ipcw) <- paste('at', colnames(at_ir_ipcw), sep = '.')
  at_ir_ipcw$row_names <- rownames(at_ir_ipcw)
  at_ir_ipcw <- at_ir_ipcw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw, by = 'row_names', suffix = c('', '.ipcw'))
  
  at_ir_ipcw_iptw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_ipcw_event) / sum(iptw_ipcw_person_time)/365, 
      IR_ref = sum(iptw_ipcw_event[trt_dummy == 0])/sum(iptw_ipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_event[trt_dummy == 1])/sum(iptw_ipcw_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_event) / sum(siptw_sipcw_person_time), 
      IR_ref_stab = sum(siptw_sipcw_event[trt_dummy == 0])/sum(siptw_sipcw_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_event[trt_dummy == 1])/sum(siptw_sipcw_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_iptw <- data.frame(t(at_ir_ipcw_iptw))
  colnames(at_ir_ipcw_iptw) <- at_ir_ipcw_iptw[1,]
  at_ir_ipcw_iptw <- at_ir_ipcw_iptw[-1,]
  colnames(at_ir_ipcw_iptw) <- paste('at', colnames(at_ir_ipcw_iptw), sep = '.')
  at_ir_ipcw_iptw$row_names <- rownames(at_ir_ipcw_iptw)
  at_ir_ipcw_iptw <- at_ir_ipcw_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_iptw, by = 'row_names', suffix = c('', '.iptw.ipcw'))
  
  # non-lagged weights
  at_ir_ipcw_nl <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(ipcw_nl_event) / sum(ipcw_nl_person_time)/365, 
      IR_ref = sum(ipcw_nl_event[trt_dummy == 0])/sum(ipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_nl_event[trt_dummy == 1])/sum(ipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_nl_event) / sum(sipcw_nl_person_time), 
      IR_ref_stab = sum(sipcw_nl_event[trt_dummy == 0])/sum(sipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_nl_event[trt_dummy == 1])/sum(sipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl <- data.frame(t(at_ir_ipcw_nl))
  colnames(at_ir_ipcw_nl) <- at_ir_ipcw_nl[1,]
  at_ir_ipcw_nl <- at_ir_ipcw_nl[-1,]
  colnames(at_ir_ipcw_nl) <- paste('at', colnames(at_ir_ipcw_nl), sep = '.')
  at_ir_ipcw_nl$row_names <- rownames(at_ir_ipcw_nl)
  at_ir_ipcw_nl <- at_ir_ipcw_nl %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl, by = 'row_names', suffix = c('', '.ipcw_nl'))
  
  at_ir_ipcw_nl_iptw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_ipcw_nl_event) / sum(iptw_ipcw_nl_person_time)/365, 
      IR_ref = sum(iptw_ipcw_nl_event[trt_dummy == 0])/sum(iptw_ipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_nl_event[trt_dummy == 1])/sum(iptw_ipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_nl_event) / sum(siptw_sipcw_nl_person_time), 
      IR_ref_stab = sum(siptw_sipcw_nl_event[trt_dummy == 0])/sum(siptw_sipcw_nl_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_nl_event[trt_dummy == 1])/sum(siptw_sipcw_nl_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_iptw <- data.frame(t(at_ir_ipcw_nl_iptw))
  colnames(at_ir_ipcw_nl_iptw) <- at_ir_ipcw_nl_iptw[1,]
  at_ir_ipcw_nl_iptw <- at_ir_ipcw_nl_iptw[-1,]
  colnames(at_ir_ipcw_nl_iptw) <- paste('at', colnames(at_ir_ipcw_nl_iptw), sep = '.')
  at_ir_ipcw_nl_iptw$row_names <- rownames(at_ir_ipcw_nl_iptw)
  at_ir_ipcw_nl_iptw <- at_ir_ipcw_nl_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_iptw, by = 'row_names', suffix = c('', '.iptw.ipcw_nl'))
  
  # pooled weights
  at_ir_ipcw_pooled <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(ipcw_pooled_event) / sum(ipcw_pooled_person_time)/365, 
      IR_ref = sum(ipcw_pooled_event[trt_dummy == 0])/sum(ipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_pooled_event[trt_dummy == 1])/sum(ipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_pooled_event) / sum(sipcw_pooled_person_time), 
      IR_ref_stab = sum(sipcw_pooled_event[trt_dummy == 0])/sum(sipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_pooled_event[trt_dummy == 1])/sum(sipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_pooled <- data.frame(t(at_ir_ipcw_pooled))
  colnames(at_ir_ipcw_pooled) <- at_ir_ipcw_pooled[1,]
  at_ir_ipcw_pooled <- at_ir_ipcw_pooled[-1,]
  colnames(at_ir_ipcw_pooled) <- paste('at', colnames(at_ir_ipcw_pooled), sep = '.')
  at_ir_ipcw_pooled$row_names <- rownames(at_ir_ipcw_pooled)
  at_ir_ipcw_pooled <- at_ir_ipcw_pooled %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_pooled, by = 'row_names', suffix = c('', '.ipcw_pooled'))
  
  at_ir_ipcw_pooled_iptw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_ipcw_pooled_event) / sum(iptw_ipcw_pooled_person_time)/365, 
      IR_ref = sum(iptw_ipcw_pooled_event[trt_dummy == 0])/sum(iptw_ipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_pooled_event[trt_dummy == 1])/sum(iptw_ipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_pooled_event) / sum(siptw_sipcw_pooled_person_time), 
      IR_ref_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 0])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_pooled_event[trt_dummy == 1])/sum(siptw_sipcw_pooled_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_pooled_iptw <- data.frame(t(at_ir_ipcw_pooled_iptw))
  colnames(at_ir_ipcw_pooled_iptw) <- at_ir_ipcw_pooled_iptw[1,]
  at_ir_ipcw_pooled_iptw <- at_ir_ipcw_pooled_iptw[-1,]
  colnames(at_ir_ipcw_pooled_iptw) <- paste('at', colnames(at_ir_ipcw_pooled_iptw), sep = '.')
  at_ir_ipcw_pooled_iptw$row_names <- rownames(at_ir_ipcw_pooled_iptw)
  at_ir_ipcw_pooled_iptw <- at_ir_ipcw_pooled_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_pooled_iptw, by = 'row_names', suffix = c('', '.iptw.ipcw_pooled'))
  
  # modified non-lagged weights
  
  # non-lagged weights
  at_ir_ipcw_nl_mod <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(ipcw_nl_mod_event) / sum(ipcw_nl_mod_person_time)/365, 
      IR_ref = sum(ipcw_nl_mod_event[trt_dummy == 0])/sum(ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(ipcw_nl_mod_event[trt_dummy == 1])/sum(ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(sipcw_nl_mod_event) / sum(sipcw_nl_mod_person_time), 
      IR_ref_stab = sum(sipcw_nl_mod_event[trt_dummy == 0])/sum(sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(sipcw_nl_mod_event[trt_dummy == 1])/sum(sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_mod <- data.frame(t(at_ir_ipcw_nl_mod))
  colnames(at_ir_ipcw_nl_mod) <- at_ir_ipcw_nl_mod[1,]
  at_ir_ipcw_nl_mod <- at_ir_ipcw_nl_mod[-1,]
  colnames(at_ir_ipcw_nl_mod) <- paste('at', colnames(at_ir_ipcw_nl_mod), sep = '.')
  at_ir_ipcw_nl_mod$row_names <- rownames(at_ir_ipcw_nl_mod)
  at_ir_ipcw_nl_mod <- at_ir_ipcw_nl_mod %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_mod, by = 'row_names', suffix = c('', '.ipcw_nl_mod'))
  
  at_ir_ipcw_nl_mod_iptw <- at_d %>%
    group_by(month_year) %>% 
    summarise (
      IR = sum(iptw_ipcw_nl_mod_event) / sum(iptw_ipcw_nl_mod_person_time)/365, 
      IR_ref = sum(iptw_ipcw_nl_mod_event[trt_dummy == 0])/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp = sum(iptw_ipcw_nl_mod_event[trt_dummy == 1])/sum(iptw_ipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR = IR_comp / IR_ref,
      IRD = IR_comp - IR_ref,
      IR_stab = sum(siptw_sipcw_nl_mod_event) / sum(siptw_sipcw_nl_mod_person_time), 
      IR_ref_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 0])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 0]) * 365,
      IR_comp_stab = sum(siptw_sipcw_nl_mod_event[trt_dummy == 1])/sum(siptw_sipcw_nl_mod_person_time[trt_dummy == 1]) * 365,
      IRR_stab = IR_comp_stab / IR_ref_stab,
      IRD_stab = IR_comp_stab - IR_ref_stab
    )
  
  at_ir_ipcw_nl_mod_iptw <- data.frame(t(at_ir_ipcw_nl_mod_iptw))
  colnames(at_ir_ipcw_nl_mod_iptw) <- at_ir_ipcw_nl_mod_iptw[1,]
  at_ir_ipcw_nl_mod_iptw <- at_ir_ipcw_nl_mod_iptw[-1,]
  colnames(at_ir_ipcw_nl_mod_iptw) <- paste('at', colnames(at_ir_ipcw_nl_mod_iptw), sep = '.')
  at_ir_ipcw_nl_mod_iptw$row_names <- rownames(at_ir_ipcw_nl_mod_iptw)
  at_ir_ipcw_nl_mod_iptw <- at_ir_ipcw_nl_mod_iptw %>% relocate(row_names)
  
  ir_results <- left_join(ir_results, at_ir_ipcw_nl_mod_iptw, by = 'row_names', suffix = c('', '.iptw.ipcw_nl_mod'))
  
  
  ### 3. RESULTS ###
  
  ir_results <- ir_results[,-1]
  
  result <- c(
    IR = as.numeric(ir_results[1,]),
    IR_ref = as.numeric(ir_results[2,]),
    IR_comp = as.numeric(ir_results[3,]),
    IRR = as.numeric(ir_results[4,]),
    IRD = as.numeric(ir_results[5,]),
    IR_stab = as.numeric(ir_results[6,]),
    IR_ref_stab = as.numeric(ir_results[7,]),
    IR_comp_stab = as.numeric(ir_results[8,]),
    IRR_stab = as.numeric(ir_results[9,]),
    IRD_stab = as.numeric(ir_results[10,])
  )
  
  save_bootstrap <<- rbind(save_bootstrap, result)
  
  return(result)
  
}

# extract columns names

# names <- names(ir_results)
# endings <- c(
#   rep('.IR', each = 490),
#   rep('.IR_ref', each = 490),
#   rep('.IR_comp', each = 490),
#   rep('.IRR', each = 490),
#   rep('.IRD', each = 490),
#   rep('.IR_stab', each = 490),
#   rep('.IR_ref_stab', each = 490),
#   rep('.IR_comp_stab', each = 490),
#   rep('.IRR_stab', each = 490),
#   rep('.IRD_stab', each = 490)
#   )
# 
# new_names <- rep(names, times = 10) 
# new_names <- paste(new_names, endings, sep = '')
# 
# saveRDS(new_names, file = paste(path_results, 'result_names.rds', sep = '/'))

#### BOOTSTRAP EXECUTION ####

result_names <- readRDS(file = paste(path_results, 'result_names.rds', sep='/'))

save_bootstrap <- data.frame(matrix(nrow = 1, ncol = length(result_names)))
names(save_bootstrap) <- result_names

bootstrap_results <- boot(data = cohort, statistic = bs, R = 1000, parallel="snow")

save_bootstrap <- save_bootstrap[-1,]
saveRDS(save_bootstrap, file = paste(path_results, 'save_bootstrap.rds', sep='/')) 
saveRDS(bootstrap_results, file = paste(path_results, 'bootstrap_results.rds', sep='/')) 

#### BUILD BOOTSTRAPPED CONFIDENCE INTERVALS ####

save_bootstrap <- readRDS(paste(path_results, 'save_bootstrap.rds', sep='/'))

bootstrap_ci <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(bootstrap_ci) <- c('estimate', 'lower_ci', 'upper_ci')

get_stat_ci <- function(stat_name) {
  
  stat_ci <- data.frame(
    estimate = NA,
    lower_ci = NA,
    upper_ci = NA
  )
  
  stat_values <- save_bootstrap[stat_name]
  
  stat_ci$estimate <- quantile(stat_values, prob = c(0.5), na.rm = TRUE)
  stat_ci$lower_ci <- quantile(stat_values, prob = c(0.025), na.rm = TRUE)
  stat_ci$upper_ci <- quantile(stat_values, prob = c(0.975), na.rm = TRUE)
  
  rownames(stat_ci)[1] <- stat_name
  bootstrap_ci <<- rbind(bootstrap_ci, stat_ci)
  
}

stat_names <- names(save_bootstrap)
invisible(lapply(stat_names, get_stat_ci))

bootstrap_ci$variable <- row.names(bootstrap_ci)

bootstrap_ci <- bootstrap_ci %>%
  select(variable, everything())
  
saveRDS(bootstrap_ci, file = paste(path_results, 'bootstrap_ci.rds', sep='/'))

# start: itt. or at. 
# middle: .overall or .y-m (ex: 2019-01)
# middle: .iptw, .iptw_ipcw
# end: .IRR, .IRD, .IRR_stab, .IRD_stab

rows_to_keep <- c(
  ## ITT
  'itt.overall.IRR',
  'itt.overall.IRD',
  'itt.overall.iptw.IRR',
  'itt.overall.iptw.IRD',
  'itt.overall.iptw.IRR_stab',
  'itt.overall.iptw.IRD_stab',
  'at.overall.IRR',
  'at.overall.IRD',
  'at.overall.iptw.IRR',
  'at.overall.iptw.IRD',
  'at.overall.iptw.IRR_stab',
  'at.overall.iptw.IRD_stab',
  
  ## AT
  # lagged weights
  'at.overall.ipcw.IRR',
  'at.overall.ipcw.IRD',
  'at.overall.ipcw.IRR_stab',
  'at.overall.ipcw.IRD_stab',
  'at.overall.ipcw.iptw.IRR',
  'at.overall.ipcw.iptw.IRD',
  'at.overall.ipcw.iptw.IRR_stab',
  'at.overall.ipcw.iptw.IRD_stab',
  
  # non-lagged weights
  'at.overall.ipcw_nl.IRR',
  'at.overall.ipcw_nl.IRD',
  'at.overall.ipcw_nl.IRR_stab',
  'at.overall.ipcw_nl.IRD_stab',
  'at.overall.ipcw_nl.iptw.IRR',
  'at.overall.ipcw_nl.iptw.IRD',
  'at.overall.ipcw_nl.iptw.IRR_stab',
  'at.overall.ipcw_nl.iptw.IRD_stab',
  
  
  # pooled weights
  'at.overall.ipcw_pooled.IRR',
  'at.overall.ipcw_pooled.IRD',
  'at.overall.ipcw_pooled.IRR_stab',
  'at.overall.ipcw_pooled.IRD_stab',
  'at.overall.ipcw_pooled.iptw.IRR',
  'at.overall.ipcw_pooled.iptw.IRD',
  'at.overall.ipcw_pooled.iptw.IRR_stab',
  'at.overall.ipcw_pooled.iptw.IRD_stab',
  
  # modified non-lagged weights
  'at.overall.ipcw_nl_mod.IRR',
  'at.overall.ipcw_nl_mod.IRD',
  'at.overall.ipcw_nl_mod.IRR_stab',
  'at.overall.ipcw_nl_mod.IRD_stab',
  'at.overall.ipcw_nl_mod.iptw.IRR',
  'at.overall.ipcw_nl_mod.iptw.IRD',
  'at.overall.ipcw_nl_mod.iptw.IRR_stab',
  'at.overall.ipcw_nl_mod.iptw.IRD_stab',
)

results <- bootstrap_ci[row.names(bootstrap_ci) %in% rows_to_keep,]

results$variable <- row.names(results)
results <- results %>%
  select(variable, everything())

write_xlsx(results, paste(path_results, 'bootstrap_ci_overall.xlsx', sep ='/'))

#### USING BOOT.CI FUNCTION - NOT WORKING ####

#RDS(save_bootstrap, file = paste(path_results, 'save_bootstrap.rds', sep='/')) 
# 
# 
# bootstrap_ci <- data.frame(matrix(nrow = 0, ncol = 3))
# colnames(bootstrap_ci) <- c('estimate', 'lower_ci', 'upper_ci')
# 
# for (i in 1:length(bootstrap_results$t0)) {
#   stat <- names(save_bootstrap[i])
#   stat
#   
#   stat_ci <- data.frame(estimate = NA,
#                         lower_ci = NA,
#                         upper_ci = NA)
#   
#   conf.int <- boot.ci(bootstrap_results, type = 'perc', index = i)
#   CI <- conf.int$percent
#   
#   stat_ci$estimate <- bootstrap_results$t0[i]
#   stat_ci$lower_ci <- t(CI[, 4])[[1]]
#   stat_ci$upper_ci <- t(CI[, 5])[[1]]
#   
#   bootstrap_ci <- rbind(bootstrap_ci, stat_ci)
#   rownames(bootstrap_ci)[i] <- stat
#   
# }





