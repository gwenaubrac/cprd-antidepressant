## ---------------------------
##
## Program: 7. IPCW
##
## Purpose: Calculate inverse probability of censoring weights at each decile of the censoring
## distribution, with updated values for the covariates. Three different methods are compared:
## 1) Lagged censoring model: probability of being uncensored is predicted based on data from the previous decile (icpw)
## 2) Non-lagged censoring model: probability of being uncensored is predicted based on data from the current decile (ipcw_nl)
## 3) Pooled censoring model: the probability of being uncensored is predicted on entire data based on pooled model
## that includes a term for the decile. 
##
##
## Author: Gwen Aubrac
##
## Date Created: 2024-07-22
##
## ---------------------------
##
## Notes:
##
## 1. Censoring models are fitted separately for each exposure group. 
## 2. A long format of the cohort dataframe ('cohort_long') is created and split by the deciles of censoring,
## with an indicator for whether or not patient was censored at the Tstop for that decile
## and an indicator for whether or not the patient experienced the event at the Tstop for that decile. 
## 3. There were too few counts for some of the covariates (hypocalcemia, hypomagnesemia)
## so these were removed from the IPCW model due to lack of causal contrast when fitting logistic models.
## ---------------------------

# analysis: flex_grace_period, 90_day_grace_period
# male, female, young, old, 2019, 2020, 2021, 2022
# depressed, not_depressed

analysis <- ''

#### LOAD PACKAGES ####

library(dplyr)
library(magrittr)
library(survival)
library(ggplot2)
library(writexl)
library(broom)

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

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))

setwd(path_results) 
ipcw_desc <- "ipcw_desc.txt"
writeLines("IPCW description:", ipcw_desc)
cat(paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n'), file = ipcw_desc, append= TRUE)

covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_main, 'dec_comorb.rds', sep = '/'))

#### PREPARE DATA ####

cohort_long <- cohort
variables <- c(covariates, comorbidities)
variables <- variables[!variables %in% c('hypocalcemia', 'hypomagnesemia')] # too few counts/violate positivity assumption

if (analysis == 'male' | analysis == 'female') {
  variables <- variables[!variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  variables <- variables[!variables %in% c('year', 'month_year')]
} else if (analysis == 'depressed' | analysis == 'not_depressed') {
  variables <- variables[!variables %in% c('depression')]
}

base_variables <- covariates

if (analysis == 'male' | analysis == 'female') {
  base_variables <- base_variables[!base_variables %in% c('sex')]
} else if (analysis == '2019' | analysis == '2020' | analysis == '2021' | analysis == '2022') {
  base_variables <- base_variables[!base_variables %in% c('year', 'month_year')]
}

base_variables

#### EXAMINE PREDICTORS OF CENSORING ####

predictors_switch <- glm(reformulate(base_variables, 'switch'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_switch), paste(path_results, 'predictors_switch.xlsx', sep ='/'))

predictors_disc <- glm(reformulate(base_variables, 'disc'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_disc), paste(path_results, 'predictors_disc.xlsx', sep ='/'))

predictors_cens <- glm(reformulate(base_variables, 'censor'), data = cohort, family = binomial())
write_xlsx(tidy(predictors_cens), paste(path_results, 'predictors_cens.xlsx', sep ='/'))

#### CONVERT TO COUNTING TIME ####

# time: time to censoring (from disc or switch only) or cohort exit (ITT_exit_date, which includes experiencing event)
# counting_time: time in days between cohort entry and censoring or cohort exit, whichever occurred first
# uncensored_at_time: indicator for whether patient was censored at time (1 if uncensored, meaning did not leave cohort due to disc or switch)
# times_dec: deciles of the censoring distribution for any reason (based on distribution of at_exit_date)

cohort_long$time <- cohort_long$at_exit_date
cohort_long$uncensored_at_time <- if_else(is.na(cohort_long$censor_date), 1, # if no censor date, uncensored
                                          if_else(cohort_long$time == cohort_long$censor_date, 0, 1)) # or if time does not correspond to censoring date, uncensored
cohort_long$counting_time <- as.numeric(difftime(cohort_long$time, cohort_long$entry_date, units = 'days'))
cohort_long$censor_counting_time <- as.numeric(difftime(cohort_long$censor_date, cohort_long$entry_date, units = 'days'))
cohort_long <- cohort_long[order(cohort_long$counting_time),]

times_dec <- readRDS(paste(path_cohort, 'times_dec.rds', sep = '/'))
times_dec <- times_dec[1:9] # remove last decile (100%)
times_dec

# note: will run into issues if someone's follow up ends on same day they enter (counting time of 0)
length(which(cohort$entry_date == cohort$switch_date))
length(which(cohort$entry_date == cohort$itt_exit_date))
length(which(cohort_long$counting_time == 0))

cohort_long %<>%
  filter(counting_time>0)
summary(cohort_long$counting_time)

#### SPLIT DATA INTO INTERVALS OF CENSORING TIMES ####

cohort_long$Tstart <- 0
cohort_long <- survSplit(cohort_long, cut = times_dec, end = 'counting_time', start = 'Tstart', event = 'uncensored_at_time')
names(cohort_long)[names(cohort_long) == 'counting_time'] <- 'Tstop' 
cohort_long$uncensored_at_tstop <- if_else(is.na(cohort_long$censor_counting_time), 1, # indicator for being uncensored for given time interval
                                           if_else(cohort_long$Tstop == cohort_long$censor_counting_time, 0, 1)) 

# retrieve covariate values at deciles
if (analysis == 'depressed' | analysis == 'not_depressed') {
  comorbidities <- comorbidities[!comorbidities %in% c('depression')]
}


for (i in 1:length(comorbidities)) {
  comorb_name <- comorbidities[i]
  cohort_long$comorb <- if_else(
    cohort_long$Tstart == 0,
    cohort_long[[paste(comorb_name, 'base', sep = '_')]],
    if_else(
      cohort_long$Tstart == times_dec[[1]],
      cohort_long[[paste(comorb_name, 'd1', sep = '_')]],
      if_else(
        cohort_long$Tstart == times_dec[[2]],
        cohort_long[[paste(comorb_name, 'd2', sep = '_')]],
        if_else(
          cohort_long$Tstart == times_dec[[3]],
          cohort_long[[paste(comorb_name, 'd3', sep = '_')]],
          if_else(
            cohort_long$Tstart == times_dec[[4]],
            cohort_long[[paste(comorb_name, 'd4', sep = '_')]],
            if_else(
              cohort_long$Tstart == times_dec[[5]],
              cohort_long[[paste(comorb_name, 'd5', sep = '_')]],
              if_else(
                cohort_long$Tstart == times_dec[[6]],
                cohort_long[[paste(comorb_name, 'd6', sep = '_')]],
                if_else(
                  cohort_long$Tstart == times_dec[[7]],
                  cohort_long[[paste(comorb_name, 'd7', sep = '_')]],
                  if_else(
                    cohort_long$Tstart == times_dec[[8]],
                    cohort_long[[paste(comorb_name, 'd8', sep = '_')]], cohort_long[[paste(comorb_name, 'd9', sep = '_')]])
                )
              )
            )
          )
        )
      )
    )
  )
  names(cohort_long)[names(cohort_long) == 'comorb'] <- comorb_name
}

# define censoring model

marginal_formula <- as.formula('uncensored_at_tstop ~ 1')
parametric_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*depression'), collapse = " + ")))
pooled_formula <- as.formula(paste('uncensored_at_tstop', "~", paste(c(variables, 'age_group*depression', 'dec'), collapse = " + ")))

#### MODEL THE CENSORING MECHANISM: LAGGED MODEL ####

# Censoring models are fitted separately for each exposure group. 
# The probability of being uncensored for the d0 to d1 is always 1 (when using this study design)
# because we are interested in the probability of remaining uncensored by the start of the interval
# and by definition all patients in the cohort are in the study at t0.

## CENSORING FOR REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

subset_d1 <- subset(cohort_long_ref, Tstart == 0)
subset_d2 <- subset(cohort_long_ref, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_ref, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_ref, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_ref, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_ref, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_ref, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_ref, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_ref, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_ref, Tstart == times_dec[[9]])

# marginal/baseline probability of remaining uncensored by decile
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_ref, family = binomial())
cohort_long_ref$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(marginal_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(marginal_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(marginal_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(marginal_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(marginal_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(marginal_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(marginal_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(marginal_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(marginal_formula, data = subset_d9, family = binomial)

predict_d2 <- predict(uncensored_in_d1, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d2, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d3, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d4, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d5, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d6, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d7, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d8, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d9, subset_d10, type = 'response')

cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens'] <- 1
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)

# calculate cumulative probability of remaining uncensored
cohort_long_ref <- cohort_long_ref %>%
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

# calculate censoring weights
cohort_long_ref$ipcw <- 1/cohort_long_ref$p_uncens_cum
cohort_long_ref$sipcw <- cohort_long_ref$p_uncens_base/cohort_long_ref$p_uncens_cum

summary(cohort_long_ref$ipcw)
summary(cohort_long_ref$sipcw)

# add weights to long cohort
cohort_long$ipcw <- 0
cohort_long$sipcw <- 0
cohort_long[cohort_long$trt_dummy == 0, 'ipcw'] <- cohort_long_ref$ipcw
cohort_long[cohort_long$trt_dummy == 0, 'sipcw'] <- cohort_long_ref$sipcw

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9)

## CENSORING FOR COMPARATOR GROUP ##

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1)

subset_d1 <- subset(cohort_long_comp, Tstart == 0)
subset_d2 <- subset(cohort_long_comp, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_comp, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_comp, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_comp, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_comp, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_comp, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_comp, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_comp, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_comp, Tstart == times_dec[[9]])

# marginal/baseline probability of remaining uncensored by decile
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_comp, family = binomial())
cohort_long_comp$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(marginal_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(marginal_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(marginal_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(marginal_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(marginal_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(marginal_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(marginal_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(marginal_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(marginal_formula, data = subset_d9, family = binomial)

predict_d2 <- predict(uncensored_in_d1, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d2, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d3, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d4, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d5, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d6, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d7, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d8, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d9, subset_d10, type = 'response')

cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens'] <- 1
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)

# calculate cumulative probability of remaining uncensored
cohort_long_comp <- cohort_long_comp %>% 
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

# calculate censoring weights
cohort_long_comp$ipcw <- 1/cohort_long_comp$p_uncens_cum
cohort_long_comp$sipcw <- cohort_long_comp$p_uncens_base/cohort_long_comp$p_uncens_cum

summary(cohort_long_comp$ipcw)
summary(cohort_long_comp$sipcw)

# add weights to long cohort
cohort_long[cohort_long$trt_dummy == 1, 'ipcw'] <- cohort_long_comp$ipcw
cohort_long[cohort_long$trt_dummy == 1, 'sipcw'] <- cohort_long_comp$sipcw

summary(cohort_long$ipcw)
summary(cohort_long$sipcw)

cat(paste('Unstabilized lagged IPCW weights: \n min:', min(cohort_long$ipcw), '\n max:', max(cohort_long$ipcw), '\n mean:', mean(cohort_long$ipcw), '\n sd:', sd(cohort_long$ipcw), '\n'), file = ipcw_desc, append = TRUE)
cat(paste('Stabilized lagged IPCW weights: \n min:', min(cohort_long$sipcw), '\n max:', max(cohort_long$sipcw), '\n mean:', mean(cohort_long$sipcw), '\n sd:', sd(cohort_long$sipcw), '\n'), file = ipcw_desc, append = TRUE)

# weights distribution by decile
summary(cohort_long[cohort_long$Tstart==0,]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[1]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[2]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[3]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[4]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[5]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[6]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[7]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[8]],]$ipcw)
summary(cohort_long[cohort_long$Tstart==times_dec[[9]],]$ipcw)

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9, uncensored_in_d10)

#### MODEL THE CENSORING MECHANISM: NON-LAGGED MODEL ####

## CENSORING FOR REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

subset_d1 <- subset(cohort_long_ref, Tstart == 0)
subset_d2 <- subset(cohort_long_ref, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_ref, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_ref, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_ref, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_ref, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_ref, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_ref, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_ref, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_ref, Tstart == times_dec[[9]])

# # marginal/baseline probability of remaining uncensored
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_ref, family = binomial())
cohort_long_ref$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(parametric_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(parametric_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(parametric_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(parametric_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(parametric_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(parametric_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(parametric_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(parametric_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(parametric_formula, data = subset_d9, family = binomial)
uncensored_in_d10 <- glm(parametric_formula, data = subset_d10, family = binomial)

predict_d1 <- predict(uncensored_in_d1, subset_d1, type = 'response')
predict_d2 <- predict(uncensored_in_d2, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d3, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d4, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d5, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d6, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d7, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d8, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d9, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d10, subset_d10, type = 'response')

cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens'] <- predict_d1
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d1, predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)

# calculate cumulative probability of remaining uncensored
cohort_long_ref <- cohort_long_ref %>%
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

# calculate censoring weights
cohort_long_ref$ipcw <- 1/cohort_long_ref$p_uncens_cum
cohort_long_ref$sipcw <- cohort_long_ref$p_uncens_base/cohort_long_ref$p_uncens_cum

summary(cohort_long_ref$ipcw)
summary(cohort_long_ref$sipcw)

# add weights to long cohort
cohort_long$ipcw_nl <- 0
cohort_long$sipcw_nl <- 0
cohort_long[cohort_long$trt_dummy == 0, 'ipcw_nl'] <- cohort_long_ref$ipcw
cohort_long[cohort_long$trt_dummy == 0, 'sipcw_nl'] <- cohort_long_ref$sipcw

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9, uncensored_in_d10)

## CENSORING FOR COMPARATOR GROUP ##

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1)

subset_d1 <- subset(cohort_long_comp, Tstart == 0)
subset_d2 <- subset(cohort_long_comp, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_comp, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_comp, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_comp, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_comp, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_comp, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_comp, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_comp, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_comp, Tstart == times_dec[[9]])

# marginal/baseline probability of remaining uncensored
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_comp, family = binomial())
cohort_long_comp$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(parametric_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(parametric_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(parametric_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(parametric_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(parametric_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(parametric_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(parametric_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(parametric_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(parametric_formula, data = subset_d9, family = binomial)
uncensored_in_d10 <- glm(parametric_formula, data = subset_d10, family = binomial)

predict_d1 <- predict(uncensored_in_d1, subset_d1, type = 'response')
predict_d2 <- predict(uncensored_in_d2, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d3, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d4, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d5, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d6, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d7, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d8, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d9, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d10, subset_d10, type = 'response')

cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens'] <- predict_d1
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d1, predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)

# calculate cumulative probability of remaining uncensored
cohort_long_comp <- cohort_long_comp %>% 
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

# calculate censoring weights
cohort_long_comp$ipcw <- 1/cohort_long_comp$p_uncens_cum
cohort_long_comp$sipcw <- cohort_long_comp$p_uncens_base/cohort_long_comp$p_uncens_cum

summary(cohort_long_comp$ipcw)
summary(cohort_long_comp$sipcw)

# add weights to long cohort
cohort_long[cohort_long$trt_dummy == 1, 'ipcw_nl'] <- cohort_long_comp$ipcw
cohort_long[cohort_long$trt_dummy == 1, 'sipcw_nl'] <- cohort_long_comp$sipcw

summary(cohort_long$ipcw_nl)
summary(cohort_long$sipcw_nl)

cat(paste('Unstabilized non-lagged IPCW weights: \n min:', min(cohort_long$ipcw_nl), '\n max:', max(cohort_long$ipcw_nl), '\n mean:', mean(cohort_long$ipcw_nl), '\n sd:', sd(cohort_long$ipcw_nl), '\n'), file = ipcw_desc, append = TRUE)

# truncate extreme stabilized weights to 50
length(which(cohort_long$sipcw_nl > 50))
cat(paste('Non-lagged IPCW weights truncated to 50:', length(which(cohort_long$sipcw_nl > 50)), '\n'), file = ipcw_desc, append = TRUE)

cohort_long %<>%
  mutate(sipcw_nl = if_else(sipcw_nl > 50, 50, sipcw_nl))

cat(paste('Stabilized non-lagged IPCW weights: \n min:', min(cohort_long$sipcw_nl), '\n max:', max(cohort_long$sipcw_nl), '\n mean:', mean(cohort_long$sipcw_nl), '\n sd:', sd(cohort_long$sipcw_nl), '\n'), file = ipcw_desc, append = TRUE)


# weights distribution by decile
summary(cohort_long[cohort_long$Tstart==0,]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[1]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[2]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[3]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[4]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[5]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[6]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[7]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[8]],]$ipcw_nl)
summary(cohort_long[cohort_long$Tstart==times_dec[[9]],]$ipcw_nl)

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9, uncensored_in_d10)

#### MODEL THE CENSORING MECHANISM: POOLED MODEL ####

## REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0) %>% 
  group_by(id) %>% 
  mutate (dec = as.factor(row_number()))

# marginal/baseline probability of remaining uncensored
baseline_uncensored <- glm(marginal_formula, data = cohort_long_ref, family = binomial())
summary(baseline_uncensored)
cohort_long_ref$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional probability of remaining uncensored
uncensored_pooled_ref <- glm(pooled_formula, data = cohort_long_ref, family = binomial) 
summary(uncensored_pooled_ref)

cohort_long_ref$p_uncens_pooled <- predict(uncensored_pooled_ref, type = 'response')

# cumulative probability 
cohort_long_ref <- cohort_long_ref %>%  
  group_by(id) %>% 
  mutate (p_uncens_cum_pooled = cumprod(p_uncens_pooled)) 

# calculate censoring weights
cohort_long_ref$ipcw <- 1/cohort_long_ref$p_uncens_cum_pooled
cohort_long_ref$sipcw <- cohort_long_ref$p_uncens_base/cohort_long_ref$p_uncens_cum_pooled

summary(cohort_long_ref$ipcw)
summary(cohort_long_ref$sipcw)

# add weights to long cohort
cohort_long$ipcw_pooled <- 0
cohort_long$sipcw_pooled <- 0

cohort_long[cohort_long$trt_dummy == 0, 'ipcw_pooled'] <- cohort_long_ref$ipcw
cohort_long[cohort_long$trt_dummy == 0, 'sipcw_pooled'] <- cohort_long_ref$sipcw

rm(baseline_uncensored, uncensored_pooled_ref, cohort_long_ref)

## COMPARATOR GROUP ##

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1) %>% 
  group_by(id) %>% 
  mutate (dec = as.factor(row_number()))

# marginal/baseline probability of remaining uncensored
baseline_uncensored <- glm(marginal_formula, data = cohort_long_comp, family = binomial())
cohort_long_comp$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional probability of remaining uncensored
uncensored_pooled_comp <- glm(pooled_formula, data = cohort_long_comp, family = binomial)
cohort_long_comp$p_uncens_pooled <- predict(uncensored_pooled_comp, type = 'response')

# cumulative probability
cohort_long_comp %<>% 
  group_by(id) %>% 
  mutate (p_uncens_cum_pooled = cumprod(p_uncens_pooled)) 

# calculate censoring weights
cohort_long_comp$ipcw <- 1/cohort_long_comp$p_uncens_cum_pooled
cohort_long_comp$sipcw <- cohort_long_comp$p_uncens_base/cohort_long_comp$p_uncens_cum_pooled

summary(cohort_long_comp$ipcw)
summary(cohort_long_comp$sipcw)

# add weights to long cohort
cohort_long[cohort_long$trt_dummy == 1, 'ipcw_pooled'] <- cohort_long_comp$ipcw
cohort_long[cohort_long$trt_dummy == 1, 'sipcw_pooled'] <- cohort_long_comp$sipcw

rm(baseline_uncensored, uncensored_pooled_comp, cohort_long_comp)

summary(cohort_long$ipcw_pooled)
summary(cohort_long$sipcw_pooled)

cat(paste('Unstabilized pooled IPCW weights: \n min:', min(cohort_long$ipcw_pooled), '\n max:', max(cohort_long$ipcw_pooled), '\n mean:', mean(cohort_long$ipcw_pooled), '\n sd:', sd(cohort_long$ipcw_pooled), '\n'), file = ipcw_desc, append = TRUE)

# truncate extreme stabilized weights to 50
length(which(cohort_long$sipcw_pooled > 50))
cat(paste('Pooled IPCW weights truncated to 50:', length(which(cohort_long$sipcw_pooled > 50)), '\n'), file = ipcw_desc, append = TRUE)

cohort_long %<>%
  mutate(sipcw_pooled = if_else(sipcw_pooled > 50, 50, sipcw_pooled))

cat(paste('Stabilized pooled IPCW weights: \n min:', min(cohort_long$sipcw_pooled), '\n max:', max(cohort_long$sipcw_pooled), '\n mean:', mean(cohort_long$sipcw_pooled), '\n sd:', sd(cohort_long$sipcw_pooled), '\n'), file = ipcw_desc, append = TRUE)

# weights distribution by decile
summary(cohort_long[cohort_long$Tstart==0,]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[1]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[2]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[3]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[4]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[5]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[6]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[7]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[8]],]$ipcw_pooled)
summary(cohort_long[cohort_long$Tstart==times_dec[[9]],]$ipcw_pooled)



#### MODEL THE CENSORING MECHANISM: NON-LAGGED MODEL MODIFIED Q1 ####

## CENSORING FOR REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

subset_d1 <- subset(cohort_long_ref, Tstart == 0)
subset_d2 <- subset(cohort_long_ref, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_ref, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_ref, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_ref, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_ref, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_ref, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_ref, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_ref, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_ref, Tstart == times_dec[[9]])

# marginal/baseline probability of remaining uncensored
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_ref, family = binomial())
cohort_long_ref$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(parametric_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(parametric_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(parametric_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(parametric_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(parametric_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(parametric_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(parametric_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(parametric_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(parametric_formula, data = subset_d9, family = binomial)
uncensored_in_d10 <- glm(parametric_formula, data = subset_d10, family = binomial)

predict_d1 <- predict(uncensored_in_d1, subset_d1, type = 'response')
predict_d2 <- predict(uncensored_in_d2, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d3, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d4, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d5, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d6, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d7, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d8, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d9, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d10, subset_d10, type = 'response')

cohort_long_ref[cohort_long_ref$Tstart==0, 'p_uncens'] <- predict_d1
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_ref[cohort_long_ref$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d1, predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)


# calculate cumulative probability of remaining uncensored
cohort_long_ref <- cohort_long_ref %>%
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

# modify weight for q1 to be 1 to account for lots of censoring happening at end of q1
cohort_long_ref <- cohort_long_ref %>%
  group_by(id) %>% 
  mutate (p_uncens_cum = if_else(row_number() == 1, 1, p_uncens_cum))

# calculate censoring weights
cohort_long_ref$ipcw <- 1/cohort_long_ref$p_uncens_cum
cohort_long_ref$sipcw <- cohort_long_ref$p_uncens_base/cohort_long_ref$p_uncens_cum

summary(cohort_long_ref$ipcw)
summary(cohort_long_ref$sipcw)

# add weights to long cohort
cohort_long$ipcw_nl_mod <- 0
cohort_long$sipcw_nl_mod <- 0
cohort_long[cohort_long$trt_dummy == 0, 'ipcw_nl_mod'] <- cohort_long_ref$ipcw
cohort_long[cohort_long$trt_dummy == 0, 'sipcw_nl_mod'] <- cohort_long_ref$sipcw

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9, uncensored_in_d10)

## CENSORING FOR COMPARATOR GROUP ##

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1)

subset_d1 <- subset(cohort_long_comp, Tstart == 0)
subset_d2 <- subset(cohort_long_comp, Tstart == times_dec[[1]])
subset_d3 <- subset(cohort_long_comp, Tstart == times_dec[[2]])
subset_d4 <- subset(cohort_long_comp, Tstart == times_dec[[3]])
subset_d5 <- subset(cohort_long_comp, Tstart == times_dec[[4]])
subset_d6 <- subset(cohort_long_comp, Tstart == times_dec[[5]])
subset_d7 <- subset(cohort_long_comp, Tstart == times_dec[[6]])
subset_d8 <- subset(cohort_long_comp, Tstart == times_dec[[7]])
subset_d9 <- subset(cohort_long_comp, Tstart == times_dec[[8]])
subset_d10 <- subset(cohort_long_comp, Tstart == times_dec[[9]])

# marginal/baseline probability of remaining uncensored
# baseline_uncensored_d1 <- glm(marginal_formula, data = subset_d1, family = binomial())
# baseline_uncensored_d2 <- glm(marginal_formula, data = subset_d2, family = binomial())
# baseline_uncensored_d3 <- glm(marginal_formula, data = subset_d3, family = binomial())
# baseline_uncensored_d4 <- glm(marginal_formula, data = subset_d4, family = binomial())
# baseline_uncensored_d5 <- glm(marginal_formula, data = subset_d5, family = binomial())
# baseline_uncensored_d6 <- glm(marginal_formula, data = subset_d6, family = binomial())
# baseline_uncensored_d7 <- glm(marginal_formula, data = subset_d7, family = binomial())
# baseline_uncensored_d8 <- glm(marginal_formula, data = subset_d8, family = binomial())
# baseline_uncensored_d9 <- glm(marginal_formula, data = subset_d9, family = binomial())
# baseline_uncensored_d10 <- glm(marginal_formula, data = subset_d10, family = binomial())
# 
# cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens_base'] <- predict(baseline_uncensored_d1, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens_base'] <- predict(baseline_uncensored_d2, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens_base'] <- predict(baseline_uncensored_d3, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens_base'] <- predict(baseline_uncensored_d4, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens_base'] <- predict(baseline_uncensored_d5, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens_base'] <- predict(baseline_uncensored_d6, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens_base'] <- predict(baseline_uncensored_d7, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens_base'] <- predict(baseline_uncensored_d8, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens_base'] <- predict(baseline_uncensored_d9, type = 'response')
# cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens_base'] <- predict(baseline_uncensored_d10, type = 'response')

# marginal probability of remaining uncensored overall
baseline_uncensored <- glm(marginal_formula, data = cohort_long_comp, family = binomial())
cohort_long_comp$p_uncens_base <- predict(baseline_uncensored, type = 'response')

# conditional/parametric probability of remaining uncensored at start of time interval
uncensored_in_d1 <- glm(parametric_formula, data = subset_d1, family = binomial)
uncensored_in_d2 <- glm(parametric_formula, data = subset_d2, family = binomial)
uncensored_in_d3 <- glm(parametric_formula, data = subset_d3, family = binomial)
uncensored_in_d4 <- glm(parametric_formula, data = subset_d4, family = binomial)
uncensored_in_d5 <- glm(parametric_formula, data = subset_d5, family = binomial)
uncensored_in_d6 <- glm(parametric_formula, data = subset_d6, family = binomial)
uncensored_in_d7 <- glm(parametric_formula, data = subset_d7, family = binomial)
uncensored_in_d8 <- glm(parametric_formula, data = subset_d8, family = binomial)
uncensored_in_d9 <- glm(parametric_formula, data = subset_d9, family = binomial)
uncensored_in_d10 <- glm(parametric_formula, data = subset_d10, family = binomial)

predict_d1 <- predict(uncensored_in_d1, subset_d1, type = 'response')
predict_d2 <- predict(uncensored_in_d2, subset_d2, type = 'response')
predict_d3 <- predict(uncensored_in_d3, subset_d3, type = 'response')
predict_d4 <- predict(uncensored_in_d4, subset_d4, type = 'response')
predict_d5 <- predict(uncensored_in_d5, subset_d5, type = 'response')
predict_d6 <- predict(uncensored_in_d6, subset_d6, type = 'response')
predict_d7 <- predict(uncensored_in_d7, subset_d7, type = 'response')
predict_d8 <- predict(uncensored_in_d8, subset_d8, type = 'response')
predict_d9 <- predict(uncensored_in_d9, subset_d9, type = 'response')
predict_d10 <- predict(uncensored_in_d10, subset_d10, type = 'response')

cohort_long_comp[cohort_long_comp$Tstart==0, 'p_uncens'] <- predict_d1
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[1]], 'p_uncens'] <- predict_d2
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[2]], 'p_uncens'] <- predict_d3
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[3]], 'p_uncens'] <- predict_d4
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[4]], 'p_uncens'] <- predict_d5
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[5]], 'p_uncens'] <- predict_d6
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[6]], 'p_uncens'] <- predict_d7
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[7]], 'p_uncens'] <- predict_d8
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[8]], 'p_uncens'] <- predict_d9
cohort_long_comp[cohort_long_comp$Tstart==times_dec[[9]], 'p_uncens'] <- predict_d10

rm(subset_d1, subset_d2, subset_d3, subset_d4, subset_d5, subset_d6, subset_d7, subset_d8, subset_d9, subset_d10)
rm(predict_d1, predict_d2, predict_d3, predict_d4, predict_d5, predict_d6, predict_d7, predict_d8, predict_d9, predict_d10)

# calculate cumulative probability of remaining uncensored
cohort_long_comp <- cohort_long_comp %>% 
  group_by(id) %>% 
  mutate (p_uncens_cum = cumprod(p_uncens)) 

cohort_long_comp <- cohort_long_comp %>%
  group_by(id) %>% 
  mutate (p_uncens_cum = if_else(row_number() == 1, 1, p_uncens_cum))

# calculate censoring weights
cohort_long_comp$ipcw <- 1/cohort_long_comp$p_uncens_cum
cohort_long_comp$sipcw <- cohort_long_comp$p_uncens_base/cohort_long_comp$p_uncens_cum

summary(cohort_long_comp$ipcw)
summary(cohort_long_comp$sipcw)

# add weights to long cohort
cohort_long[cohort_long$trt_dummy == 1, 'ipcw_nl_mod'] <- cohort_long_comp$ipcw
cohort_long[cohort_long$trt_dummy == 1, 'sipcw_nl_mod'] <- cohort_long_comp$sipcw

summary(cohort_long$ipcw_nl_mod)
summary(cohort_long$sipcw_nl_mod)

cat(paste('Unstabilized modified non-lagged IPCW weights: \n min:', min(cohort_long$ipcw_nl_mod), '\n max:', max(cohort_long$ipcw_nl_mod), '\n mean:', mean(cohort_long$ipcw_nl_mod), '\n sd:', sd(cohort_long$ipcw_nl_mod), '\n'), file = ipcw_desc, append = TRUE)

# truncate extreme stabilized weights to 50
length(which(cohort_long$sipcw_nl_mod > 50))
cat(paste('Modified IPCW weights truncated to 50:', length(which(cohort_long$sipcw_nl_mod > 50)), '\n'), file = ipcw_desc, append = TRUE)

cohort_long %<>%
  mutate(sipcw_nl_mod = if_else(sipcw_nl_mod > 50, 50, sipcw_nl_mod))

cat(paste('Stabilized modified non-lagged IPCW weights: \n min:', min(cohort_long$sipcw_nl_mod), '\n max:', max(cohort_long$sipcw_nl_mod), '\n mean:', mean(cohort_long$sipcw_nl_mod), '\n sd:', sd(cohort_long$sipcw_nl_mod), '\n'), file = ipcw_desc, append = TRUE)

# weights distribution by decile
summary(cohort_long[cohort_long$Tstart==0,]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[1]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[2]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[3]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[4]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[5]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[6]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[7]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[8]],]$ipcw_nl_mod)
summary(cohort_long[cohort_long$Tstart==times_dec[[9]],]$ipcw_nl_mod)

rm(uncensored_in_d1, uncensored_in_d2, uncensored_in_d3, uncensored_in_d4, uncensored_in_d5, uncensored_in_d6,
   uncensored_in_d7, uncensored_in_d8, uncensored_in_d9, uncensored_in_d10)
rm(variables, times_dec)

saveRDS(cohort_long, file = paste(path_cohort, 'antidepressant_cohort_ipcw.rds', sep='/'))

