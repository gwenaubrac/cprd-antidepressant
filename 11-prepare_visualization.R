## ---------------------------
##
## Program: 10. Prepare visualization
##
## Purpose: Prepare the data for visualization. 
##
## Author: Gwen Aubrac
##
## Date Created: 2014-07-22
##
## ---------------------------
##
## Notes: 
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

path_vis <- "Z:/EPI/Protocol 24_004042/Gwen/results/visualization"

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))

covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
quart_comorb <- readRDS(file = paste(path_main, 'quart_comorb.rds', sep = '/'))

itt_person_time <- readRDS(file = paste(path_cohort, 'itt_person_time.rds', sep = '/'))
itt_events <- readRDS(file = paste(path_cohort, 'itt_events.rds', sep = '/'))

at_person_time <- readRDS(file = paste(path_cohort, 'at_person_time.rds', sep = '/'))
at_events <- readRDS(file = paste(path_cohort, 'at_events.rds', sep = '/'))

cohort_vis <- merge(cohort, itt_person_time, by.x = 'id', by.y = 'id', all.x = TRUE)
cohort_vis <- merge(cohort_vis, at_person_time, by.x= 'id', by.y = 'id', all.x = TRUE)

cohort_vis$region <- 'CPRD'

#### All X BY YEAR ####

allxbyyear <- cohort_vis %>% 
  group_by(year, trt) %>% 
  summarize(
    num_events = sum(itt_event)
  )

total_events <- sum(cohort_vis$itt_event)

allxbyyear %<>%
  mutate(prop_overall = num_events / total_events) %>%
  group_by(year) %>%
  mutate(prop_by_year = num_events / sum(num_events)) %>%
  ungroup() %>%
  group_by(trt) %>%
  mutate(prop_by_trt = num_events / sum(num_events))

allxbyyear <- cbind(region = 'CPRD', allxbyyear, row.names = NULL)

saveRDS(allxbyyear, file = paste(path_vis, 'allxbyyear.xlsx', sep = '/'))

#### ALL Y BY YEAR ####

cohort_vis %<>%
  filter(at_follow_up>0)

allybyyear <- cohort_vis %>% 
  group_by(year) %>% 
  summarize (
    itt_num_events = sum(itt_event),
    itt_person_days = sum(itt_pt),
    itt_person_years = itt_person_days / 365,
    itt_IR_p100y = itt_num_events/itt_person_years * 100,
    at_num_events = sum(at_event),
    at_person_days = sum(at_pt),
    at_person_years = at_person_days / 365,
    at_IR_p100y = at_num_events/at_person_years * 100
  )

allybyyear <- cbind(region = 'CPRD', allybyyear, row.names = NULL)

saveRDS(allybyyear, file = paste(path_vis, 'allybyyear.xlsx', sep = '/'))

# read all excel with IR for at
test <- read_excel(paste(path_results,'at_ir_ipcw_iptw_time.xlsx', sep = '/')) %>% 
  mutate(at_IR_iptw_ipcwlagged = 100*IR,
         at_IR_siptw_sipcwlagged = 100*IR_stab) %>% 
  select(month_year, at_IR_iptw_ipcwlagged, at_IR_siptw_sipcwlagged)

allybymonth <- cohort_vis %>% 
  group_by(month_year) %>% 
  summarize (
    itt_num_events = sum(itt_event),
    itt_person_days = sum(itt_pt),
    itt_person_years = itt_person_days / 365,
    itt_IR_p100y = itt_num_events/itt_person_years * 100,
    at_num_events = sum(at_event),
    at_person_days = sum(at_pt),
    at_person_years = at_person_days / 365,
    at_IR_p100y = at_num_events/at_person_years * 100
  )

# merge with results
allybymonth <- merge(allybymonth, test, by = 'month_year', all.x = TRUE)

saveRDS(allybymonth, file = paste(path_vis, 'allybymonth.xlsx', sep = '/'))

#### COVARIATES BY REGION ####

covsbyreg <- cohort_vis %>%
  summarise(across(all_of(base_comorb), ~ sum(. == 1, na.rm = TRUE)))

covsbyreg <- covsbyreg %>% 
  pivot_longer(cols = -c(), names_to = "comorbidity", values_to = "count")

covsbyreg <- cbind(region = 'CPRD', covsbyreg, row.names = NULL)

saveRDS(covscoef, file = paste(path_vis, 'covsbyreg.xlsx', sep = '/'))

#### COVARIATE COEFFICIENTS ####

base_variables <- c(covariates, base_comorb) 
grouping_var <- 'region' 
base_variables <- base_variables[!base_variables %in% c('year', 'month_year', 'age_at_entry')]
base_variables <- c(base_variables, 'age_group')

trt_col_name <- 'trt_dummy'

covscoef <- data.frame(group = c(unique(cohort_vis[,grouping_var]))) %>% 
  arrange(group)

for (i in 1:length(base_variables)) {
  var <- base_variables[i]
  
  cols <- ((colnames(cohort_vis) == grouping_var) |
             (colnames(cohort_vis) == trt_col_name) |
             (colnames(cohort_vis) == var))
  
  # extract variable info from cohort_vis
  var_df <- cohort_vis[cols]
  var_i <- colnames(var_df) == var
  
  # check variable type
  if (is.character(var_df[,var_i])) {
    var_df[,var_i] <- as.factor(var_df[,var_i])
  }
  
  if (is.factor(var_df[,var_i]) & nlevels(var_df[,var_i])>2) { # if factor with >2 levels, create dummies
    var_df <- dummy_cols(var_df, select_columns = var, remove_first_dummy = TRUE)
    cols_keep <- colnames(var_df) != var
    var_df <- var_df[cols_keep]
  }
  
  names(var_df)[colnames(var_df) == trt_col_name] <- 'trt'
  names(var_df)[colnames(var_df) == grouping_var] <- 'group'
  
  # iterate through variables (may be multiple if var was a factor)
  for (j in 1:length(var_df)) {
    if (colnames(var_df)[j] == 'trt' | colnames(var_df)[j] == 'group') {
      next
    }
    
    name <- colnames(var_df[j])
    var_sum <- var_df 
    names(var_sum)[j] <- 'var'
    
    covscoef$var <- NA
    
    # calculate association with exposure through logistic regression
    for (i in 1:nrow(covscoef)) {
      group_i <- covscoef[i,1]
      covscoef[covscoef$'group' == group_i, 'var'] <- glm(trt ~ var, data = subset(var_sum, group == group_i), family = 'binomial')$coefficients[[2]]
    }
    
    names(covscoef)[names(covscoef) == 'var'] <- name
  }
}

rm(var_df, var_sum, cols, group_i, i, j, name, var, var_i)

covscoef <- covscoef[,-1] %>%
  pivot_longer(cols = everything(), names_to = "covariate", values_to = "exp(b)")

covscoef <- cbind(region = 'CPRD', covscoef, row.names = NULL)

saveRDS(covscoef, file = paste(path_vis, 'covscoef.xlsx', sep = '/'))

#### MODEL RESULTS ####

# stabilized weights

itt_model <- readRDS(paste(path_results, 'cox_no_weight.rds', sep = '/'))
itt_siptw_model <- readRDS(paste(path_results, 'cox_siptw.rds', sep = '/'))
at_model <- readRDS(paste(path_results, 'cox_at.rds', sep = '/'))
at_siptw_model <- readRDS(paste(path_results, 'cox_at_siptw.rds', sep = '/'))
at_sipcw_model <- readRDS(paste(path_results, 'cox_stab_lagged.rds', sep = '/'))
at_siptw_sipcw_model <- readRDS(paste(path_results, 'cox_stab_lagged_iptw.rds', sep = '/'))
at_sipcw_nl_model <- readRDS(paste(path_results, 'cox_nl_stab.rds', sep = '/'))
at_sipcw_nl_siptw_model <- readRDS(paste(path_results, 'cox_nl_stab_iptw.rds', sep = '/'))
at_sipcw_pooled_model <- readRDS(paste(path_results, 'cox_pooled_stab.rds', sep = '/'))
at_sipcw_pooled_siptw_model <- readRDS(paste(path_results, 'cox_pooled_stab_iptw.rds', sep = '/'))
at_sipcw_nl_mod_model <- readRDS(paste(path_results, 'cox_nl_mod_stab.rds', sep = '/'))
at_sipcw_nl_mod_siptw_model <- readRDS(paste(path_results, 'cox_nl_mod_stab_iptw.rds', sep = '/'))

result_chart_stab <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(result_chart_stab) <- c('estimate', 'lower_ci', 'upper_ci')
rownames(result_chart_stab) <- c('ITT', 'ITT (sIPTW)', 'AT', 'AT (sIPTW)', 'AT (lagged sIPCW)', 'AT (lagged sIPCW + sIPTW)',
                            'AT (non-lagged sIPCW)', 'AT (non-lagged sIPCW + sIPTW)', 'AT (pooled sIPCW)',
                            'AT (pooled sIPCW + sIPTW)', 'AT (mod non-lagged sIPCW)', 'AT (mod non-lagged sIPCW + sIPTW)')

result_chart_stab[1, 'estimate'] <- exp(itt_model$coef)
result_chart_stab[1, 'lower_ci'] <- exp(confint(itt_model))[1]
result_chart_stab[1, 'upper_ci'] <- exp(confint(itt_model))[2]

result_chart_stab[2, 'estimate'] <- exp(itt_siptw_model$coef)
result_chart_stab[2, 'lower_ci'] <- exp(confint(itt_siptw_model))[1]
result_chart_stab[2, 'upper_ci'] <- exp(confint(itt_siptw_model))[2]

result_chart_stab[3, 'estimate'] <- exp(at_model$coef)
result_chart_stab[3, 'lower_ci'] <- exp(confint(at_model))[1]
result_chart_stab[3, 'upper_ci'] <- exp(confint(at_model))[2]

result_chart_stab[4, 'estimate'] <- exp(at_siptw_model$coef)
result_chart_stab[4, 'lower_ci'] <- exp(confint(at_siptw_model))[1]
result_chart_stab[4, 'upper_ci'] <- exp(confint(at_siptw_model))[2]

result_chart_stab[5, 'estimate'] <- exp(at_sipcw_model$coef)
result_chart_stab[5, 'lower_ci'] <- exp(confint(at_sipcw_model))[1]
result_chart_stab[5, 'upper_ci'] <- exp(confint(at_sipcw_model))[2]

result_chart_stab[6, 'estimate'] <- exp(at_siptw_sipcw_model$coef)
result_chart_stab[6, 'lower_ci'] <- exp(confint(at_siptw_sipcw_model))[1]
result_chart_stab[6, 'upper_ci'] <- exp(confint(at_siptw_sipcw_model))[2]

result_chart_stab[7, 'estimate'] <- exp(at_sipcw_nl_model$coef)
result_chart_stab[7, 'lower_ci'] <- exp(confint(at_sipcw_nl_model))[1]
result_chart_stab[7, 'upper_ci'] <- exp(confint(at_sipcw_nl_model))[2]

result_chart_stab[8, 'estimate'] <- exp(at_sipcw_nl_siptw_model$coef)
result_chart_stab[8, 'lower_ci'] <- exp(confint(at_sipcw_nl_siptw_model))[1]
result_chart_stab[8, 'upper_ci'] <- exp(confint(at_sipcw_nl_siptw_model))[2]

result_chart_stab[9, 'estimate'] <- exp(at_sipcw_pooled_model$coef)
result_chart_stab[9, 'lower_ci'] <- exp(confint(at_sipcw_pooled_model))[1]
result_chart_stab[9, 'upper_ci'] <- exp(confint(at_sipcw_pooled_model))[2]

result_chart_stab[10, 'estimate'] <- exp(at_sipcw_pooled_siptw_model$coef)
result_chart_stab[10, 'lower_ci'] <- exp(confint(at_sipcw_pooled_siptw_model))[1]
result_chart_stab[10, 'upper_ci'] <- exp(confint(at_sipcw_pooled_siptw_model))[2]

result_chart_stab[11, 'estimate'] <- exp(at_sipcw_nl_mod_model$coef)
result_chart_stab[11, 'lower_ci'] <- exp(confint(at_sipcw_nl_mod_model))[1]
result_chart_stab[11, 'upper_ci'] <- exp(confint(at_sipcw_nl_mod_model))[2]

result_chart_stab[12, 'estimate'] <- exp(at_sipcw_nl_mod_siptw_model$coef)
result_chart_stab[12, 'lower_ci'] <- exp(confint(at_sipcw_nl_mod_siptw_model))[1]
result_chart_stab[12, 'upper_ci'] <- exp(confint(at_sipcw_nl_mod_siptw_model))[2]

result_chart_stab <- cbind(model = row.names(result_chart_stab), result_chart_stab, row.names = NULL)
write_xlsx(result_chart_stab, paste(path_results, 'result_chart_stab.xlsx', sep ='/'))

#summary_cox_results <- result_chart_stab[c(2,4,6,8,10),]
summary_cox_results <- result_chart_stab[c(2,4,6,8,10,12),]
summary_cox_results$model <- factor(summary_cox_results$model, levels = rev(summary_cox_results$model))

setwd(path_results)

# non-stabilized weights

itt_model <- readRDS(paste(path_results, 'cox_no_weight.rds', sep = '/'))
itt_siptw_model <- readRDS(paste(path_results, 'cox_siptw.rds', sep = '/'))
at_model <- readRDS(paste(path_results, 'cox_at.rds', sep = '/'))
at_siptw_model <- readRDS(paste(path_results, 'cox_at_siptw.rds', sep = '/'))
at_ipcw_model <- readRDS(paste(path_results, 'cox_lagged.rds', sep = '/'))
at_iptw_ipcw_model <- readRDS(paste(path_results, 'cox_lagged_iptw.rds', sep = '/'))
at_ipcw_nl_model <- readRDS(paste(path_results, 'cox_nl.rds', sep = '/'))
at_ipcw_nl_iptw_model <- readRDS(paste(path_results, 'cox_nl_iptw.rds', sep = '/'))
at_ipcw_pooled_model <- readRDS(paste(path_results, 'cox_pooled.rds', sep = '/'))
at_ipcw_pooled_iptw_model <- readRDS(paste(path_results, 'cox_pooled_iptw.rds', sep = '/'))
at_ipcw_nl_mod_model <- readRDS(paste(path_results, 'cox_nl_mod.rds', sep = '/'))
at_ipcw_nl_mod_iptw_model <- readRDS(paste(path_results, 'cox_nl_mod_iptw.rds', sep = '/'))

result_chart <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(result_chart) <- c('estimate', 'lower_ci', 'upper_ci')

rownames(result_chart) <- c('ITT', 'ITT (IPTW)', 'AT', 'AT (IPTW)', 'AT (lagged IPCW)', 'AT (lagged IPCW + IPTW)',
                                 'AT (non-lagged IPCW)', 'AT (non-lagged IPCW + IPTW)', 'AT (pooled IPCW)',
                                 'AT (pooled IPCW + IPTW)', 'AT (mod non-lagged IPCW)', 'AT (mod non-lagged IPCW + IPTW)')

result_chart[1, 'estimate'] <- exp(itt_model$coef)
result_chart[1, 'lower_ci'] <- exp(confint(itt_model))[1]
result_chart[1, 'upper_ci'] <- exp(confint(itt_model))[2]

result_chart[2, 'estimate'] <- exp(itt_siptw_model$coef)
result_chart[2, 'lower_ci'] <- exp(confint(itt_siptw_model))[1]
result_chart[2, 'upper_ci'] <- exp(confint(itt_siptw_model))[2]

result_chart[3, 'estimate'] <- exp(at_model$coef)
result_chart[3, 'lower_ci'] <- exp(confint(at_model))[1]
result_chart[3, 'upper_ci'] <- exp(confint(at_model))[2]

result_chart[4, 'estimate'] <- exp(at_siptw_model$coef)
result_chart[4, 'lower_ci'] <- exp(confint(at_siptw_model))[1]
result_chart[4, 'upper_ci'] <- exp(confint(at_siptw_model))[2]

result_chart[5, 'estimate'] <- exp(at_ipcw_model$coef)
result_chart[5, 'lower_ci'] <- exp(confint(at_ipcw_model))[1]
result_chart[5, 'upper_ci'] <- exp(confint(at_ipcw_model))[2]

result_chart[6, 'estimate'] <- exp(at_iptw_ipcw_model$coef)
result_chart[6, 'lower_ci'] <- exp(confint(at_iptw_ipcw_model))[1]
result_chart[6, 'upper_ci'] <- exp(confint(at_iptw_ipcw_model))[2]

result_chart[7, 'estimate'] <- exp(at_ipcw_nl_model$coef)
result_chart[7, 'lower_ci'] <- exp(confint(at_ipcw_nl_model))[1]
result_chart[7, 'upper_ci'] <- exp(confint(at_ipcw_nl_model))[2]

result_chart[8, 'estimate'] <- exp(at_ipcw_nl_iptw_model$coef)
result_chart[8, 'lower_ci'] <- exp(confint(at_ipcw_nl_iptw_model))[1]
result_chart[8, 'upper_ci'] <- exp(confint(at_ipcw_nl_iptw_model))[2]

result_chart[9, 'estimate'] <- exp(at_ipcw_pooled_model$coef)
result_chart[9, 'lower_ci'] <- exp(confint(at_ipcw_pooled_model))[1]
result_chart[9, 'upper_ci'] <- exp(confint(at_ipcw_pooled_model))[2]

result_chart[10, 'estimate'] <- exp(at_ipcw_pooled_iptw_model$coef)
result_chart[10, 'lower_ci'] <- exp(confint(at_ipcw_pooled_iptw_model))[1]
result_chart[10, 'upper_ci'] <- exp(confint(at_ipcw_pooled_iptw_model))[2]

result_chart[11, 'estimate'] <- exp(at_ipcw_nl_mod_model$coef)
result_chart[11, 'lower_ci'] <- exp(confint(at_ipcw_nl_mod_model))[1]
result_chart[11, 'upper_ci'] <- exp(confint(at_ipcw_nl_mod_model))[2]

result_chart[12, 'estimate'] <- exp(at_ipcw_nl_mod_iptw_model$coef)
result_chart[12, 'lower_ci'] <- exp(confint(at_ipcw_nl_mod_iptw_model))[1]
result_chart[12, 'upper_ci'] <- exp(confint(at_ipcw_nl_mod_iptw_model))[2]

result_chart <- cbind(model = row.names(result_chart), result_chart, row.names = NULL)
write_xlsx(result_chart, paste(path_results, 'result_chart.xlsx', sep ='/'))

#summary_cox_results <- result_chart[c(2,4,6,8,10),]
summary_cox_results <- result_chart[c(2,4,6,8,10,12),]
summary_cox_results$model <- factor(summary_cox_results$model, levels = rev(summary_cox_results$model))

setwd(path_results)


#### IR RESULTS ####

bootstrap_ci <- readRDS(paste(path_results, 'bootstrap_ci.rds', sep = '/'))
time <- data.frame(time = c(unique(cohort$month_year))) %>% 
  arrange(time)

# ITT IPTW: itt.2019-01.iptw.IR
itt_iptw_IR <- bootstrap_ci %>%
  filter(grepl("^itt.*iptw\\.IR$", variable)) %>% 
  select(-variable)

itt_iptw_IR <- itt_iptw_IR[-1,]
itt_iptw_IR <- cbind(time, variable = 'ITT (IPTW)', itt_iptw_IR)

# AT IPTW: att.2019-01.iptw.IR
at_iptw_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_IR <- at_iptw_IR[-(1:4),]
at_iptw_IR <- cbind(time, variable = 'AT (IPTW)', at_iptw_IR)

# AT IPTW + IPCW: att.2019-01.iptw.ipcw.IR
at_iptw_ipcw_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_IR <- cbind(time, variable = 'AT (lagged IPCW + IPTW)', at_iptw_ipcw_IR)

# AT IPTW + nl IPCW: att.2019-01.iptw.ipcw_nl.IR
at_iptw_ipcw_nl_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_nl\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_nl_IR <- cbind(time, variable = 'AT (non-lagged IPCW + IPTW)', at_iptw_ipcw_nl_IR)

# AT IPTW + pooled IPCW: att.2019-01.iptw.ipcw_pooled.IR
at_iptw_ipcw_pooled_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_pooled\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_pooled_IR <- cbind(time, variable = 'AT (pooled IPCW + IPTW)', at_iptw_ipcw_pooled_IR)

# AT IPTW + mod nl IPCW: att.2019-01.iptw.ipcw_nl_mod.IR
at_iptw_ipcw_nl_mod_IR <- bootstrap_ci %>%
  filter(grepl("^at.*iptw\\.ipcw_nl_mod\\.IR$", variable)) %>% 
  select(-variable)

at_iptw_ipcw_nl_mod_IR <- cbind(time, variable = 'AT (mod non-lagged IPCW + IPTW)', at_iptw_ipcw_nl_mod_IR)

# Summary plot

IR_results <- rbind(itt_iptw_IR, at_iptw_IR, at_iptw_ipcw_IR, at_iptw_ipcw_nl_IR, at_iptw_ipcw_pooled_IR, at_iptw_ipcw_nl_mod_IR)

png(filename="IR_over_time_by_model.png", width=6000, height=4000, res=500)

ggplot(IR_results, aes(x = time, y = estimate*100, color = variable, group = variable)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "calendar time", y = "incidence rate (per 100 person-years)", title = "Incidence Rate over Time by Model") +
  scale_color_viridis_d(option = "H", begin = 0.1, end = 0.9, direction = -1, name = 'Model') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  ggtitle("Incidence Rate over Time by Model")

dev.off()

# Plot with CI for one

ggplot(itt_iptw_IR, aes(x = time, y = estimate*100, color = variable, group = variable)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "calendar time", y = "incidence rate (per 100 person-years)", title = "Incidence Rate over Time by Model") +
  scale_color_viridis_d(option = "H", begin = 0.1, end = 0.9, direction = -1, name = 'Model') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  ggtitle("Incidence Rate over Time by Model")


# Create the plot with confidence intervals
save <- itt_iptw_IR
itt_iptw_IR <- save

itt_iptw_IR$time <- as.Date(paste0(itt_iptw_IR$time, "-01"), format = "%Y-%m-%d")

ggplot(itt_iptw_IR, aes(x = time, y = estimate * 100, color = variable, group = variable)) +
  geom_point(size = 2) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci * 100, ymax = upper_ci * 100, fill = variable), alpha = 0.2, color = NA) +
  labs(x = "Calendar Time", y = "Incidence Rate (per 100 person-years)", title = "Incidence Rate over Time by Model") +
  scale_color_viridis_d(option = "H", begin = 0.1, end = 0.9, direction = -1, name = 'Model') +
  scale_fill_viridis_d(option = "H", begin = 0.1, end = 0.9, direction = -1, name = 'Model', guide = "none") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") +  # Adjust date labels and breaks
  scale_y_continuous(labels = scales::label_comma()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  ggtitle("Incidence Rate over Time by Model")


