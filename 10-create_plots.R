## ---------------------------
##
## Program: 10. Create plots
##
## Purpose: Create plots and tables for thesis/publication. 
## 
##
## Author: Gwen Aubrac
##
## Date Created: 2024-08-06
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------


#### LOAD PACKAGES ####

# for plots
library(lubridate)
library(dplyr)
library(magrittr)
library(fastDummies)
library(ggplot2)
library(tidyr)
library(cobalt)
library(scales)
library(survival)
library(cowplot)
library(table1)
library(writexl)

# for analyses

colors_trt = c('#0d0887', '#fdbd22')
colors_adj = c('#f56864', '#6a00a8')

'#6a00a8'
'#9a00b1'
'#d53e4f'
'#f56864'
'#f9a463'
'#fdbd22'
'#f9c43d'
'#f8e4a7'
'#f7f4f9'

#### DEFINE PATHS ####

path_cohort <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main"
path_main <- "Z:/EPI/Protocol 24_004042/Gwen/data/cohort/main"
path_results <- "Z:/EPI/Protocol 24_004042/Gwen/results/main" 

covariates <- readRDS(file = paste(path_main, 'covariates.rds', sep = '/'))
comorbidities <- readRDS(file = paste(path_main, 'comorbidities.rds', sep = '/'))
base_comorb <- readRDS(file = paste(path_main, 'base_comorb.rds', sep = '/'))
dec_comorb <- readRDS(file = paste(path_main, 'dec_comorb.rds', sep = '/'))

path_plots <- "Z:/EPI/Protocol 24_004042/Gwen/results/plots"
setwd(path_plots) 

#### TABLE 1 ####

# data formatting
cohort <- readRDS(paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))

cohort %<>%
  mutate(disc_one_year = as.factor(if_else(disc == 1 & disc_date < ymd(entry_date+365), 1, 0)))

tab_cohort <- cohort
tab_cohort %<>% mutate (censor = as.factor(censor))

label(tab_cohort$disc_one_year) <- 'Discontinued treatment within one year'
label(tab_cohort$sex) <- 'Sex'
label(tab_cohort$age_at_entry) <- 'Age'
label(tab_cohort$age_group) <- 'Age group'
label(tab_cohort$itt_follow_up) <- 'ITT follow-up'
label(tab_cohort$at_follow_up) <- 'AT follow-up'
label(tab_cohort$disc) <- 'Discontinued treatment'
label(tab_cohort$switch) <- 'Switched treatment'
label(tab_cohort$depression_base) <- 'Depression'
label(tab_cohort$hepatic_disease_base) <- 'Hepatic disease'
label(tab_cohort$hyperlipidemia_base) <- 'Hyperlipidemia'
label(tab_cohort$hypertension_base) <- 'Hypertension'
label(tab_cohort$chronic_renal_disease_base) <- 'Chronic renal disease'
label(tab_cohort$heart_failure_base) <- 'Heart failure'
label(tab_cohort$ischemic_heart_disease_base) <- 'Ischemic heart disease'
label(tab_cohort$lvh_base) <- 'Left ventricular hypertrophy'
label(tab_cohort$cerebrovascular_disease_base) <- 'Cerebrovascular disease'
label(tab_cohort$arrhythmia_base) <- 'Arrhythmia'
label(tab_cohort$pvd_base) <- 'Peripheral vascular disease'
label(tab_cohort$epilepsy_base) <- 'Epilepsy'
label(tab_cohort$valvular_heart_disease_base) <- 'Valvular heart disease'
label(tab_cohort$pacemaker_base) <- 'Pacemaker'
label(tab_cohort$cardiomyopathy_base) <- 'Cardiomyopathy'
label(tab_cohort$hypokalemia_base) <- 'Hypokalemia'
label(tab_cohort$hypomagnesemia_base) <- 'Hypomagnesemia'
label(tab_cohort$acute_renal_disease_base) <- 'Acute renal disease'
label(tab_cohort$stroke_base) <- 'Stroke'
label(tab_cohort$suicidal_ideation_self_harm_base) <- 'Suicidal ideation and self-harm'
label(tab_cohort$myocardial_infarction_base) <- 'Myocardial infarction'
label(tab_cohort$copd_base) <- 'COPD'
label(tab_cohort$hypocalcemia_base) <- 'Hypocalcemia'
label(tab_cohort$deprivation) <- 'Deprivation Index'
label(tab_cohort$ethnicity) <- 'Ethnicity'

units(tab_cohort$age_at_entry) <- 'years'
units(tab_cohort$itt_follow_up) <- 'days'
units(tab_cohort$at_follow_up) <- 'days'

caption <- 'Baseline Characteristics of Cohort'

# order by most common baseline comorb
order_comorb <- cohort
order_comorb[base_comorb] <- lapply(order_comorb[base_comorb], function(x) as.numeric(as.character(x)))

comorb_freq <- order_comorb %>%
  summarise(across(all_of(base_comorb), ~ mean(.x, na.rm = TRUE) * 100))

comorb_freq_long <- comorb_freq %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Frequency") %>%
  arrange(desc(Frequency))

ordered_base_comorb <- comorb_freq_long$Variable

# create table
table1_formula <- as.formula(paste("~", paste(c("sex", "age_at_entry", "age_group", "deprivation", "ethnicity", "disc_one_year", "itt_follow_up", "at_follow_up", "disc", "switch", ordered_base_comorb), collapse = "+"), "|trt"))
table1 <- table1(
  table1_formula,
  data = tab_cohort,
  overall = c(right = 'Total'),
  caption = caption
)

table1
write.table (table1 , "table1.csv", col.names = T, row.names=F, append= F, sep=',')

#### COVARIATE ASSOCIATION OVER TIME ####

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))

base_variables <- c('hypertension_base', 'depression_base', 'hyperlipidemia_base', 'suicidal_ideation_self_harm_base') # most common comorb
base_variables <- c(covariates) # specified covariates
grouping_var <- 'month_year' 
trt_col_name <- 'trt_dummy'

cov_desc <- data.frame(group = c(unique(cohort[,grouping_var]))) %>% 
  arrange(group)

for (i in 1:length(base_variables)) {
  var <- base_variables[i]
  
  cols <- ((colnames(cohort) == grouping_var) |
             (colnames(cohort) == trt_col_name) |
             (colnames(cohort) == var))
  
  # extract variable info from cohort
  var_df <- cohort[cols]
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
    
    cov_desc$var <- NA
    
    # calculate association with exposure through logistic regression
    for (i in 1:nrow(cov_desc)) {
      group_i <- cov_desc[i,1]
      cov_desc[cov_desc$'group' == group_i, 'var'] <- glm(trt ~ var, data = subset(var_sum, group == group_i), family = 'binomial')$coefficients[[2]]
    }
    
    names(cov_desc)[names(cov_desc) == 'var'] <- name
  }
}

cov_desc_long <- cov_desc %>% 
  pivot_longer(cols = -group, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

cov_desc_long <- cov_desc_long %>%
  mutate(
    variable = variable %>%
      sub("_base$", "", .) %>%
      gsub("_", " ", .)
  )

cov_desc_long$group <- as.character(cov_desc_long$group)
cov_desc_long$group <- as.Date(paste0(cov_desc_long$group, "-01"), format = "%Y-%m-%d")

p <- ggplot(cov_desc_long, aes(x = group, y = value, color = variable, group = variable)) +
  geom_point(size = 1) +
  geom_line(size = 0.5) +
  labs(x = "calendar time", y = "coefficient of association") +
  scale_color_viridis_d(option = "C", begin = 0.1, end = 0.9, direction = -1, name = 'Covariate') +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'), 
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "bottom", 
    legend.title = element_text(size = 8, face = 'bold'),
    legend.text = element_text(size = 8),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  scale_x_date(
    breaks = date_breaks("6 months"), 
    labels = date_format("%b %Y")
  ) + 
  ggtitle("Association Between Covariates and Treatment over Time")

p

ggsave("covariate_association_plot.png", plot = p, width = 6, height = 4, units = "in", bg = 'white')

dev.off()

rm(cov_desc, cov_desc_long, p, var_df, var_sum)
rm(cols, cols_keep, group_i, grouping_var, i, j, name, trt_col_name, var, var_i)
rm(cohort)

#### PROPENSITY SCORE DENSITY ####

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_covariates.rds', sep = '/'))

cohort %<>% mutate(trt_dummy = as.factor(trt_dummy), trt = as.factor(trt))

cohort <- cohort %>%
  filter (!is.na(sex))

model <- readRDS(paste(path_results, 'iptw_model.rds', sep = '/'))

weights_out <- WeightIt::weightit(
  formula = model,
  data = cohort,
  method = 'glm')

weights_out$treat <- if_else(weights_out$treat == 0, 'SNRI', 'SSRI')

p <- bal.plot(
  weights_out,
  var.name = "prop.score",
  which = "both",
  type = "density",
  mirror = FALSE,
  colors = colors_trt, 
  grid = FALSE,
  legend.labels = c('SNRI', 'SSRI'),
  position = 'bottom'
)

p <- p +
  ggtitle("Propensity Score Distributional Balance") +
  xlab("propensity score") +
  ylab("density")  

p

ggsave("propensity_score_plot.png", plot = p, width = 6, height = 4, units = "in", bg = 'white')

dev.off()

#### LOVE PLOT ####

cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))
covs <- readRDS(file = paste(path_cohort, 'iptw_covs.rds', sep = '/'))

bal_tab <- bal.tab(trt_dummy ~ covs, data = cohort,
                   weights = 'iptw',
                   binary = "std", continuous = "std", 
                   stats = c('mean.diffs'),
                   s.d.denom = 'pooled',
                   un = TRUE, 
                   thresholds = c(m = 0.1)
) 


iptw_balance_tab <- bal_tab$Balance %>% 
  dplyr::rename('SMD (Unadjusted)' = Diff.Un,
                'SMD (IPTW)' = Diff.Adj,
                'Balance' = M.Threshold)

iptw_balance_tab$Variable <- row.names(iptw_balance_tab)

write.table (iptw_balance_tab , "iptw_balance_tab.csv", col.names = T, row.names=T, append= F, sep=',')

var_names <- c(iptw_balance_tab$Variable)
loveplot_var_names <- sub("_base$", "", var_names)
loveplot_var_names <- gsub("_", " ", loveplot_var_names)
loveplot_var_names[loveplot_var_names=='sex Male'] <- 'male'
loveplot_var_names[loveplot_var_names=='pvd'] <- 'perivascular disease'
loveplot_var_names[loveplot_var_names=='lvh'] <- 'left ventricular heart disease'
loveplot_var_names[loveplot_var_names=='copd'] <- 'chronic obstructive pulmonary disease'
loveplot_var_names <- setNames(loveplot_var_names, var_names)
loveplot_var_names

p <- love.plot(
  bal_tab,
  binary = "std",
  thresholds = c(m = .1),
  title = 'Covariate Balance with IPTW',
  colors = colors_adj,
  var.names = loveplot_var_names,
  position = 'bottom',
  var.order = 'unadjusted',
  size = 2
)


p <- p + theme(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title.x = element_text(size = 8, face = 'bold'),
  plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
  
  legend.position = "bottom", 
  legend.title = element_text(size = 8, face = 'bold'),
  legend.text = element_text(size = 8),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(0, 0, 0, 0),
  legend.key.size = unit(0.5, "lines")
) 

p

ggsave("loveplot.png", plot = p, width = 6, height = 6, units = "in", bg = 'white')

dev.off()

rm(covs, p, iptw_balance_tab, bal_tab)
rm(cohort)

#### CUMULATIVE INCIDENCE BY TRT ####

cohort_long <- readRDS(file = paste(path_results, 'cohort_analytic_at.rds', sep = '/'))
cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))

## REFERENCE GROUP ##

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

cohort_ref <- cohort %>% 
  filter(trt_dummy == 0)

ref_itt <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref)
ref_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date - entry_date), itt_event) ~ 1, data = cohort_ref, weights = siptw)
ref_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref)
ref_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = iptw)

ref_at_iptw_ipcw_lagged <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw*iptw)
ref_at_iptw_ipcw_nl <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = sipcw_nl*siptw)
ref_at_iptw_ipcw_pooled <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = sipcw_pooled*siptw)
ref_at_iptw_ipcw_nl_mod <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = sipcw_nl_mod*siptw)

# ITT vs AT

y_limit <- c(0, 0.03) 

png("snri_incidence_by_analysis.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SNRIs",
     col = '#f56864',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

lines(ref_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(ref_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8'),
       lty = c(2, 1, 2, 1), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

y_limit <- c(0, 0.03) 

png(filename = "snri_incidence_by_model.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_at_iptw_ipcw_lagged,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SNRIs",
     col = '#9a00b1',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw_ipcw_nl,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pooled,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_nl_mod,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + lagged IPCW)", "AT (IPTW + non-lagged IPCW)", 
                  "AT (IPTW + pooled IPCW)", "AT (IPTW + modified nl IPCW"),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'),
       lty = c(2, 1, 1, 1, 1), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


## COMPARATOR GROUP ## 

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1)

cohort_comp <- cohort %>% 
  filter(trt_dummy == 1)

comp_itt <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp)
comp_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp, weights = iptw)
comp_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp)
comp_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw)

comp_at_iptw_ipcw_lagged <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw*iptw)
comp_at_iptw_ipcw_nl <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_nl*iptw)
comp_at_iptw_ipcw_pooled <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_pooled*iptw)
comp_at_iptw_ipcw_nl_mod <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = ipcw_nl_mod*iptw)

# ITT vs AT

y_limit <- c(0, 0.03) 

png("ssri_incidence_by_analysis.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SSRIs",
     col = '#f56864',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

lines(comp_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(comp_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8'),
       lty = c(2, 1, 2, 1), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

y_limit <- c(0, 0.03)  

png(filename = "ssri_incidence_by_model.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_at_iptw_ipcw_lagged,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SSRIs",
     col = '#9a00b1',
     ylim = y_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw_ipcw_nl,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pooled,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_nl_mod,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + lagged IPCW)", "AT (IPTW + non-lagged IPCW)", 
                  "AT (IPTW + pooled IPCW)", "AT (IPTW + modified nl IPCW"),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'),
       lty = c(2, 1, 1, 1, 1), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()

#### CUMULATIVE INCIDENCE ZOOM IN ####

cohort_long <- readRDS(file = paste(path_results, 'cohort_analytic_at.rds', sep = '/'))
cohort <- readRDS(file = paste(path_cohort, 'antidepressant_cohort_iptw.rds', sep = '/'))
times_dec <- readRDS(paste(path_main, 'times_dec.rds', sep = '/'))

## REFERENCE GROUP ## 

cohort_long_ref <- cohort_long %>% 
  filter(trt_dummy == 0)

cohort_ref <- cohort %>% 
  filter(trt_dummy == 0)

ref_itt <- survfit(Surv(itt_follow_up, itt_event) ~ 1, data = cohort_ref)
ref_itt_iptw <- survfit(Surv(itt_follow_up, itt_event) ~ 1, data = cohort_ref, weights = iptw)
ref_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref)
ref_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = iptw)

ref_at_iptw_ipcw_lagged <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw*iptw)
ref_at_iptw_ipcw_nl <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_nl*iptw)
ref_at_iptw_ipcw_pooled <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_pooled*iptw)
ref_at_iptw_ipcw_nl_mod <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_ref, weights = ipcw_nl_mod*iptw)

# ITT vs AT

y_limit <- c(0, 0.02) 
x_limit <- c(0, 500)

png("snri_incidence_by_analysis_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SNRIs",
     col = '#f56864',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(ref_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(ref_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)", "intervals"),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8', '#f8e4a7'),
       lty = c(2, 1, 2, 1, 3), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

png(filename = "snri_incidence_by_model_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(ref_at_iptw_ipcw_lagged,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SNRIs",
     col = '#9a00b1',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(ref_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(ref_at_iptw_ipcw_nl,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_pooled,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(ref_at_iptw_ipcw_nl_mod,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + lagged IPCW)", "AT (IPTW + non-lagged IPCW)", 
                  "AT (IPTW + pooled IPCW)", "AT (IPTW + modified nl IPCW", 'intervals'),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d', '#f8e4a7'),
       lty = c(2, 1, 1, 1, 1, 3), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


## COMPARATOR GROUP ## 

cohort_long_comp <- cohort_long %>% 
  filter(trt_dummy == 1) 

cohort_comp <- cohort %>% 
  filter(trt_dummy == 1)

comp_itt <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp)
comp_itt_iptw <- survfit(Surv(as.numeric(itt_exit_date-entry_date), itt_event) ~ 1, data = cohort_comp, weights = iptw)
comp_at <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp)
comp_at_iptw <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw)

comp_at_iptw_ipcw_lagged <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw*ipcw)
comp_at_iptw_ipcw_nl <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw*ipcw_nl)
comp_at_iptw_ipcw_pooled <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw*ipcw_pooled)
comp_at_iptw_ipcw_nl_mod <- survfit(Surv(Tstart, Tstop, event_at_tstop) ~ 1, data = cohort_long_comp, weights = iptw*ipcw_nl_mod)

# ITT vs AT

y_limit <- c(0, 0.02) 
x_limit <- c(0, 500)

png("ssri_incidence_by_analysis_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_itt,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SSRIs",
     col = '#f56864',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 2,
     lwd = 2
)

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(comp_itt_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f56864',
      lty = 1,
      lwd = 2)

lines(comp_at,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 1,
      lwd = 2)

legend("topleft", 
       legend = c("ITT", "ITT (IPTW)", "AT", "AT (IPTW)", 'intervals'),
       col = c('#f56864', '#f56864', '#6a00a8', '#6a00a8', '#f8e4a7'),
       lty = c(2, 1, 2, 1, 3), 
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()


# comparing censoring models

y_limit <- c(0, 0.02) 
x_limit <- c(0, 500)

png(filename = "ssri_incidence_by_model_zoom.png", width = 1800, height = 1800, res = 300)

par(lab = c(10, 10, 7), 
    font.lab = 2, 
    bty = 'n'
)

plot(comp_at_iptw_ipcw_lagged,
     fun = function(x) 1 - x,
     ylab = "cumulative probability of death",
     xlab = "days since cohort entry",
     main = "Cumulative Incidence Rate for SSRIs",
     col = '#9a00b1',
     ylim = y_limit,
     xlim = x_limit,
     conf.int = FALSE, 
     lty = 1,
     lwd = 2
) 

abline(v = times_dec, col = '#f8e4a7', lty = 3, lwd = 1)

lines(comp_at_iptw,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#6a00a8',
      lty = 2,
      lwd = 2)

lines(comp_at_iptw_ipcw_nl,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#d53e4f', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_pooled,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9a463', 
      lwd = 2
)

lines(comp_at_iptw_ipcw_nl_mod,
      fun = function(x) 1 - x,
      conf.int = FALSE, 
      col = '#f9c43d', 
      lwd = 2
)

legend("topleft", 
       legend = c("AT (IPTW)", "AT (IPTW + lagged IPCW)", "AT (IPTW + non-lagged IPCW)", 
                  "AT (IPTW + pooled IPCW)", "AT (IPTW + modified nl IPCW", 'quartiles'),
       col = c( '#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d', '#f8e4a7'),
       lty = c(2, 1, 1, 1, 1, 3), 
       cex = 0.8,
       bty = 'n', 
       lwd = 2, 
       ncol = 1) 

dev.off()

#### FOREST PLOT - STABILIZED ####

itt_model <- readRDS(paste(path_results, 'cox_no_weight.rds', sep = '/'))
itt_iptw_model <- readRDS(paste(path_results, 'cox_siptw.rds', sep = '/'))
at_model <- readRDS(paste(path_results, 'cox_at.rds', sep = '/'))
at_iptw_model <- readRDS(paste(path_results, 'cox_at_siptw.rds', sep = '/'))
at_ipcw_model <- readRDS(paste(path_results, 'cox_stab_lagged.rds', sep = '/'))
at_iptw_ipcw_model <- readRDS(paste(path_results, 'cox_stab_lagged_iptw.rds', sep = '/'))
at_ipcw_nl_model <- readRDS(paste(path_results, 'cox_nl_stab.rds', sep = '/'))
at_ipcw_nl_iptw_model <- readRDS(paste(path_results, 'cox_nl_stab_iptw.rds', sep = '/'))
at_ipcw_pooled_model <- readRDS(paste(path_results, 'cox_pooled_stab.rds', sep = '/'))
at_ipcw_pooled_iptw_model <- readRDS(paste(path_results, 'cox_pooled_stab_iptw.rds', sep = '/'))
at_ipcw_nl_mod_model <- readRDS(paste(path_results, 'cox_nl_mod_stab.rds', sep = '/'))
at_ipcw_nl_mod_iptw_model <- readRDS(paste(path_results, 'cox_nl_mod_stab_iptw.rds', sep = '/'))

result_chart <- data.frame(matrix(nrow = 12, ncol = 3))
colnames(result_chart) <- c('estimate', 'lower_ci', 'upper_ci')
rownames(result_chart) <- c('ITT', 'ITT (IPTW)', 'AT', 'AT (IPTW)', 'AT (lagged IPCW)', 'AT (lagged IPCW + IPTW)',
                                 'AT (non-lagged IPCW)', 'AT (non-lagged IPCW + IPTW)', 'AT (pooled IPCW)',
                                 'AT (pooled IPCW + sIPTW)', 'AT (mod non-lagged IPCW)', 'AT (mod non-lagged IPCW + IPTW)')

result_chart[1, 'estimate'] <- exp(itt_model$coef)
result_chart[1, 'lower_ci'] <- exp(confint(itt_model))[1]
result_chart[1, 'upper_ci'] <- exp(confint(itt_model))[2]

result_chart[2, 'estimate'] <- exp(itt_iptw_model$coef)
result_chart[2, 'lower_ci'] <- exp(confint(itt_iptw_model))[1]
result_chart[2, 'upper_ci'] <- exp(confint(itt_iptw_model))[2]

result_chart[3, 'estimate'] <- exp(at_model$coef)
result_chart[3, 'lower_ci'] <- exp(confint(at_model))[1]
result_chart[3, 'upper_ci'] <- exp(confint(at_model))[2]

result_chart[4, 'estimate'] <- exp(at_iptw_model$coef)
result_chart[4, 'lower_ci'] <- exp(confint(at_iptw_model))[1]
result_chart[4, 'upper_ci'] <- exp(confint(at_iptw_model))[2]

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
write_xlsx(result_chart, paste(path_results, 'cox_ratios.xlsx', sep ='/'))

summary_cox_results <- result_chart[c(2,4,6,8,10,12),]
summary_cox_results$model <- factor(summary_cox_results$model, levels = rev(summary_cox_results$model))

model_order <- rev(c("ITT (IPTW)", "AT (IPTW)", "AT (lagged sIPCW + IPTW)", "AT (non-lagged IPCW + IPTW)", "AT (pooled IPCW + sIPTW)", "AT (mod non-lagged IPCW + IPTW)"))
summary_cox_results$model <- factor(summary_cox_results$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(summary_cox_results, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "estimate (95% CI)", y = "model", title = "Forest Plot of Hazard Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 2)),
  ' (',round(summary_cox_results$lower_ci, 2),', ',
  round(summary_cox_results$upper_ci, 2), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()



#### FOREST PLOT - UNSTABILIZED ####

itt_model <- readRDS(paste(path_results, 'cox_no_weight.rds', sep = '/'))
itt_iptw_model <- readRDS(paste(path_results, 'cox_iptw.rds', sep = '/'))
at_model <- readRDS(paste(path_results, 'cox_at.rds', sep = '/'))
at_iptw_model <- readRDS(paste(path_results, 'cox_at_iptw.rds', sep = '/'))
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
                            'AT (pooled IPCW + sIPTW)', 'AT (mod non-lagged IPCW)', 'AT (mod non-lagged IPCW + IPTW)')

result_chart[1, 'estimate'] <- exp(itt_model$coef)
result_chart[1, 'lower_ci'] <- exp(confint(itt_model))[1]
result_chart[1, 'upper_ci'] <- exp(confint(itt_model))[2]

result_chart[2, 'estimate'] <- exp(itt_iptw_model$coef)
result_chart[2, 'lower_ci'] <- exp(confint(itt_iptw_model))[1]
result_chart[2, 'upper_ci'] <- exp(confint(itt_iptw_model))[2]

result_chart[3, 'estimate'] <- exp(at_model$coef)
result_chart[3, 'lower_ci'] <- exp(confint(at_model))[1]
result_chart[3, 'upper_ci'] <- exp(confint(at_model))[2]

result_chart[4, 'estimate'] <- exp(at_iptw_model$coef)
result_chart[4, 'lower_ci'] <- exp(confint(at_iptw_model))[1]
result_chart[4, 'upper_ci'] <- exp(confint(at_iptw_model))[2]

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
write_xlsx(result_chart, paste(path_results, 'cox_ratios_unstab.xlsx', sep ='/'))

summary_cox_results <- result_chart[c(2,4,6,8,10,12),]
summary_cox_results$model <- factor(summary_cox_results$model, levels = rev(summary_cox_results$model))

model_order <- rev(c("ITT (IPTW)", "AT (IPTW)", "AT (lagged sIPCW + IPTW)", "AT (non-lagged IPCW + IPTW)", "AT (pooled IPCW + sIPTW)", "AT (mod non-lagged IPCW + IPTW)"))
summary_cox_results$model <- factor(summary_cox_results$model, levels = model_order)
model_colors <- rev(c('#f56864','#6a00a8', '#9a00b1', '#d53e4f', '#f9a463', '#f9c43d'))

p <- ggplot(summary_cox_results, aes(x = estimate, y = model, color = model)) +
  geom_point(size = 3) + 
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") + 
  theme_minimal() +  
  labs(x = "estimate (95% CI)", y = "model", title = "Forest Plot of Hazard Ratios by Model") +  
  scale_color_manual(values = model_colors) +  
  theme(
    axis.text.y = element_text(size = 8), 
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "none", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(25, 0, 25, 0) 
  )

p

y_labels <- rev(paste(
  c(round(summary_cox_results$estimate, 2)),
  ' (',round(summary_cox_results$lower_ci, 2),', ',
  round(summary_cox_results$upper_ci, 2), ')',
  sep = ''
))

y_labels

dummy_data <- data.frame(
  model = factor(y_labels, levels = y_labels),
  dummy = 1 
)

p2 <- ggplot(dummy_data, aes(x = dummy, y = model)) +
  scale_y_discrete(labels = y_labels) +
  theme_void() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(40, 0, 50, 10),
    legend.position = "none" 
  )

p2

combined_plot <- plot_grid(p, p2, ncol = 2, rel_widths = c(3, 1))
combined_plot

ggsave("forest_plot_unstab.png", plot = combined_plot, width = 6, height = 3, units = "in", bg = 'white')

dev.off()




