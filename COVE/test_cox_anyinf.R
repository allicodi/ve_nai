# Test reproducing El Sahly VE Cox proportional-hazards model

# Numbers to replicate
# Placebo - 14,164; mRNA-1273 - 14,287
# SARS-CoV-2 infection: 1339 placebo, 280 mRNA-1273
# 82.0 % efficacy (79.5 - 84.2%)

here::i_am("ve_nai/test_cox_anyinf.R")

library(survival)

el_sahly_data <- readRDS(here::here("analytic_data/P3001_BlindedPhasePartA_Seconary2EfficacyData.Rds"))

# Subset to "TTINFB1"
# Per-protocol

any_inf_data <- el_sahly_data %>%
  filter(PARAMCD == "TTINFB1") %>%
  filter(PPROTFL == "Y") 

# CHECK MATCHES REPORTED N -- yes
#table(any_inf_data$TRT01P)
# mRNA-1273   Placebo 
# 14287     14164 

# PDV + Blinded phase matches numbers in table - NOT early infections
# table(any_inf_data$EVNTDESC, any_inf_data$TRT01P)
#                                                                                  mRNA-1273 Placebo
# Early SARS-CoV-2 Infection Regardless of Symptomatology or Severity                   33      95
# No SARS-CoV-2 Infection Regardless of Symptomatology or Severity                   13974   12730
# SARS-CoV-2 Infection Regardless of Symptomatology or Severity at PDV                 203     491
# SARS-CoV-2 Infection Regardless of Symptomatology or Severity in Blinded Phase        77     848

# Unclear if early infections included because this is 'SARS-CoV-2 infection' - not 'COVID-19'. Exclude for now?
any_inf_data <- any_inf_data %>%
  filter(EVNTDESC != "Early SARS-CoV-2 Infection Regardless of Symptomatology or Severity")

# ------------------------------------------------------------------------------

# "Vaccine efficacy was defined as 1 minus the hazard ratio (mRNA-1273 vs. placebo), 
# and 95% confidence intervals were estimated using a stratified Cox proportional-hazards model
# with Efron’s method of tie handling and with the treatment group as a covariate, 
# adjusted for stratification factor."

case <- c("SARS-CoV-2 Infection Regardless of Symptomatology or Severity at PDV",
          "SARS-CoV-2 Infection Regardless of Symptomatology or Severity in Blinded Phase")

analysis_data <- any_inf_data %>%
  mutate(ftype = if_else(EVNTDESC %in% case, 1, 0),
         vax = if_else(TRT01P == "mRNA-1273", 1, 0)) %>%
  select(ftype,
         vax,
         AVAL,
         STRATAR) %>%
  rename("ftime" = AVAL,
         "strat_rand" = STRATAR)

cox_model <- survival::coxph(survival::Surv(ftime, ftype == 1 ) ~ vax + strata(strat_rand), 
                         data = analysis_data, 
                         ties = "efron")

summary(cox_model)

ve <- 1 - exp(coef(cox_model))
ci <- 1 - exp(confint(cox_model))

# VE = 0.8199723 
# 95% CI = 0.7951979, 0.8417497

# so that matches too yay
