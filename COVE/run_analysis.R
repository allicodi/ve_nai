# --------------------------------------------------------------------------
# Run VE NAI analysis on COVE data with bootstrap confidence intervals
# --------------------------------------------------------------------------

here::i_am("COVE/run_analysis.R")

# Source VE functions, bootstrap functions
# source(here::here("COVE/R/ve_nai.R"))
# source(here::here("COVE/R/ve_from_fits.R"))
# source(here::here("COVE/R/bootstrap.R"))

# Test w package
devtools::load_all(here::here("../VEnai/"))

# Load data
cove_data <- readRDS(here::here("COVE/analytic_data/final_asymp_data.Rds"))

results <- ve_nai(data = cove_data,
                 vax_name = "vax",
                 time_name = "ftime",
                 event_name = "ftype",
                 symp_level = 1, 
                 asymp_level = 2,
                 covariate_names = c("Age", "Sex", "HighRiskInd",
                                     "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
                                     "Black", "Asian", "NatAmer", 
                                     "PacIsl", "Multiracial", "RaceOther",
                                     "RaceNotreported", "RaceUnknown"),
                 t0 = c(175, 180, 185, 190, 195, 200),
                 n_boot = 10,
                 return_models = FALSE)


# -------------------------------------------------------------------------------
# Old code

# Load data
cove_data <- readRDS(here::here("COVE/analytic_data/final_asymp_data.Rds"))

# Fit cox models for symptomatic & asymptomatic infections, symptomatic GLM
symp_cox_fit <- survival::coxph(survival::Surv(ftime, ftype == 1) ~ vax + 
                                  Age + Sex + HighRiskInd + 
                                  EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
                                  Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + 
                                  RaceNotreported + RaceUnknown, data = cove_data)

asymp_cox_fit <- survival::coxph(survival::Surv(ftime, ftype == 2) ~ vax + 
                                   Age + Sex + HighRiskInd + 
                                   EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
                                   Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + 
                                   RaceNotreported + RaceUnknown, data = cove_data)
symp_lr_fit <- glm(
  symp_ind ~ vax + 
    Age + Sex + HighRiskInd + 
    EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
    Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + RaceNotreported + RaceUnknown, 
  # subset to infected ptcpt only
  data = cove_data[cove_data$ftype > 0, ],
  family = stats::binomial()
)

# Get point estimate
ve_fit <- ve_from_fits(
  symp_cox_fit = symp_cox_fit, 
  asymp_cox_fit = asymp_cox_fit, 
  symp_lr_fit = symp_lr_fit, 
  data = cove_data, 
  t0 = c(175, 180, 185, 190, 195, 200)
)

# Get bootstrap estimates
pt_est <- ve_nai(data = cove_data,
                                vax_name = "vax",
                                time_name = "ftime",
                                event_name = "ftype",
                                symp_level = 1, 
                                asymp_level = 2,
                                covariate_names = c("Age", "Sex", "HighRiskInd",
                                                    "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown",
                                                    "Black", "Asian", "NatAmer", 
                                                    "PacIsl", "Multiracial", "RaceOther",
                                                    "RaceNotreported", "RaceUnknown"))
