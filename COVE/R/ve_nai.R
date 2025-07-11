# ----------------------------------------------------------------------------
# Main function for VE NAI analyses
# ----------------------------------------------------------------------------

ve_nai <- function(data,
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
                   n_boot = n_boot,
                   return_models = FALSE,
                   id_name = NULL, 
                   seed = 12345){
  
  # set seed for reproducibility
  set.seed(seed)
  
  # Create formulas for symptomatic and asymptomatic Cox models, LR model
  symp_formula <- as.formula(paste0("survival::Surv(", time_name, ", ", event_name, " == ", symp_level, ") ~ ",
                                    vax_name, " + ", paste(covariate_names, collapse = "+")))
  
  asymp_formula <- as.formula(paste0("survival::Surv(", time_name, ", ", event_name, " == ", asymp_level, ") ~ ",
                                    vax_name, " + ", paste(covariate_names, collapse = "+")))
  
  symp_lr_formula <-  as.formula(paste0("symp_ind ~ ", vax_name, " + ", paste(covariate_names, collapse = "+")))
  
  # Fit cox models for symptomatic & asymptomatic infections, symptomatic GLM
  symp_cox_fit <- survival::coxph(symp_formula, data = data)
  
  asymp_cox_fit <- survival::coxph(asymp_formula, data = data)

  symp_lr_fit <- glm(symp_lr_formula, 
    # subset to infected participants only (given uninfected == 0)
    data = data[data[[event_name]] > 0, ],
    family = stats::binomial()
  )
  
  # Get point estimates
  ve_fit <- ve_from_fits(
    symp_cox_fit = symp_cox_fit, 
    asymp_cox_fit = asymp_cox_fit, 
    symp_lr_fit = symp_lr_fit, 
    data = data, 
    t0 = t0
  )
  
  # Bootstrap estimates
  boot_est <- bootstrap_estimates(data = data,
                                  asymp_formula = asymp_formula,
                                  symp_formula = symp_formula,
                                  symp_lr_formula = symp_lr_formula,
                                  n_boot = n_boot,
                                  t0 = t0,
                                  id_name = id_name,
                                  event_name = event_name)
  
  if(return_models){
    out <- list(ve_fit = ve_fit,
                boot_est = boot_est,
                symp_cox_fit = symp_cox_fit,
                asymp_cox_fit = asymp_cox_fit,
                symp_lr_fit = symp_lr_fit)
  } else{
    out <- list(ve_fit = ve_fit,
                boot_est = boot_est,
                symp_cox_fit = NULL,
                asymp_cox_fit = NULL,
                symp_lr_fit =  NULL)
  }
  
  class(out) <- "ve_nai"
  
  return(out)

}
