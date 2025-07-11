# ---------------------------------------------------------------------------
# Functions to get bootstrap confidence intervals from VE NAI code
# ---------------------------------------------------------------------------

# ----------------------------------------
# Function to do one bootstrap replicate
# ----------------------------------------
one_boot <- function(data, asymp_formula, symp_formula, symp_lr_formula, n_boot, t0, id_name, event_name){
  
  # If no ID, assume all unique enrollments and just sample by idx
  # Else ID name provided, sample by ID name (in case of multiple enrollments, this makes sure all are in the bootstrap sample)
  if(is.null(id_name)){
    boot_id <- sample(1:nrow(data), replace = TRUE)
    boot_data <- data[boot_id, ,drop = FALSE]
  } else{
    boot_id <- sample(unique(data[[id_name]]), replace = TRUE)
    boot_id_rows_list <- sapply(boot_id, function(id){
      which(data[[id_name]] == id)
    })
    boot_id_rows_vec <- Reduce(c, boot_id_rows_list)
    boot_data <- data[boot_id_rows_vec, , drop = FALSE]
  }
  
  # Fit cox models for symptomatic & asymptomatic infections, symptomatic GLM
  symp_cox_fit <- survival::coxph(symp_formula, data = boot_data, model = TRUE)
  
  asymp_cox_fit <- survival::coxph(asymp_formula, data = boot_data, model = TRUE)
  
  symp_lr_fit <- glm(symp_lr_formula, 
                     # subset to infected participants only (given uninfected == 0)
                     data = boot_data[boot_data[[event_name]] > 0, ],
                     family = stats::binomial()
  )
  
  # Get point estimates
  ve_fit <- ve_from_fits(
    symp_cox_fit = symp_cox_fit, 
    asymp_cox_fit = asymp_cox_fit, 
    symp_lr_fit = symp_lr_fit, 
    data = boot_data, 
    t0 = t0
  )
  
  # Reformat results to get SE for VE_X and Ps more easily
  results_df <- data.frame(t0 = ve_fit$t0,
                           ve_i = ve_fit$ve_i,
                           ve_s = ve_fit$ve_s,
                           ve_ai = ve_fit$ve_ai,
                           ve_nai = ve_fit$ve_nai,
                           p_immune = ve_fit$p_immune,
                           p_doomed = ve_fit$p_doomed,
                           p_alwaysinf = ve_fit$p_alwaysinf,
                           p_converted = ve_fit$p_converted,
                           p_helped = ve_fit$p_helped,
                           p_helpedplus = ve_fit$p_helpedplus)
  
  return(results_df)
}

# -------------------------------------------------------
# Function to do n_boot bootstrap replicates and get SE
# -------------------------------------------------------
bootstrap_estimates <- function(data, asymp_formula, symp_formula, symp_lr_formula, n_boot, t0, id_name, event_name){
  
  # Do n_boot bootstrap replicates
  boot_estimates <- replicate(n_boot, one_boot(data = data, 
                                               asymp_formula = asymp_formula, 
                                               symp_formula = symp_formula, 
                                               symp_lr_formula = symp_lr_formula,
                                               n_boot = n_boot,
                                               t0 = t0,
                                               id_name = id_name,
                                               event_name = event_name), simplify = FALSE) 
  
  boot_res <- data.frame(do.call(rbind, boot_estimates))
  
  # Full results list
  out <- list()
  
  # For each threshold:
  for(t in unique(boot_res$t0)){
    boot_res_t <- boot_res[boot_res$t0 == t,]
    
    t0_out <- list()
    
    # Get se, lower bound, upper bound for each column in boot_res for the given threshold
    for(col in colnames(boot_res_t)[-1]){
      se <- sd(boot_res_t[,col])
      lower <- quantile(boot_res_t[,col], p = 0.025, names = FALSE)
      upper <- quantile(boot_res_t[,col], p = 0.975, names = FALSE)
      
      t0_out[[col]] <- list(se = se, lower = lower, upper = upper)
    }
    
    # Add results for given threshold to full results list
    out[[paste0("t0_", t)]] <- t0_out
    
  }
  
  return(out)
  
}

