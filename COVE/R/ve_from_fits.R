#' Compute VE for naturally asymptomatic infections and other related parameters
#' 
#' @param symp_cox_fit A \code{coxph} fit for symptomatic infections
#' @param asymp_cox_fit A \code{coxph} fit for asymptomatic infections
#' @param symp_lr_fit A \code{glm} fit for probability of symptoms given infection
#' @param data The \code{data} used to fit the above models. Should have columns \code{ftime}, \code{ftype}, 
#' 						 \code{vax}, and other columns for covariates.
#' @param t0 The time(s) of interest. The code computes cumulative incidence at these prespecified times 
#' 					 and reports VE estimates (and other estimates) at each specified time
ve_from_fits <- function(
    symp_cox_fit,
    asymp_cox_fit,
    symp_lr_fit,
    data,
    t0 = quantile(data$ftime, c(0.2, 0.5, 0.9)),
    ...
){
  n <- dim(data)[1]
  symp_chaz <- survival::basehaz(symp_cox_fit)
  asymp_chaz <- survival::basehaz(asymp_cox_fit)
  # if the unique times are not the same, then it may indicate a problem with
  # how the models were fit
  stopifnot(all(unique(symp_chaz$time) == unique(asymp_chaz$time)))
  
  # find times closest to requested t0
  closest_indices <- sapply(t0, function(t) {
    which.min(abs(symp_chaz$time - t))
  })
  subset_symp_chaz <- symp_chaz[closest_indices,]
  subset_asymp_chaz <- asymp_chaz[closest_indices,]
  
  surv_data <- merge(subset_symp_chaz, subset_asymp_chaz, by = "time", suffixes = c("_symp", "_asymp"))
  
  surv_data$total_chaz <- surv_data$hazard_symp + surv_data$hazard_asymp
  # total baseline cumulative incidence of any infection
  surv_data$total_ci <- 1 - exp(-surv_data$total_chaz)
  
  # predict from models under vaccine/no vaccine
  data0 <- data; data0$vax <- 0
  lp_symp_vax0 <- predict(symp_cox_fit, newdata = data0, type = "lp")
  lp_asymp_vax0 <- predict(asymp_cox_fit, newdata = data0, type = "lp")
  # P(Y_S = 1 | V = 0,  X = x)
  p_symp__inf_vax0 <- predict(
    symp_lr_fit, newdata = data0, type = "response"
  )
  
  data1 <- data; data1$vax <- 1
  lp_symp_vax1 <- predict(symp_cox_fit, newdata = data1, type = "lp")
  lp_asymp_vax1 <- predict(asymp_cox_fit, newdata = data1, type = "lp")
  p_symp__inf_vax1 <- predict(
    symp_lr_fit, newdata = data1, type = "response"
  )
  
  # conditional cumulative incidence any infection under no vaccine/vaccine
  # P(Y_I = 1 | V = v, X = x)
  ci_inf_vax0 <- mapply(
    hazard_symp = surv_data$hazard_symp,
    hazard_asymp = surv_data$hazard_asymp,
    FUN = function(hazard_symp, hazard_asymp){
      1 - exp(- ( hazard_symp * exp(lp_symp_vax0) + hazard_asymp * exp(lp_asymp_vax0) ))
    }
  )
  ci_inf_vax1 <- mapply(
    hazard_symp = surv_data$hazard_symp,
    hazard_asymp = surv_data$hazard_asymp,
    FUN = function(hazard_symp, hazard_asymp){
      1 - exp(- ( hazard_symp * exp(lp_symp_vax1) + hazard_asymp * exp(lp_asymp_vax1) ))
    }
  )
  
  # P(Y_I = 1 | V = 0, X = X_i), i = 1,...,n = ci_inf_vax0
  # P(Y_I = 1 | V = 1, X = X_i), i = 1,...,n = ci_inf_vax1
  # P(Y_S = 1 | V = 0, X = X_i), i = 1,...,n = p_symp__inf_vax0
  # P(Y_S = 1 | V = 1, X = X_i), i = 1,...,n = p_symp__inf_vax1
  
  # marginal parameters any infection
  # E[P(Y_I = 1 | V = 0, X)]
  marg_ci_inf_vax0 <- colMeans(ci_inf_vax0)
  # E[P(Y_I = 1 | V = 1, X)]
  marg_ci_inf_vax1 <- colMeans(ci_inf_vax1)
  ve_any_inf <- 1 - marg_ci_inf_vax1 / marg_ci_inf_vax0
  # VE_I(X_i) = 1 - ci_inf_vax1[i,] / ci_inf_vax0[i,]
  
  # conditional cumulative incidence symptomatic infection under no vaccine/vaccine
  ci_symp_inf_vax0 <- apply(ci_inf_vax0, 2, function(ci_inf){
    ci_inf * p_symp__inf_vax0
  })
  ci_symp_inf_vax1 <- apply(ci_inf_vax1, 2, function(ci_inf){
    ci_inf * p_symp__inf_vax1
  })
  
  # marignal parameters symptomatic infection
  marg_ci_symp_vax0 <- colMeans(ci_symp_inf_vax0)
  marg_ci_symp_vax1 <- colMeans(ci_symp_inf_vax1)
  ve_symp_inf <- 1 - marg_ci_symp_vax1 / marg_ci_symp_vax0
  
  # conditional cumulative incidence asymptomatic infection under no vaccine/vaccine
  # denominator for VE_NAI (and VE_AI)
  ci_asymp_inf_vax0 <- apply(ci_inf_vax0, 2, function(ci_inf){
    ci_inf * ( 1 - p_symp__inf_vax0 )
  })
  ci_asymp_inf_vax1_for_veai <- apply(ci_inf_vax1, 2, function(ci_inf){
    ci_inf * ( 1 - p_symp__inf_vax1 )
  })
  # usual numerator for VE_NAI
  ci_asymp_inf_vax1_for_venai <- apply(ci_inf_vax1, 2, function(ci_inf){
    ci_inf * ( 1 - p_symp__inf_vax0 )
  })
  
  # marginal parameters asymptomatic infection
  marg_ci_asymp_vax0 <- colMeans(ci_asymp_inf_vax0)
  marg_ci_asymp_vax1_for_veai <- colMeans(ci_asymp_inf_vax1_for_veai)
  marg_ci_asymp_vax1_for_venai <- colMeans(ci_asymp_inf_vax1_for_venai)
  ve_ai <- 1 - marg_ci_asymp_vax1_for_veai / marg_ci_asymp_vax0
  ve_nai <- 1 - marg_ci_asymp_vax1_for_venai / marg_ci_asymp_vax0
  
  
  # principal strata estimates
  p_immune__covariates <- 1 - ci_inf_vax0
  
  p_doomed__covariates <- ci_symp_inf_vax1
  
  p_alwaysinf__covariates <- ci_asymp_inf_vax1_for_venai
  
  # p_converted__covariates <- mapply(p1 = as.list(ci_symp_inf_vax1_for_veai), p2 = as.list(ci_symp_inf_vax1_for_venai),
  # 	function(p1, p2){
  # 		p1 - ( 1 - p_symp__inf_vax0 ) * p2
  # })
  p_converted__covariates <- matrix(NA, nrow = nrow(p_alwaysinf__covariates), ncol = length(t0))
  for(j in seq_along(t0)){
    p_converted__covariates[,j]  <- ci_asymp_inf_vax1_for_veai[,j] - ( 1 - p_symp__inf_vax0 ) * ci_inf_vax1[,j]
  }
  
  # p_helped__covariates <- mapply(p1 = as.list(ci_symp_inf_vax0), p2 = as.list(ci_symp_inf_vax1_for_venai),
  # 	function(p1, p2){
  # 		p1 - ( 1 - p_symp__inf_vax0 ) * p2
  # 	})
  p_helped__covariates <- matrix(NA, nrow = nrow(p_alwaysinf__covariates), ncol = length(t0))
  for(j in seq_along(t0)){
    p_helped__covariates[,j]  <- ci_asymp_inf_vax0[,j] - ( 1 - p_symp__inf_vax0 ) * ci_inf_vax1[,j]
  }
  
  
  # p_helpedplus__covariates <- mapply(
  # 	p1 = as.list(p_immune__covariates), p2 = as.list(p_doomed__covariates),
  # 	p3 = as.list(p_alwaysinf__covariates), p4 = as.list(p_converted__covariates), 
  # 	p5 = as.list(p_helped__covariates),
  # 	function(p1, p2, p3, p4, p5){
  # 		1 - (p1 + p2 + p3 + p4 + p5)
  # 	}
  # )
  p_helpedplus__covariates <- matrix(NA, nrow = nrow(p_alwaysinf__covariates), ncol = length(t0))
  for(j in seq_along(t0)){
    p_helpedplus__covariates[,j] <- 1 - (
      p_immune__covariates[,j] + 
        p_doomed__covariates[,j] + 
        p_alwaysinf__covariates[,j] + 
        p_converted__covariates[,j] + 
        p_helped__covariates[,j]
    )
  }
  
  out <- list(
    t0 = t0, ve_i = ve_any_inf, ve_s = ve_symp_inf, 
    ve_ai = ve_ai, ve_nai = ve_nai, 
    p_immune = colMeans(p_immune__covariates),
    p_immune_X = p_immune__covariates,
    p_doomed = colMeans(p_doomed__covariates),
    p_doomed_X = p_doomed__covariates,
    p_alwaysinf = colMeans(p_alwaysinf__covariates),
    p_alwaysinf_X = p_alwaysinf__covariates,
    p_converted = colMeans(p_converted__covariates),
    p_converted_X = p_converted__covariates,
    p_helped = colMeans(p_helped__covariates),
    p_helped_X = p_helped__covariates,
    p_helpedplus = colMeans(p_helpedplus__covariates),
    p_helpedplus_X = p_helpedplus__covariates
  )
  class(out) <- "ps_ve"
  return(out)
}

print.ps_ve <- function(object, what = "ve"){
  if(what == "ve"){
    print(data.frame(
      object[c("t0", "ve_i", "ve_s", "ve_ai", "ve_nai")]
    ))
  }else if(what == "ps"){
    print(data.frame(
      object[c("t0", paste0("p_", c("immune", "doomed", "alwaysinf", "converted", "helped", "helpedplus")))]
    ))
  }
}

fit_ps_msm <- function(
    object,
    data,
    t0 = object$t0[1],
    which_ps_numerator = "converted",
    which_ps_denominator = c("immune", "doomed", "helped", "helpedplus", "converted", "alwaysinf"),
    msm_covariates = "Age",
    msm_formula = "splines::ns(Age, 3)",
    msm_family = binomial(),
    ...
){
  stopifnot(t0 %in% object$t0)
  stopifnot(all(msm_covariates %in% colnames(data)))
  stopifnot(all(which_ps_numerator %in% c("immune", "doomed", "helped", "helpedplus", "converted", "alwaysinf")))
  
  if(which_ps_denominator == "naturally_symptomatic"){
    which_ps_denominator <- c("doomed", "helpedplus", "converted")
  }else if(which_ps_denominator == "naturally_asymptomatic"){
    which_ps_denominator <- c("alwaysinf", "helped")
  }else if(which_ps_denominator == "naturally_infected"){
    which_ps_denominator <- c("doomed", "helped", "helpedplus", "converted", "alwaysinf")
  }
  stopifnot(all(which_ps_denominator %in% c("immune", "doomed", "helped", "helpedplus", "converted", "alwaysinf")))
  
  msm_outcome_names_numerator <- paste0("p_", which_ps_numerator, "_X")
  msm_outcome_names_denominator <- paste0("p_", which_ps_denominator, "_X")
  
  t0_idx <- which(object$t0 == t0)
  
  msm_outcome_numerator <- Reduce("+", lapply(object[msm_outcome_names_numerator], function(x){
    x[,t0_idx]
  }))
  msm_outcome_denominator <- Reduce("+", lapply(object[msm_outcome_names_denominator], function(x){
    x[,t0_idx]
  }))
  msm_outcome <- msm_outcome_numerator / msm_outcome_denominator
  
  full_msm_formula <- paste0("msm_outcome ~ ", msm_formula)
  msm_data <- data.frame(data[,msm_covariates,drop = FALSE], msm_outcome = msm_outcome)
  
  msm_fit <- glm(
    as.formula(full_msm_formula),
    data = msm_data,
    family = msm_family
  )
  out <- list(
    msm_fit = msm_fit,
    which_ps_numerator = which_ps_numerator, 
    which_ps_denominator = which_ps_denominator
  )
  return(out)
}
