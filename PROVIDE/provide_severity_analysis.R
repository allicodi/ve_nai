# ---------------------------------------------------------------------------
# VE NAI analysis with severity 
# ---------------------------------------------------------------------------

here::i_am("PROVIDE/provide_severity_analysis.R")

devtools::load_all("../VEnai/")

data <- readRDS(here::here("PROVIDE/analysis_data/provide_severity_pp.Rds"))

# drop obs missing wk10haz
data <- data[-which(is.na(data$wk10_haz)),]

results <-  ve_nai(data = data,
                   vax_name = "rotaarm",
                   time_name = "ftime52",
                   symp_ind_name = "severe_RVD", # == 1 if severe
                   event_name = "inf_level",  # 2 = severe, 1 = mild, 0 = none
                   symp_level = 2,   # symptomatic == severe
                   asymp_level = 1,  # asymptomatic == mild
                   covariate_names = c("wk10_haz", "num_hh_lt_5", "toil_bin",
                                       "gender", "num_hh_sleep", "fedu_bin", 
                                       "medu_bin", "inco", "gas", "tv",
                                       "food_avail_bin"),
                   t0 = 365,
                   n_boot = 1000)
