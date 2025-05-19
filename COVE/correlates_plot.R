# ------------------------------------------------------------------------------
# Plot correlates VE NAI
# ------------------------------------------------------------------------------

here::i_am("ve_nai/correlates_plot.R")

library(tidyverse)
library(wCorr)

source(here::here("ve_nai/ve_from_fits.R"))

# Read correlates data
# corr_data <- read.csv("/Volumes/trials/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv")
corr_data <- readRDS(here::here("raw_data/P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.Rds")) %>%
  rename('SUBJID' = Ptid) %>%
  select(-c("Age", "Sex", "HighRiskInd", "EthnicityHispanic", "EthnicityNotreported", "EthnicityUnknown", "Black", "Asian", "NatAmer", "PacIsl", "Multiracial"))

# Read in prepped analytic dataset
cove_data <- readRDS(here::here("analytic_data/final_asymp_data.Rds"))

# Merge on ID
data <- left_join(cove_data, corr_data, by = "SUBJID")

# Get VE FIT results ------------------------------------------------------------

symp_cox_fit <- survival::coxph(survival::Surv(ftime, ftype == 1) ~ vax + 
                                  Age + Sex + HighRiskInd + 
                                  EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
                                  Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + RaceNotreported + RaceUnknown, data = data)

asymp_cox_fit <- survival::coxph(survival::Surv(ftime, ftype == 2) ~ vax + 
                                   Age + Sex + HighRiskInd + 
                                   EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
                                   Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + RaceNotreported + RaceUnknown, data = data)

# Note - 168 people missing BMI and being dropped from models
# causing issues later on so leave out BMI for now

symp_lr_fit <- glm(
  symp_ind ~ vax + 
    Age + Sex + HighRiskInd + 
    EthnicityHispanic + EthnicityNotreported + EthnicityUnknown +
    Black + Asian + NatAmer + PacIsl + Multiracial + RaceOther + RaceNotreported + RaceUnknown, 
  # subset to infected ptcpt only
  data = data[data$ftype > 0, ],
  family = stats::binomial()
)

ve_fit <- ve_from_fits(
  symp_cox_fit = symp_cox_fit, 
  asymp_cox_fit = asymp_cox_fit, 
  symp_lr_fit = symp_lr_fit, 
  data = data, 
  t0 = c(175, 180, 185, 190, 195, 200)
)

id_df <- data %>%
  select(SUBJID, 
         vax, 
         symp_ind, 
         asymp_ind,
         Age,
         # Day29bindSpike,
         # Day29bindRBD,
         # Day29bindN,
         Day57bindSpike,
         Day57bindRBD,
         Day57bindN,
         Day57pseudoneutid50,
         Day29pseudoneutid50,
         wt.D29,
         wt.D57) %>%
  mutate(vax = factor(vax, levels = 0:1, labels = c("Placebo", "Vaccine")),
         case_type = if_else(symp_ind == 1, "Symptomatic", 
                             if_else(asymp_ind == 1, "Asymptomatic", "No case")),
         case_type = factor(case_type, levels = c("No case", "Asymptomatic", "Symptomatic")))

p_immune_X <- cbind(id_df, ve_fit$p_immune_X) %>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`) 

p_doomed_X <- cbind(id_df, ve_fit$p_doomed_X)%>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`) 

p_alwaysinf_X <- cbind(id_df, ve_fit$p_alwaysinf_X) %>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`)  

p_converted_X <- cbind(id_df, ve_fit$p_converted_X) %>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`)  

p_helped_X <- cbind(id_df, ve_fit$p_helped_X) %>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`) 

p_helpedplus_X <- cbind(id_df, ve_fit$p_helpedplus_X) %>%
  rename("t0_175" = `1`,
         "t0_180" = `2`,
         "t0_185" = `3`,
         "t0_190" = `4`,
         "t0_195" = `5`,
         "t0_200" = `6`)  


####################### Plotting function #######################
plot_correlates <- function(data, x_vars, y_vars, case_type_var, data_type = "vaccine", p_type = "P(Immune)", cor_method = "Spearman"){
  
  plot_list <- vector("list", length = length(x_vars) * length(y_vars))
  i <- 1
  
  if(data_type == "vaccine"){
    data <- data[data$vax == "Vaccine",]
  } else if(data_type == "Placebo"){
    data <- data[data$vax == "Placebo",]
  } 
  
  for(x_name in x_vars){
    for(y_name in y_vars){
      
      wt_var_name <- if (grepl("29", x_name)) "wt.D29" else "wt.D57"
      
      data_filtered <- data %>%
        filter(!is.na(.data[[x_name]])) %>%
        filter(!is.na(.data[[wt_var_name]]))
      
      p_cor <- wCorr::weightedCorr(
        data_filtered[[x_name]], 
        data_filtered[[y_name]], 
        method = cor_method, 
        weights = data_filtered[[wt_var_name]]
      )
      
      p_cor_plot <- ggplot(data = data_filtered, aes(x = .data[[x_name]], y = .data[[y_name]], color = .data[[case_type_var]])) +
        geom_point(alpha = 0.6) +
        labs(title = paste0(p_type, " vs ", x_name, " at ", y_name),
             subtitle = paste0("Weighted ", cor_method, " r = ", round(p_cor, 2)), 
             x = x_name,
             y = p_type,
             color = "Case status"
        ) + 
        theme_minimal()
      
      plot_list[[i]] <- p_cor_plot
      i <- i + 1
      
    }
  }
  
  return(plot_list)
  
}

###############################################################################


# Immune ---------------------------------------------------------------------

# x_vars <- c("t0_175", "t0_180", "t0_185", "t0_190", "t0_195", "t200")
# Just use 180 for now so I don't have 1000 plots??? could go back and generalize later
y_vars <- c("t0_180")

x_vars <- c("Day29pseudoneutid50") #,
            #"Day57bindSpike",
            #"Day57bindRBD",
            #"Day57bindN")


immune_plots <- plot_correlates(p_immune_X, 
                                x_vars = x_vars, 
                                y_vars = y_vars, 
                                case_type_var = "case_type",
                                data_type = "vaccine", 
                                p_type = "P(Immune)")

# Doomed ---------------------------------------------------------------------

doomed_plots <- plot_correlates(p_doomed_X, 
                                x_vars = x_vars, 
                                y_vars = y_vars, 
                                case_type_var = "case_type",
                                data_type = "vaccine", 
                                p_type = "P(Doomed)")

test_regression <- glm(as.formula("t0_180 ~ Day29pseudoneutid50 + Age"), data = p_doomed_X, weights = data$wt.D29)


# Always Infected ------------------------------------------------------------

always_plots <- plot_correlates(p_alwaysinf_X, 
                                x_vars = x_vars, 
                                y_vars = y_vars, 
                                case_type_var = "case_type",
                                data_type = "vaccine", 
                                p_type = "P(Always Infected)")

# Converted ------------------------------------------------------------------

converted_plots <- plot_correlates(p_converted_X,
                                x_vars = x_vars, 
                                y_vars = y_vars, 
                                case_type_var = "case_type",
                                data_type = "vaccine", 
                                p_type = "P(Converted)")

# Helped ---------------------------------------------------------------------

helped_plots <- plot_correlates(p_helped_X,
                                   x_vars = x_vars, 
                                   y_vars = y_vars, 
                                   case_type_var = "case_type",
                                   data_type = "vaccine", 
                                   p_type = "P(Helped)")

# Helped Plus ----------------------------------------------------------------

helpedplus_plots <- plot_correlates(p_helpedplus_X,
                                x_vars = x_vars, 
                                y_vars = y_vars, 
                                case_type_var = "case_type",
                                data_type = "vaccine", 
                                p_type = "P(Helped Plus)")
