
here::i_am("PROVIDE/prep_data_with_weights.R")

library(tidyverse)
library(labelled)

# Per-protocol definition (in word doc):
# A) Volunteers assigned to Rotarix arm were vaccinated per protocol
# B) At least 365 days of study surveillance 
# C) No rotavirus diarrhea prior to week 18
# D) Wk18 and Wk53 plasma samples available (≥20µl)

# ^ this is what we used for asymptomatic data

# Per-protocol definition in severity paper says:
# "Secondary per-protocol (PP) analyses were performed including all children who had 365 days of
# follow-up and, if assigned to the vaccine arm, received both doses
# of rotavirus vaccine within the protocol-specified window."

# So just A & B from the above?
# Except we're also interested in it being post-vax so still going to do C

prep_provide_data_severe <- function(){
  
  # Load data
  load(here::here("PROVIDE/raw_data/PROVIDE diarrhea.RData"))
  load(here::here("PROVIDE/raw_data/provide_f10_r1-3bk7-.RData"))
  
  # Year 1 outcome data:
  # select id, per-protocol flag, flag dropped from study <365 days, at least one episode of rota diarrhea,
  # number of missing diarrheal stools
  data <- mgmt_rto_rotatrial_outcome_f10 %>%
    select(sid,
           rotaarm,
           rotaprotd,
           drop365,
           rotaepi,
           nsvrtepi)
  
  #### A) Volunteers assigned to Rotarix arm were vaccinated per protocol ######
  data <- data %>%
    filter(rotaprotd == "PP")
  
  # Post A - N = 678 (350 placebo, 328 vax) (499 no rota, 179 rota)
  
  # table(data$rotaarm, data$rotaepi)
  # 0   1
  # 0 236 114
  # 1 263  65
  
  #### B) At least 365 days of study surveillance ##############################
  data <- data %>%
    filter(drop365 == 0)
  
  # Post B - N = 596 (304 placebo, 292 vax) (420 no rota, 176 rota)
  
  # table(data$rotaarm, data$rotaepi)
  # 
  # 0   1
  # 0 193 111
  # 1 227  65
  
  # NOTE: Paper reports N = 593 - we have an extra 3 placebo
  
  # Get enrollment date (to use later to derive ftime)
  enr_date <- danth %>%
    filter(vstwk == "Enrollment") %>%
    select(sid, vstdate) %>%
    rename('enr_date' = vstdate)
  
  data <- left_join(data, enr_date, by = 'sid')
  
  # Growth & visit days (baseline, week 10, week 18, week 52, week 91)
  sub_growth <- danth %>%
    select(sid,
           vstwk,
           haz,
           ageday,
           gender) %>%
    filter(vstwk %in% c("Enrollment", "Week 10", "Week 18", "Week 52")) %>%
    pivot_wider(names_from = vstwk,
                values_from = c(haz, ageday),
                names_sep = c("_")) %>%
    rename(enr_haz = haz_Enrollment,
           wk10_haz = `haz_Week 10`,
           wk18_haz = `haz_Week 18`,
           wk52_haz = `haz_Week 52`,
           enr_ageday = ageday_Enrollment,
           wk10_ageday = `ageday_Week 10`,
           wk18_ageday = `ageday_Week 18`,
           wk52_ageday = `ageday_Week 52`,
    ) %>%
    mutate(across(c(enr_haz, wk10_haz, wk18_haz, wk52_haz, 
                    enr_ageday, wk10_ageday, wk18_ageday, wk52_ageday), ~ na_if(.x, -9))) # -9 = Missing

   data <- left_join(data, sub_growth, by = "sid")
  
  #### C) No rotavirus diarrhea prior to week 18 ###############################
  # ddep is dataset of diarrhea episodes
  # identify episodes before week 18 (need to get wk 18 age from growth data) that are positive for rota & remove those ids from main data
  early <- ddep %>%
    left_join(data, by = "sid") %>%
    filter(daysbirth < wk18_ageday) %>%
    filter(rotads == 1)

  data <- data %>%
    filter(!(sid %in% unique(early$sid)))
  
  # Post C - N = 569 (295 placebo, 274 vax) (420 no rota, 149 rota)
  
  #table(data$rotaarm, data$rotaepi)
  # 0   1
  # 0 193 102
  # 1 227  47
  
  # Rotavirus - Week 52
  rota_ds_52 <- ddep %>%
    left_join(data, by = "sid") %>% # add in growth data
    filter(daysbirth < wk52_ageday) %>% 
    select(sid,
           epinum,
           enr_date,
           epidate,
           daysenc,
           lstvstwk,
           daysbirth,
           rotads,
           ruuska,
           wk52_ageday
    ) %>%
    #mutate(across(c(daysbirth, rotads), ~ na_if(.x, -9))) %>%
    # Get if any rotads for each id
    group_by(sid) %>%
    mutate(
      # either day of first infection or day of week 52 visit (censoring)
      ftime52 = if (any(rotads == 1)) min(daysbirth[rotads == 1]) else wk52_ageday
    ) %>%
    summarise(
      any_rotads_elisa_wk52 = if_else(any(rotads == 1), 1, 0),
      non_rotads_elisa_wk52 = if_else(any(rotads == 1), 0, 1),
      ftime52 = first(ftime52),
      sev_score = first(ruuska[which(daysbirth == first(ftime52))])
    ) %>%
    mutate(
      severe_RVD = if_else(is.na(sev_score), 0, # If NA, did not have any rotavirus (censored at week 52 visit)
                           if_else(sev_score >= 11, 1, 0)) # If else, check their severity to see if severe diarrhea. If >= 11 severe, if <11 (including -9 missing) assume not severe
    )
  
  data <- left_join(data, rota_ds_52, by = "sid")
  
  # If any NAs, assume no diarrhea; censoring time wk52
  data$any_rotads_elisa_wk52 <- ifelse(is.na(data$any_rotads_elisa_wk52), 0, data$any_rotads_elisa_wk52)
  data$non_rotads_elisa_wk52 <- ifelse(is.na(data$non_rotads_elisa_wk52), 0, data$non_rotads_elisa_wk52)
  data$severe_RVD <- ifelse(is.na(data$severe_RVD), 0, data$severe_RVD)
  data$ftime52 <- ifelse(is.na(data$ftime52), data$wk52_ageday, data$ftime52)
  
  # Make final variable for infection level and time
  # time_name = "ftime_all_30",
  # symp_ind_name = "symp_ind_30",
  # event_name = "all_inf_30",
  # symp_level = 2,
  # asymp_level = 1,
  
  data$inf_level <- ifelse(data$severe_RVD == 1, 2,
                           ifelse(data$any_rotads_elisa_wk52 == 1, 1, 0))
  
  # ---------------------------------------------------------------------------
  
  # Baseline SES related covariates
  
  # Number of siblings <5years
  # Number of people sleeping in household
  # Mother education
  # Father education
  # Total monthly income
  # Electricity
  # Cooking gas
  # TV
  # Drinking water source
  # Toilet
  # Food availability
  
  # NEW:
  # continuous IgA @ week 6
  # binary IgA < 20 @ week 6
  
  sub_ses <- bv_ses_water_f10 %>%
    select(sid,
           sb5y,
           pepl,
           fedu,
           medu,
           inco,
           elec,
           gas,
           tv,
           watr,
           toil,
           hhfd)
  
  sub_iga <- lab_bv_rpliga_f10 %>%
    filter(vstnum == 2) %>%   # filter to week 6 visit (2 = week 6)
    select(sid, uperml2) %>%  # select sid and endpoint titer highlighted for analysis
    mutate(iga_cont = uperml2,
           iga_bin = if_else(iga_cont >= 20, 1, 0)) %>%
    select(-uperml2)
  
  data <- left_join(data, sub_ses, by = "sid") %>%
    left_join(sub_iga, by = "sid")
  
  # Make missingness var for IgA
  data <- data %>%
    mutate(I_iga_measured = if_else(is.na(iga_cont), 0, 1),
           I_iga_measured_x_cont = if_else(is.na(iga_cont), 0, iga_cont),
           I_iga_measured_x_bin = if_else(is.na(iga_bin), 0, iga_bin))
  
  # Transform factors 
  data$gender <- factor(data$gender)
  data$fedu <- factor(data$fedu, levels = 1:18, labels = c("No formal education",
                                                           "1st year",
                                                           "2nd year",
                                                           "3rd year",
                                                           "4th year",
                                                           "5th year",
                                                           "6th year",
                                                           "7th year",
                                                           "8th year",
                                                           "9th year",
                                                           "SSC_Dakhil passed",
                                                           "HSC_Fazil passed",
                                                           "Vocational_diploma_homeopathy_LMF_etc",
                                                           "Degree_Alim passed",
                                                           "Hons passed_3 or 4yr hons",
                                                           "Master_Kamil passed",
                                                           "MBBS_MD_FCPS_FRCP",
                                                           "BSC Engineer_MSc_PhD_etc"))
  data$medu <- factor(data$medu, levels = 1:18, labels = c("No formal education",
                                                           "1st year",
                                                           "2nd year",
                                                           "3rd year",
                                                           "4th year",
                                                           "5th year",
                                                           "6th year",
                                                           "7th year",
                                                           "8th year",
                                                           "9th year",
                                                           "SSC_Dakhil passed",
                                                           "HSC_Fazil passed",
                                                           "Vocational_diploma_homeopathy_LMF_etc",
                                                           "Degree_Alim passed",
                                                           "Hons passed_3 or 4yr hons",
                                                           "Master_Kamil passed",
                                                           "MBBS_MD_FCPS_FRCP",
                                                           "BSC Engineer_MSc_PhD_etc"))
  
  # Make edu vars with fewer categories
  # No formal education
  # Primary school = 1st - 6th
  # Some secondary = 7th - 9th
  # Dakhil = secondary school? / Degree Alim = upper secondary
  # Fazil = bachelors? add Hons to that too?
  # Vocational
  # MBBS = medical degree, Engineer, Kamil? = higher than bachelors? 
  data <- data %>%
    mutate(
      medu_cat = case_when(
        medu %in% c("No formal education") ~ 1,
        medu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year") ~ 2,
        medu %in% c("7th year", "8th year", "9th year") ~ 3,
        medu %in% c("SSC_Dakhil passed", "Degree_Alim passed") ~ 4,
        medu %in% c("Vocational_diploma_homeopathy_LMF_etc") ~ 5,
        medu %in% c("HSC_Fazil passed", "Hons passed_3 or 4yr hons") ~ 6,
        medu %in% c("Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 7,
        TRUE ~ NA_real_
      ),
      medu_cat2 = case_when(
        medu %in% c("No formal education") ~ 1,
        medu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year", 
                    "7th year", "8th year", "9th year") ~ 2,
        medu %in% c("SSC_Dakhil passed", "Degree_Alim passed") ~ 3,
        medu %in% c("Vocational_diploma_homeopathy_LMF_etc") ~ 4,
        medu %in% c("HSC_Fazil passed", "Hons passed_3 or 4yr hons", 
                    "Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 5,
        TRUE ~ NA_real_
      ),
      medu_bin = case_when(
        medu %in% c("No formal education") ~ 0,
        medu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year", 
                    "7th year", "8th year", "9th year",
                    "SSC_Dakhil passed", "Degree_Alim passed",
                    "Vocational_diploma_homeopathy_LMF_etc",
                    "HSC_Fazil passed", "Hons passed_3 or 4yr hons", 
                    "Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 1,
        TRUE ~ NA_real_
      ),
      fedu_cat = case_when(
        fedu %in% c("No formal education") ~ 1,
        fedu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year") ~ 2,
        fedu %in% c("7th year", "8th year", "9th year") ~ 3,
        fedu %in% c("SSC_Dakhil passed", "Degree_Alim passed") ~ 4,
        fedu %in% c("Vocational_diploma_homeopathy_LMF_etc") ~ 5,
        fedu %in% c("HSC_Fazil passed", "Hons passed_3 or 4yr hons") ~ 6,
        fedu %in% c("Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 7,
        TRUE ~ NA_real_
      ),
      fedu_cat2 = case_when(
        fedu %in% c("No formal education") ~ 1,
        fedu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year", 
                    "7th year", "8th year", "9th year") ~ 2,
        fedu %in% c("SSC_Dakhil passed", "Degree_Alim passed") ~ 3,
        fedu %in% c("Vocational_diploma_homeopathy_LMF_etc") ~ 4,
        fedu %in% c("HSC_Fazil passed", "Hons passed_3 or 4yr hons", 
                    "Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 5,
        TRUE ~ NA_real_
      ),
      fedu_bin = case_when(
        fedu %in% c("No formal education") ~ 0,
        fedu %in% c("1st year", "2nd year", "3rd year", 
                    "4th year", "5th year", "6th year", 
                    "7th year", "8th year", "9th year",
                    "SSC_Dakhil passed", "Degree_Alim passed",
                    "Vocational_diploma_homeopathy_LMF_etc",
                    "HSC_Fazil passed", "Hons passed_3 or 4yr hons", 
                    "Master_Kamil passed", "MBBS_MD_FCPS_FRCP", "BSC Engineer_MSc_PhD_etc") ~ 1,
        TRUE ~ NA_real_
      ),
      medu_cat = factor(medu_cat, levels = 1:7, labels = c("No formal education",
                                                           "Some primary school",
                                                           "Some secondary school",
                                                           "High school",
                                                           "Vocational degree",
                                                           "Bachelors degree",
                                                           "Graduate degree")),
      fedu_cat = factor(fedu_cat, levels = 1:7, labels = c("No formal education",
                                                           "Some primary school",
                                                           "Some secondary school",
                                                           "High school",
                                                           "Vocational degree",
                                                           "Bachelors degree",
                                                           "Graduate degree")),
      medu_cat2 = factor(medu_cat2, levels = 1:5, labels = c("No formal education",
                                                             "Less than high school",
                                                             "High school diploma",
                                                             "Vocational degree",
                                                             "College degree")),
      fedu_cat2 = factor(fedu_cat2, levels = 1:5, labels = c("No formal education",
                                                             "Less than high school",
                                                             "High school diploma",
                                                             "Vocational degree",
                                                             "College degree")),
    )
  
  data$elec <- ifelse(data$elec == 1, 1, 0)
  data$gas <- ifelse(data$gas == 1, 1, 0)
  data$tv <- ifelse(data$tv == 1, 1, 0)
  
  data$watr <- factor(data$watr, levels = 1:3, labels = c("Municipality supply_piped",
                                                          "Own arrangement by pump",
                                                          "Tube well"))
  
  data$watr_bin <- ifelse(data$watr == "Municipality supply_piped", 1, 0)
  
  data$toil <- factor(data$toil, levels = 1:5, labels = c("Septic tank or toilet",
                                                          "Water sealed or slap latrine",
                                                          "Pit latrine",
                                                          "Open latrine",
                                                          "Hanging latrine"))
  
  data$toil_bin <- ifelse(data$toil %in% c("Septic tank or toilet",
                                           "Water sealed or slap latrine"), 1, 0)
  
  data$hhfd <- factor(data$hhfd, levels = 1:4, labels = c("Deficit whole year",
                                                          "Sometimes deficit",
                                                          "Neither deficit nor surplus",
                                                          "Surplus"))
  
  data$hhfd_bin <- ifelse(data$hhfd %in% c("Neither deficit nor surplus",
                                           "Surplus"), 1, 0)
  
  data <- data %>%
    rename('num_hh_lt_5' = sb5y,
           'num_hh_sleep' = pepl,
           'food_avail' = hhfd,
           'food_avail_bin' = hhfd_bin) %>%
    set_variable_labels(sid = "Participant ID",
                        rotaarm = "Rotavirus vaccine",
                        gender = "Sex",
                        enr_haz = "Baseline HAZ",
                        wk52_haz = "Week 52 HAZ",
                        enr_date = "Enrollment date",
                        enr_ageday = "Enrollment age days",
                        wk52_ageday = "Week 52 age days",
                        non_rotads_elisa_wk52 = "Dirrhea episode not attributed to rotavirus ELISA - week 52",
                        any_rotads_elisa_wk52 = "Any rotavirus diarrhea episode ELISA - week 52",
                        sev_score = "RUUSKA severity score",
                        severe_RVD = "Severe rotavirus (RUUSKA >= 11)",
                        inf_level = "Indicator of rotavirus & severity (2=severe, 1=mild, 0=no RVD)",
                        num_hh_lt_5 = "Number of siblings under 5 years old",
                        num_hh_sleep = "Number of people usually sleeping in household",
                        fedu = "Father's education level",
                        fedu_cat = "Father's education level - 7 categories",
                        fedu_cat2 = "Father's education level - 5 categories",
                        fedu_bin = "Father's education level binary (No education = 0, Any formal education = 1)",
                        medu = "Mother's education level",
                        medu_cat = "Mother's education level - 7 categories",
                        medu_cat2 = "Mother's education level - 5 categories",
                        fedu_bin = "Mother's education level binary (No education = 0, Any formal education = 1)",
                        inco = "Household income",
                        elec = "Household electricity",
                        gas = "Household gas",
                        tv = "Household TV",
                        watr = "Source of drinking water",
                        watr_bin = "Municipality piped water (=1 yes, =0 no)",
                        toil = "Primary toilet facility",
                        toil_bin = "Improved toilet (=1 yes, =0 no)",
                        food_avail = "Household food availability",
                        food_avail_bin = "Food availability (=0 deficit, =1 non-deficit)",
                        iga_cont = "IgA endpoint titer U/mL",
                        iga_bin = "I(IgA endpoint titer >= 20U/mL)",
                        I_iga_measured = "Indicator IgA measured (non-missing)",
                        I_iga_measured_x_cont = "Indicator not missing * continuous",
                        I_iga_measured_x_bin = "Indicator not missing * binary")
  
  return(data)
  
}

data <- prep_provide_data_severe()

saveRDS(data, file = here::here("PROVIDE/analysis_data/provide_severity_pp.Rds"))
