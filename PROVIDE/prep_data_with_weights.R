# ---------------------------------------------------
# Script to make dataset for PROVIDE analysis
# ---------------------------------------------------

here::i_am("prep_data_with_weights.R")

library(tidyverse)
library(labelled)

# Per-protocol definition:
# A) Volunteers assigned to Rotarix arm were vaccinated per protocol
# B) At least 365 days of study surveillance 
# C) No rotavirus diarrhea prior to week 18
# D) Wk18 and Wk53 plasma samples available (≥20µl)

prep_provide_data <- function(){
  
  # Load data
  load(here::here("raw_data/PROVIDE diarrhea.RData"))
  load(here::here("raw_data/provide_f10_r1-3bk7-.RData"))
  
  # Year 1 outcome data:
  # select id, per-protocol flag, flag dropped from study <365 days, at least one episode of rota diarrhea,
  # number of missing diarrheal stools
  data <- mgmt_rto_rotatrial_outcome_f10 %>%
    select(sid,
           rotaarm,
           rotaprotd,
           drop365,
           rotaepi)
  
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
  
  # Growth (baseline, week 10, week 18, week 52, week 91)
  sub_growth <- danth %>%
    select(sid,
           vstwk,
           haz,
           ageday,
           gender) %>%
    filter(vstwk %in% c("Enrollment", "Week 10", "Week 18", "Week 52", "Week 53", "Week 91")) %>%
    pivot_wider(names_from = vstwk,
                values_from = c(haz, ageday),
                names_sep = c("_")) %>%
    rename(enr_haz = haz_Enrollment,
           wk10_haz = `haz_Week 10`,
           wk18_haz = `haz_Week 18`,
           wk52_haz = `haz_Week 52`,
           wk53_haz = `haz_Week 53`,
           wk91_haz = `haz_Week 91`,
           enr_ageday = ageday_Enrollment,
           wk10_ageday = `ageday_Week 10`,
           wk18_ageday = `ageday_Week 18`,
           wk52_ageday = `ageday_Week 52`,
           wk53_ageday = `ageday_Week 53`,
           wk91_ageday = `ageday_Week 91`
           ) %>%
    mutate(across(c(enr_haz, wk10_haz, wk18_haz, wk91_haz, wk52_haz, wk53_haz,
                    enr_ageday, wk10_ageday, wk18_ageday, wk91_ageday, wk52_ageday, wk53_ageday), ~ na_if(.x, -9))) # -9 = Missing
  
  data <- left_join(data, sub_growth, by = "sid")
  
  # Get enrollment date (to use later with asymptomatic data/derive ftime)
  enr_date <- danth %>%
    filter(vstwk == "Enrollment") %>%
    select(sid, vstdate) %>%
    rename('enr_date' = vstdate)
  
  data <- left_join(data, enr_date, by = 'sid')
  
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
  
  #### D) Wk18 and Wk53 plasma samples available (≥20µl) #######################
  
  # ^ this does not seem relevant for us but using as part of per-protocol def
  
  # Wk18 and Wk53 samples available
  wk18_plasma_avail <- bv_f18_week_18_visit_f10 %>%
    filter(f18bld == 1)
  
  wk53_plasma_avail <- bv_f53_week_53_visit_f10 %>%
    filter(f53bld == 1)
  
  data <- data %>%
    filter(sid %in% unique(wk18_plasma_avail$sid)) %>%
    filter(sid %in% unique(wk53_plasma_avail$sid))
  
  # Post D - N = 524 (264 placebo, 260 vax) (386 no rota, 138 rota)
  # table(data$rotaarm, data$rotaepi)
  # 0   1
  # 0 171  93
  # 1 215  45
  
  # Rotavirus - Week 91
  rota_ds_91 <- ddep %>%
    left_join(data, by = "sid") %>% # add in growth data
    filter(daysbirth < wk91_ageday & daysbirth >= wk18_ageday) %>% # episodes between week 18 and 91
    select(sid,
           epinum,
           epidate,
           lstvstwk,
           # epiexlbf,
           daysbirth,
           rotads,
           wk91_ageday
    ) %>%
    #mutate(across(c(daysbirth, rotads), ~ na_if(.x, -9))) %>%
    # Get if any rotads for each id
    group_by(sid) %>%
    mutate(
      # either day of first infection or day of week 91 visit (censoring)
      ftime91 = if (any(rotads == 1)) min(daysbirth[rotads == 1]) else wk91_ageday
    ) %>%
    summarise(
      any_rotads_elisa_wk91 = if_else(any(rotads == 1), 1, 0),
      non_rotads_elisa_wk91 = if_else(any(rotads == 1), 0, 1),
      ftime91 = first(ftime91)
    )
  
  # Rotavirus - Week 52
  rota_ds_52 <- ddep %>%
    left_join(data, by = "sid") %>% # add in growth data
    filter(daysbirth < wk52_ageday & daysbirth >= wk18_ageday) %>% 
    select(sid,
           epinum,
           epidate,
           lstvstwk,
           # epiexlbf,
           daysbirth,
           rotads,
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
      ftime52 = first(ftime52)
    )
  
  data <- left_join(data, rota_ds_52, by = "sid")
  data <- left_join(data, rota_ds_91, by = "sid")
  
  # If any NAs, assume no diarrhea; censoring time wk52 / wk91
  data$any_rotads_elisa_wk52 <- ifelse(is.na(data$any_rotads_elisa_wk52), 0, data$any_rotads_elisa_wk52)
  data$non_rotads_elisa_wk52 <- ifelse(is.na(data$non_rotads_elisa_wk52), 0, data$non_rotads_elisa_wk52)
  data$ftime52 <- ifelse(is.na(data$ftime52), data$wk52_ageday, data$ftime52)
  
  data$any_rotads_elisa_wk91 <- ifelse(is.na(data$any_rotads_elisa_wk91), 0, data$any_rotads_elisa_wk91)
  data$non_rotads_elisa_wk91 <- ifelse(is.na(data$non_rotads_elisa_wk91), 0, data$non_rotads_elisa_wk91)
  data$ftime91 <- ifelse(is.na(data$ftime91), data$wk91_ageday, data$ftime91)
  
  # Diarrhea TAC
  rota_ds_afe_wk91 <- diarrhea4 %>%
    left_join(data, by = "sid") %>%
    filter(daysbirth < wk91_ageday & daysbirth >= wk18_ageday) %>%  # - 28 ) %>% # filter only episodes that occured > 28 days before week 91 growth measurement
    select(sid,
           daysbirth,
           wk91_ageday,
           rotads,
           astrovirus_afe,
           norovirus_gii_afe,
           rotavirus_afe,
           sapovirus_afe,
           adenovirus_40_41_afe,
           aeromonas_afe,
           c_jejuni_coli_afe,
           cryptosporidium_afe,
           e_histolytica_afe,
           h_pylori_afe,
           salmonella_afe,
           shigella_eiec_afe,
           v_cholerae_afe,
           ST_ETEC_afe,
           TEPEC_afe,
           severe) %>% # note- assuming this is severity of diarrhea episode?? can't find dictionary for diarrhea4
    group_by(sid) %>%
    mutate(ftime91_tac = if (any(rotads == 1)) min(daysbirth[rotads == 1]) else wk91_ageday) %>%
    summarise(any_rotads_afe_wk91 = if_else(any(rotads == 1), 1, 0),
              ftime91_tac = first(ftime91_tac))
  
  data <- left_join(data, rota_ds_afe_wk91, by = c("sid"))
  
  rota_ds_afe_wk52 <- diarrhea4 %>%
    left_join(data, by = "sid") %>%
    filter(daysbirth < wk52_ageday & daysbirth >= wk18_ageday) %>% # - 28 ) %>% # filter only episodes that occured > 28 days before week 91 growth measurement
    select(sid,
           daysbirth,
           wk52_ageday,
           rotads,
           astrovirus_afe,
           norovirus_gii_afe,
           rotavirus_afe,
           sapovirus_afe,
           adenovirus_40_41_afe,
           aeromonas_afe,
           c_jejuni_coli_afe,
           cryptosporidium_afe,
           e_histolytica_afe,
           h_pylori_afe,
           salmonella_afe,
           shigella_eiec_afe,
           v_cholerae_afe,
           ST_ETEC_afe,
           TEPEC_afe,
           severe) %>% # note- assuming this is severity of diarrhea episode?? can't find dictionary for diarrhea4
    group_by(sid) %>%
    mutate(ftime52_tac = if (any(rotads == 1)) min(daysbirth[rotads == 1]) else wk52_ageday) %>%
    summarise(any_rotads_afe_wk52 = if_else(any(rotads == 1), 1, 0),
              ftime52_tac = first(ftime52_tac))
  
  data <- left_join(data, rota_ds_afe_wk52, by = c("sid"))
  
  # If any NAs, assume no diarrhea; censoring time week 91/52
  data$any_rotads_afe_wk91 <- ifelse(is.na(data$any_rotads_afe_wk91), 0, data$any_rotads_afe_wk91)
  data$ftime91_tac <- ifelse(is.na(data$ftime91_tac), data$wk91_ageday, data$ftime91_tac)
  
  data$any_rotads_afe_wk52 <- ifelse(is.na(data$any_rotads_afe_wk52), 0, data$any_rotads_afe_wk52)
  data$ftime52_tac <- ifelse(is.na(data$ftime52_tac), data$wk52_ageday, data$ftime52_tac)
  
  # Antibiotic use
  abx_key <- c(
    'azithromycin',
    'erythromycin',
    'metronidazole',
    'ciprofloxacin',
    'nalidixic acid',
    'pivmecillinam',
    'cephradine',
    'co-trimoxazole',
    'amoxycillin',
    'cefixime',
    'ceftazidime',
    'cefuroxime',
    'ceftriaxone',
    'fluclocaxillin',
    'other'
  )
  
  # get birthday to merge in
  dob <- bv_enc_infant_enrollment_f10 %>%
    select(sid, dob)
  
  abx_wk91 <- bv_dep_diarrheal_episode_f10 %>%
    left_join(dob, by = "sid") %>%
    left_join(data, by = "sid") %>%
    mutate(dedt = as.Date(dedt),
           dob = as.Date(dob),
           agedays = dedt - dob) %>%
    filter(agedays < wk91_ageday & agedays >= wk18_ageday) %>%  # - 28) %>%
    select(sid, 
           dewt,
           dedt,
           antb,
           antb1, 
           antb2) %>% # note nobody has antb3
    mutate(antb = factor(antb, levels = 1:3, labels = c("None", "Pre-visit", "At-visit")),
           antb1 = if_else(antb1 == 99, NA, antb1),
           antb2 = if_else(antb2 == 99, NA, antb2), 
           antb1 = factor(antb1, levels = 1:15, labels = abx_key),
           antb2 = factor(antb2, levels = 1:15, labels = abx_key)) %>%
    group_by(sid) %>%
    summarise(any_abx_wk91 = if_else(any(antb == "Pre-visit" | antb == "At-visit"), 1, 0))
  
  data <- left_join(data, abx_wk91, by = "sid")
  
  abx_wk52 <- bv_dep_diarrheal_episode_f10 %>%
    left_join(dob, by = "sid") %>%
    left_join(data, by = "sid") %>%
    mutate(dedt = as.Date(dedt),
           dob = as.Date(dob),
           agedays = dedt - dob) %>%
    filter(agedays < wk52_ageday & agedays >= wk18_ageday) %>% # - 28) %>%
    select(sid, 
           dewt,
           dedt,
           antb,
           antb1, 
           antb2) %>% # note nobody has antb3
    mutate(antb = factor(antb, levels = 1:3, labels = c("None", "Pre-visit", "At-visit")),
           antb1 = if_else(antb1 == 99, NA, antb1),
           antb2 = if_else(antb2 == 99, NA, antb2), 
           antb1 = factor(antb1, levels = 1:15, labels = abx_key),
           antb2 = factor(antb2, levels = 1:15, labels = abx_key)) %>%
    group_by(sid) %>%
    summarise(any_abx_wk52 = if_else(any(antb == "Pre-visit" | antb == "At-visit"), 1, 0))
  
  data <- left_join(data, abx_wk52, by = "sid")
  
  # Assume if any_abx = NA, did not get abx
  data$any_abx_wk91 <- ifelse(is.na(data$any_abx_wk91), 0 , data$any_abx_wk91)
  data$any_abx_wk52 <- ifelse(is.na(data$any_abx_wk52), 0 , data$any_abx_wk52)
  
  # NOTE bv_ptreat_f10 is any antibiotic or other med for any indication, this above im assuming is antibiotics for diarrhea specifically
  
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
  
  data <- left_join(data, sub_ses, by = "sid")
  
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
                        wk91_haz = "Week 91 HAZ",
                        wk52_haz = "Week 52 HAZ",
                        enr_date = "Enrollment date",
                        enr_ageday = "Enrollment age days",
                        wk91_ageday = "Week 91 age days",
                        wk52_ageday = "Week 52 age days",
                        non_rotads_elisa_wk91 = "Dirrhea episode not attributed to rotavirus ELISA - week 91",
                        any_rotads_elisa_wk91 = "Any rotavirus diarrhea episode ELISA - week 91",
                        any_rotads_afe_wk91 = "Any rotavirus diarrhea episode TAC - week 91",
                        any_abx_wk91 = "Antibiotics for any diarrhea episode - week 91",
                        non_rotads_elisa_wk52 = "Dirrhea episode not attributed to rotavirus ELISA - week 52",
                        any_rotads_elisa_wk52 = "Any rotavirus diarrhea episode ELISA - week 52",
                        any_rotads_afe_wk52 = "Any rotavirus diarrhea episode TAC - week 52",
                        any_abx_wk52 = "Antibiotics for any diarrhea episode - week 52",
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
                        food_avail_bin = "Food availability (=0 deficit, =1 non-deficit)")
  
  return(data)
  
}

data <- prep_provide_data()

saveRDS(data, file = here::here("analysis_data/per_protocol_data.Rds"))

# "PER-PROTOCOL" NUMBERS:

# 264 placebo, 260 vaccine
# 386 controls, 138 cases

#        control.  case
# placebo 171     93
# vaccine 215.    45

# -----------------------------------------------------------------------------

# down selected sample (N = 405 in word doc but really N = 407): 
asymp <- readxl::read_xlsx(here::here("raw_data/Copy of Merged Data_with dates_asymptomatic RV.xlsx"), sheet = "merged")

# asymp summary dataset with any asymp infection after week 18 
asymp_18_52 <- asymp %>%
  filter(Week >= 18) %>%
  group_by(SID) %>%
  mutate(
    fdate_less40 = if (any(Ct_less40 == 1)) min(Visit_Date[Ct_less40 == 1]) else as.POSIXct(NA),
    fdate_less36 = if (any(Ct_less36 == 1)) min(Visit_Date[Ct_less36 == 1]) else as.POSIXct(NA),
    fdate_less34 = if (any(Ct_less34 == 1)) min(Visit_Date[Ct_less34 == 1]) else as.POSIXct(NA)
  ) %>%
  summarise(any_asymp_Ct_less40 = if_else(any(Ct_less40 == 1), 1 , 0),
            any_asymp_Ct_less36 = if_else(any(Ct_less36 == 1), 1 , 0),
            any_asymp_Ct_less34 = if_else(any(Ct_less34 == 1), 1 , 0),
            fdate_less40 = fdate_less40[1],
            fdate_less36 = fdate_less36[1],
            fdate_less34 = fdate_less34[1]) %>%
  ungroup() %>%
  rename("sid" = SID)
# this is 1:1 when filtering to episodes week 18 to 52 (N = 270)

# table(asymp_filtered_merged$rotaarm, asymp_filtered_merged$rotaepi)
# 
#         control case
# placebo  90     90
# vaccine  45     45

# -----------------------------------------------------------------------------

# Make weights based on case status, treatment arm, sex
# (ignore middle part of word doc w/ control down selecting based on completeness of specimen collection, assume missing at random)

sampled_ids <- asymp_18_52$sid

data$sampled <- ifelse(data$sid %in% sampled_ids, 1, 0)

fit <- glm(sampled ~ rotaepi*rotaarm*gender, data = data, family = binomial())

data$weight <- 1 / fit$fitted.values

asymp_filtered_merged <- left_join(asymp_18_52, data, by = "sid")

# make asymptomatic ftime variable
asymp_filtered_merged <- asymp_filtered_merged %>%
  mutate(ftime_asymp_40 = as.integer(difftime(fdate_less40, enr_date)),
         ftime_asymp_36 = as.integer(difftime(fdate_less36, enr_date)),
         ftime_asymp_34 = as.integer(difftime(fdate_less34, enr_date)))

# Add single variable for asymptomatic vs symptomatic vs no infection and ftime
final_data <- asymp_filtered_merged %>%
  mutate(
    all_inf_40 = case_when(
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less40 == 0) ~ 0,        # No infection
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less40 == 1) ~ 1,        # Only asymptomatic
      (any_rotads_elisa_wk52 == 1) & (any_asymp_Ct_less40 == 0) ~ 2,        # Only symptomatic
      ftime_asymp_40 < ftime52 ~ 1,                                         # Both: asymptomatic occurred first
      ftime52 <= ftime_asymp_40 ~ 2                                         # Both: symptomatic occurred first or same day
    ),
    ftime_all_40 = ifelse(all_inf_40 == 0, ftime52, # censoring age (week 52 visit time) stored in ftime52 if no infection
                          ifelse(all_inf_40 == 1, ftime_asymp_40, ftime52)), 
    all_inf_36 = case_when(
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less36 == 0) ~ 0,        # No infection
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less36 == 1) ~ 1,        # Only asymptomatic
      (any_rotads_elisa_wk52 == 1) & (any_asymp_Ct_less36 == 0) ~ 2,        # Only symptomatic
      ftime_asymp_36 < ftime52 ~ 1,                                         # Both: asymptomatic occurred first
      ftime52 <= ftime_asymp_36 ~ 2                                         # Both: symptomatic occurred first or same day
    ),
    ftime_all_36 = ifelse(all_inf_36 == 0, ftime52, # censoring age (week 52 visit time) stored in ftime52 if no infection
                          ifelse(all_inf_36 == 1, ftime_asymp_36, ftime52)), 
    all_inf_34 = case_when(
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less34 == 0) ~ 0,        # No infection
      (any_rotads_elisa_wk52 == 0) & (any_asymp_Ct_less34 == 1) ~ 1,        # Only asymptomatic
      (any_rotads_elisa_wk52 == 1) & (any_asymp_Ct_less34 == 0) ~ 2,        # Only symptomatic
      ftime_asymp_34 < ftime52 ~ 1,                                         # Both: asymptomatic occurred first
      ftime52 <= ftime_asymp_34 ~ 2                                         # Both: symptomatic occurred first or same day
    ),
    ftime_all_34 = ifelse(all_inf_34 == 0, ftime52, # censoring age (week 52 visit time) stored in ftime52 if no infection
                          ifelse(all_inf_34 == 1, ftime_asymp_34, ftime52))
  )

saveRDS(final_data, file = here::here("analysis_data/ve_nai_provide.Rds"))

# --------------------------------------------------------------------------------

# test weighted VE

cox_model <- survival::coxph(survival::Surv(ftime52, rotaepi == 1) ~ rotaarm, 
                             data = final_data,
                             weights = final_data$weight,
                             ties = "efron")

ve <- 1 - exp(coef(cox_model))
ci <- 1 - exp(confint(cox_model))

# VE 45.6% (14.4%, 65.4%)
# ^^ THIS WAS BEFORE I ADDED CENSORING TIMES INTO FTIME52 (was NA before/getting dropped)
# now 56.9 (32.2 to 72.6) -- still close is

# The paper says they used formula (ARunvax - ARvax) / ARunvax * 100
# Wilson 95% CIs

AR_unvax <- sum(data$rotaepi[data$rotaarm == 0]) / nrow(data[data$rotaarm == 0,])
AR_vax <- sum(data$rotaepi[data$rotaarm == 1]) / nrow(data[data$rotaarm == 1,])
VE_fromAR <- (AR_unvax - AR_vax) / (AR_unvax) * 100

# VE_fromAR = 50.9
# paper reports 51.0 (33.8-63.7) but idk how they did the CI?? 


# "Vaccine efficacy was defined as 1 minus the hazard ratio (mRNA-1273 vs. placebo), 
# and 95% confidence intervals were estimated using a stratified Cox proportional-hazards model
# with Efron’s method of tie handling and with the treatment group as a covariate, 
# adjusted for stratification factor."

# case <- c("SARS-CoV-2 Infection Regardless of Symptomatology or Severity at PDV",
#           "SARS-CoV-2 Infection Regardless of Symptomatology or Severity in Blinded Phase")
# 
# analysis_data <- any_inf_data %>%
#   mutate(ftype = if_else(EVNTDESC %in% case, 1, 0),
#          vax = if_else(TRT01P == "mRNA-1273", 1, 0)) %>%
#   select(ftype,
#          vax,
#          AVAL,
#          STRATAR) %>%
#   rename("ftime" = AVAL,
#          "strat_rand" = STRATAR)
# 
# cox_model <- survival::coxph(survival::Surv(ftime, ftype == 1 ) ~ vax + strata(strat_rand), 
#                              data = analysis_data, 
#                              ties = "efron")
# 
# summary(cox_model)
# 
# ve <- 1 - exp(coef(cox_model))
# ci <- 1 - exp(confint(cox_model))
