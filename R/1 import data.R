
# extract baseline data --------------------------------------------------------


data_baseline <-
  full_join(
    readxl::read_excel(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx"), sheet = 1, guess_max = 20000, na = c("-9", "-1", "NA")),
    readxl::read_excel(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx"), sheet = 2, guess_max = 20000, na = c("-9", "-1", "NA")),
    by = "StudyID"
  ) %>%
  transmute(
    # admin
    ptid = StudyID,
    fup_month = 0,
    # hospid = `Hospital number`,
    date = ymd(`Study start date`),
    date_birth = ymd(`Date of Birth`),

    # demographics
    # sex = case_when(Gender==1 ~ "male", Gender==2 ~ "female", TRUE ~ NA_character_),
    sex = factor(Gender, levels = 1:2, labels = c("male", "female")),
    ethnicity = factor(`Ethnic group`, levels = 1:5, labels = c("white", "asian", "black", "chinese", "other")),
    # physical
    height = Height0,
    weight = Wgt0,

    SBP = SBP01,
    DBP = DBP01,

    PRD = `Coded/Verified PrimaryRenalDisease`,

    # lifestyle

    karnofsky = `Karnofsky performance`,
    alcohol_intake = `Units of alcohol per week`,
    smoking = `Smoking status`,
    smoking_cigs = `Cigarettes per day`,
    smoking_years = `NoOfYrsSmoked`
  ) %>%
  mutate(
    height = case_when(
      as.numeric(height) > 10 ~ as.numeric(height) / 100,
      height == "4.11" ~ 4 * 0.3048 + 11 * 0.0254,
      height == "5" ~ 5 * 0.3048,
      height == "5.1" ~ 5 * 0.3048 + 1 * 0.0254,
      height == "5.2" ~ 5 * 0.3048 + 2 * 0.0254,
      height == "5.3" ~ 5 * 0.3048 + 3 * 0.0254,
      height == "5.35" ~ 5 * 0.3048 + 3.5 * 0.0254,
      height == "5.4" ~ 5 * 0.3048 + 4 * 0.0254,
      height == "5.5" ~ 5 * 0.3048 + 5 * 0.0254,
      height == "5.6" ~ 5 * 0.3048 + 6 * 0.0254,
      height == "5.7" ~ 5 * 0.3048 + 7 * 0.0254,
      height == "5.75" ~ 5 * 0.3048 + 7.5 * 0.0254,
      height == "5.8" ~ 5 * 0.3048 + 8 * 0.0254,
      height == "5.9" ~ 5 * 0.3048 + 9 * 0.0254,
      height == "5.10" ~ 5 * 0.3048 + 10 * 0.0254,
      height == "5.11" ~ 5 * 0.3048 + 11 * 0.0254,
      height == "5.105" ~ 5 * 0.3048 + 10.5 * 0.0254,
      height == "6" ~ 6 * 0.3048,
      height == "6.1" ~ 6 * 0.3048 + +1 * 0.0254,
      height == "6.35" ~ 6 * 0.3048 + +3.5 * 0.0254,
      # height==5.1 ~ 5*0.3048 + 10*0.0254 OR height==5.1 ~ 5*0.3048 + 1*0.0254
      TRUE ~ as.numeric(height)
    ),

    date_birth =
      case_when(
        ptid == 696 & date_birth == ymd("2950-12-19") ~ ymd("1950-12-19"),
        TRUE ~ date_birth
      ),

    date = case_when(
      ptid == 2497 & date == ymd("2103-02-06") ~ ymd("2013-02-06"),
      ptid == 3036 & date == ymd("1929-12-10") ~ ymd("2014-12-10"),
      TRUE ~ date
    ),
  ) %>%
  filter(!is.na(ptid))



data_physical <- read_csv(file = here::here("data-raw", "james", "data_weight_BP.csv")) %>%
  rename(fup_month = fup) %>%
  mutate(
    weight = weight %>% na_if(-9) %>% na_if(-1),
    SBP = SBP1 %>% na_if(-9) %>% na_if(-1),
    DBP = DBP1 %>% na_if(-9) %>% na_if(-1),
  ) %>%
  select(-SBP1, -DBP1, -SBP2, -DBP2) %>%
  bind_rows(data_baseline %>% select(ptid, fup_month, weight, SBP, DBP)) %>%
  arrange(ptid, fup_month)




data_fup <- local({
    comorb_enrolment <-
      readxl::excel_sheets(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx")) %>%
      `[`(str_detect(., "comorbidity") & str_detect(., "^0M")) %>%
      set_names(str_extract(., "^[:digit:]*")) %>%
      map_dfr(
        ~ readxl::read_excel(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx"), sheet = .x, guess_max = 200000) %>%
          set_names(c("ptid", "diabetes", "has_CFF", "had_MI", "had_IHD", "has_PVD", "had_CVA", "has_COPD", "has_LD", "has_ST")),
        .id = "fup_month"
      ) %>%
      mutate(fup_month = as.numeric(fup_month)) %>%
      mutate_at(vars(has_CFF:has_ST), ~ if_else(. == 2, 0, .)) %>% # this is needed because these variables are coded differently at enrolment than at follow-up. what an absolute shitshow this database is.
      filter(!is.na(ptid))

    comorb_followup <-
      readxl::excel_sheets(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx")) %>%
      `[`(str_detect(., "comorbidity") & str_detect(., "^0M", negate = TRUE)) %>%
      set_names(str_extract(., "^[:digit:]*")) %>%
      map_dfr(
        ~ readxl::read_excel(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx"), sheet = .x, guess_max = 200000) %>%
          set_names(c("ptid", "diabetes", "has_CFF", "had_MI", "had_IHD", "has_PVD", "had_CVA", "has_COPD", "has_LD", "has_ST")),
        .id = "fup_month"
      ) %>%
      mutate(fup_month = as.numeric(fup_month)) %>%
      filter(!is.na(ptid))

    comorb <- bind_rows(comorb_enrolment, comorb_followup)

    fup <-
      readxl::excel_sheets(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx")) %>%
      `[`(str_detect(., "FU")) %>%
      set_names(str_extract(., "^[:digit:]*")) %>%
      map_dfr(
        ~ readxl::read_excel(path = here::here("data-raw", "mike", "Baseline and FU Data.xlsx"), sheet = .x, guess_max = 200000) %>%
          set_names(c("ptid", "date")),
        .id = "fup_month"
      ) %>%
      mutate(
        fup_month = as.numeric(fup_month),
        date = ymd(date)
      ) %>%
      bind_rows(data_baseline %>% select(ptid, fup_month, date), .) %>%
      filter(!is.na(ptid))


    full_join(fup, comorb, by = c("ptid", "fup_month") ) %>%
    left_join(
      data_physical,
      by = c("ptid", "fup_month")
    ) %>%
    select(ptid, fup_month, date, everything()) %>%
    arrange(ptid, fup_month, date) %>%
    filter(!is.na(date)) %>%
    mutate(
      had_MI = if_else(fup_month == 144 & had_MI == 2, 0, had_MI), # because MI is incorrectly coded for this follow up period
      has_ST = if_else(fup_month == 24 & has_ST == 4, 1, has_ST), # because ST is incorrectly coded for this follow up period

      # correct some dodgy dates
      date = case_when(
        ptid == 877 & fup_month == 60 & date == ymd("2012-05-09") ~ ymd("2013-05-09"),
        ptid == 325 & fup_month == 84 & date == ymd("2011-03-11") ~ ymd("2012-03-11"),
        ptid == 583 & fup_month == 36 & date == ymd("2007-02-19") ~ ymd(NA), # remove
        ptid == 604 & fup_month == 36 & date == ymd("2009-03-13") ~ ymd(NA), # remove this obs as post-death
        ptid == 1108 & fup_month == 60 & date == ymd("2007-10-16") ~ ymd("2008-10-16"),

        ptid == 1303 & fup_month == 12 & date == ymd("2003-06-25") ~ ymd("2004-06-25"),
        ptid == 1565 & fup_month == 12 & date == ymd("2008-10-05") ~ ymd("2009-10-05"),

        ptid == 1419 & fup_month == 108 & date == ymd("2011-10-03") ~ ymd("2012-10-03"),
        ptid == 2239 & fup_month == 84 & date == ymd("2014-04-23") ~ ymd(NA), # remove this obs as post-death
        TRUE ~ date
      ),
    ) %>%
    filter(!is.na(date))
})



# import hospital admissions data -----------------------------------------



data_CVadmissions <- read_csv(file = here::here("data-raw", "james", "data_CVevents.csv")) %>%
  transmute(
    ptid,
    fup_month = fup,
    index,
    date = admission_date,
    date_admission = admission_date,
    date_discharge = discharge_date,
    event_code,
    event = case_when(
      event_code == 1 ~ "CCF",
      event_code == 2 ~ "MI",
      event_code == 3 ~ "PVD",
      event_code == 4 ~ "CVA",
      event_code == 5 ~ "CVP",
    )
  ) %>%
  arrange(ptid, fup_month, date) %>%
  filter(ptid %in% data_baseline$ptid)




# import outcomes data (RRT, tx and death) --------------------------------



data_outcomes <-
  readxl::read_excel(path = here::here("data-raw", "mike", "Outcome RRT&RIP.xlsx"), sheet = 1, guess_max = 20000) %>%
  transmute(
    ptid = `CRISIS Study ID`,
    date_enrol = date(`Study start date`),
    date_death = date(`Date of Death`),
    date_RRTstart = date(`RRT Start Date`),
    RRT_modality = `RRT Start Type`,
    date_tx = date(`TX Date`),
    lastop_date = date(`Last OP Date`),
    lastdischarge_date = date(`Last Discharge Renal Date`)
  ) %>%
  left_join(
    data_baseline %>% select(ptid, date_birth),
    by = "ptid"
  ) %>%
  mutate(
    date_tx = if_else(date_tx <= date_enrol, ymd(NA), date_tx)
  )


# calculate all non-clinic events (birth, admissions, dialysis, tx --------

set.seed(420)
data_events <- local({

    CVadmissions <- data_CVadmissions %>%
      select(ptid, fup_month, date, event) %>%
      mutate(event = str_c("CV_", event))

    death <- data_outcomes %>%
      filter(!is.na(date_death)) %>%
      transmute(
        ptid,
        date = date_death,
        event = "death",
        fup_month = 9992
      )

    RRT <- data_outcomes %>%
      filter(!is.na(RRT_modality)) %>%
      transmute(
        ptid,
        date = date_RRTstart,
        event = str_c("RRT_", RRT_modality),
        fup_month = 9991
      )

    TX <- data_outcomes %>%
      filter(!is.na(date_tx)) %>%
      transmute(
        ptid,
        date = date_tx,
        event = "RRT_TX",
        fup_month = 9990
      )

    FUP <- data_fup %>%
      transmute(
        ptid,
        fup_month,
        date = date,
        event = "clinic"
      )


    bind_rows(
      CVadmissions, death, RRT, TX, FUP
    ) %>%
      # this chunk imputes hospitalisation date as somewhere after the previous follow up and before the current follow up, and repeats where there is more than one episode
      arrange(ptid, fup_month, date) %>%
      group_by(ptid) %>%
      mutate(
        date_new1 = if_else(
          is.na(date) & str_detect(event, "CV"),
          lag(date, 2) + (lag(date, 1) - lag(date, 2)) * runif(1, min = 0.05, max = 0.95),
          date
        )
      ) %>%
      arrange(ptid, fup_month, date_new1) %>%
      group_by(ptid) %>%
      mutate(
        date_new2 = if_else(
          is.na(date_new1) & str_detect(event, "CV"),
          lag(date_new1, 2) + (lag(date_new1, 1) - lag(date_new1, 2)) * runif(1, min = 0.05, max = 0.95),
          date_new1
        )
      ) %>%
      arrange(ptid, fup_month, date_new2) %>%
      group_by(ptid) %>%
      mutate(
        date_new3 = if_else(
          is.na(date_new2) & str_detect(event, "CV"),
          lag(date_new2, 2) + (lag(date_new2, 1) - lag(date_new2, 2)) * runif(1, min = 0.05, max = 0.95),
          date_new2
        )
      ) %>%
      ungroup() %>%
      arrange(ptid, fup_month, date_new3) %>%
      filter(event != "clinic") %>%
      mutate(date = date_new3) %>%
      select(-date_new1, -date_new2, -date_new3) %>%
      filter(ptid %in% data_baseline$ptid)

})


# extract creatinine data -------------------------------------------------------



data_creatinine <-
  bind_rows(
    readxl::read_xlsx(path = here::here("data-raw","mike","All creatinine upper half1.xlsx"), sheet=1),
    readxl::read_xlsx(path = here::here("data-raw","mike","All creatinine lower half1.xlsx"), sheet=1)
  ) %>%
  transmute(
    ptid=crisis_id,
    location=Location,
    date=date(TestDate),
    creatinine=Value
  ) %>%
  arrange(ptid, date) %>%
  filter(
    ptid %in% data_baseline$ptid,
    !is.na(creatinine)
    )

# this fuzzy matches dates of creatinine values (taken from EHR) to dates of SBP, DBP and weight (taken specifically for SKS)
# and keeps only those within a month of SKS clinic date
# also assumes date in SKS is true date (which is likely incorrect but makes later merge easier and is only off by max 1 month)
data_creatinine_cliniconly <-
  data_fup %>% select(ptid, date) %>% distinct() %>%
  fuzzyjoin::fuzzy_left_join(
    data_creatinine,
    by=c("ptid"="ptid", "date"="date"),
    match_fun=list(ptid=`==`, date=distday)#,
    #maxday=7
  ) %>%
  mutate(
    datediff = abs(lubridate::interval(date.x, date.y)/days(1))
  ) %>%
  arrange(ptid.x, date.x, datediff) %>%
  distinct(ptid.x, date.x, .keep_all=TRUE) %>%
  select(ptid=ptid.x, date=date.x, date_EHR=date.y, creatinine, location) %>%
  filter(!is.na(creatinine))

# similar to above but takes creatinine values for non-follow-up event
data_creatinine_eventsonly <-
  data_events %>% select(ptid, date) %>% distinct() %>%
  fuzzyjoin::fuzzy_left_join(
    data_creatinine,
    by=c("ptid"="ptid", "date"="date"),
    match_fun=list(ptid=`==`, date=distday)#,
    #maxday=7
  ) %>%
  mutate(
    datediff = abs(lubridate::interval(date.x, date.y)/days(1))
  ) %>%
  arrange(ptid.x, date.x, datediff) %>%
  distinct(ptid.x, date.x, .keep_all=TRUE) %>%
  select(ptid=ptid.x, date=date.x, date_EHR=date.y, creatinine, location) %>%
  filter(!is.na(creatinine))

data_creatinine_nrevents <-
  bind_rows(data_creatinine_cliniconly, data_creatinine_eventsonly) %>%
  arrange(ptid, date) %>%
  distinct()

write_rds(data_creatinine_nrevents, path=here::here("data-processed", "creatinine_nrevents.rds"))
#data_creatinine_nrevents=read_rds(here::here("data-processed", "creatinine_nrevents.rds"))


data_longitudinal <-
  data_fup %>%
  mutate(event = "clinic") %>%
  bind_rows(
    data_events %>%
      filter(event %ni% c("birth")) %>%
      select(ptid, event, date)#
  ) %>%
  left_join(data_creatinine_nrevents, by = c("ptid", "date")) %>%
  right_join(
    data_baseline %>%
      select(ptid, date_birth, date_enrol = date, sex, ethnicity, PRD, height),
    by = "ptid"
  ) %>%
  arrange(ptid, date) %>%
  group_by(ptid) %>%
  mutate(
    is_enrolment = (fup_month %in% 0),
    event_num = case_when(
      event %in% "death" ~ 4L,
      str_detect(event, "^RRT") ~ 3L,
      str_detect(event, "^CV") ~ 2L,
      event %in% "clinic" ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  arrange(ptid, date, event_num) %>%
  group_by(ptid) %>%
  mutate(
    has_RRT = cumany(str_detect(event, "^RRT") * 1),
    has_CV = cumany(str_detect(event, "^CV") * 1),
    has_died = cumany(str_detect(event, "death") * 1),

   lag_date = lag.skipna(if_else(creatinine == creatinine, date, ymd(NA))),
    date_diff_years = lubridate::interval(lag_date, date) / years(1),

    lag_creatinine = lag.skipna(creatinine),
    d_creatinine = creatinine - lag_creatinine,
    r_creatinine = d_creatinine / date_diff_years,

    lncreat = log(creatinine),
    lag_lncreat = lag(lncreat),
    d_lncreat = lncreat - lag_lncreat,
    r_lncreat = d_lncreat / date_diff_years,
  ) %>%
  ungroup() %>%
  mutate(

    age = lubridate::interval(start = date_birth, end = date) / years(1),
    fup_years = lubridate::interval(start = date_enrol, end = date) / years(1),
    fup_days = lubridate::interval(start = date_enrol, end = date) / days(1),
    BMI = weight / (height^2),
    MAP = (SBP + 2 * DBP) / 3,
    RAP = SBP - DBP,
    eGFR = 175 * ((creatinine / 88.4)^(-1.154)) * (age^(-0.203)) * (0.742^(sex == "female")) * (1.21^(ethnicity == "black")),
    is_egfrlt10 = eGFR < 10
  ) %>%
  filter( # remove events recorded after RRT or death
     !is.na(date)
  ) %>%
  filter(
    #if RRT or death occur on same day as follow-up, then assume no accual follow up occurred (unless SBP or DBP are measured) and remove fup row
    !(ptid == 2 & date == date("2010-04-27") & fup_month %in% 72 & event=="clinic"),
    !(ptid == 71 & date == date("2005-06-02") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 159 & date == date("2004-08-24") & fup_month %in% NA & event=="CV_CVP"),
    !(ptid == 194 & date == date("2011-05-10") & fup_month %in% NA & event=="death"),
    !(ptid == 256 & date == date("2009-03-16") & fup_month %in% 60 & event=="clinic"),
    !(ptid == 315 & date == date("2008-02-27") & fup_month %in% NA & event=="CV_CVP"),
    !(ptid == 908 & date == date("2003-03-14") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 920 & date == date("2003-05-06") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 923 & date == date("2003-11-18") & fup_month %in% NA & event=="death"),
    !(ptid == 940 & date == date("2005-08-05") & fup_month %in% 36 & event=="clinic"),
    !(ptid == 941 & date == date("2003-02-07") & fup_month %in% NA & event=="death"),
    !(ptid == 944 & date == date("2003-03-06") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 948 & date == date("2003-01-22") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 953 & date == date("2003-06-18") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 969 & date == date("2003-08-13") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 980 & date == date("2008-10-10") & fup_month %in% 60 & event=="clinic"),
    !(ptid == 996 & date == date("2003-07-28") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1002 & date == date("2003-09-10") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1006 & date == date("2003-09-23") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1020 & date == date("2003-11-23") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1021 & date == date("2003-03-09") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1035 & date == date("2003-08-06") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1043 & date == date("2003-07-21") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1051 & date == date("2003-09-23") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1065 & date == date("2003-07-07") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1098 & date == date("2006-04-19") & fup_month %in% 36 & event=="clinic"),
    !(ptid == 1139 & date == date("2003-04-07") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1140 & date == date("2003-05-13") & fup_month %in% NA & event=="death"),
    !(ptid == 1167 & date == date("2003-09-02") & fup_month %in% NA & event=="death"),
    !(ptid == 1213 & date == date("2003-07-23") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1225 & date == date("2003-12-25") & fup_month %in% NA & event=="death"),
    !(ptid == 1277 & date == date("2007-12-06") & fup_month %in% 60 & event=="clinic"),
    !(ptid == 1375 & date == date("2004-03-18") & fup_month %in% NA & event=="death"),
    !(ptid == 1387 & date == date("2004-04-10") & fup_month %in% 12 & event=="clinic"),
    !(ptid == 1425 & date == date("2005-10-11") & fup_month %in% 36 & event=="clinic"),
    !(ptid == 1906 & date == date("2012-06-19") & fup_month %in% 36 & event=="clinic"),
  )

# check no events occur on the same date and remove if necessary
#data_longitudinal %>% group_by(ptid) %>% filter(lag(date)==date | lead(date) == date) %>% view()

write_rds(data_longitudinal, path=here::here("data-processed", "data_longitudinal.rds"))


# canonical form of (anonymised, processed) dataset that will be used through the rest of the project
data_canon <- data_longitudinal %>%
  filter(
    fup_years>=0,
    # remove events recorded after first RRT or death
    !(has_died & event %ni% "death"), # remove anytihng occurring after death
    !(has_RRT & !str_detect(event,"^RRT")), # remove any non-RRTs occurring after RRT
    !(has_RRT & lag(has_RRT)), # remove RRTs occurring after RRT
  ) %>%
  # tsibble::as_tsibble(
  #   key = ptid,
  #   index = date,
  #   regular=FALSE
  # ) %>%
  select(
    #select all variables that will be used in analysis - add or remove here if necessary
    ptid,
    fup_month,
    fup_years,
    fup_days,
    event,
    is_enrolment,
    sex,
    ethnicity,
    PRD,
    age,
    diabetes:has_ST,
    height,weight,BMI,
    SBP,DBP,MAP, RAP,
    creatinine,lncreat,eGFR
  )

write_rds(data_canon, path=here::here("data", "data_canon.rds"))



################
