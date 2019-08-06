
library('msm')
library('minqa') #use if using opt.method = "bobyqa" in msm


data_all <- read_rds(path = here::here("data", "data_canon.rds")) %>%
  filter(!is.na(age), !is.na(sex), !is.na(height)) %>%
  arrange(ptid, date) %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)


data_all %>% group_by(ptid) %>% summarise(firstevent=first(event)) %>% with(table(firstevent))
table(data_all$event)
msm::statetable.msm(event, ptid, data=data_all)
print(msm::statetable.msm(event, ptid, data=data_all))
with(data_all, table(ptid, event))

## Need to think about censoring the observation period before deaths are observed, since panel-observation process may have stopped long before death occurs.
## Need to also think about hospitalisations occurring after the last observation date - they won't be recorded as they're self-reported at the following out-patient apmt
## Death information is taken from ?? in SKS

## Fudge: censor all death/RRT events occurring beyond 12 months after the last observation date.

data_censored <- data_all %>%
  group_by(ptid) %>%
  mutate(
    censored = (event!="clinic" & !str_detect(event,"^CV") & row_number(ptid)==length(ptid) & date > lag(date)+365.25),
    date = if_else(censored, pmin(date, lag(date)+365.25), date), # this censors any absorbing states so that it's
    event = if_else(censored, "censored", event)
  ) %>%
  ungroup() %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1) %>% select(-n_obs_pt)


# MODEL 1 -----------------------------------------------------------------



# multi-state model (no hidden states) --------

data_msm1 <- data_censored %>%
  group_by(ptid) %>%
  fill(weight, BMI, SBP, DBP, lncreat) %>% #this is done to stop msm removing rows where these values are missing - LOCF is used for panel states, and values are not used for event-states (see value passed to covariates argument)
  ungroup() %>%
  mutate(
    months = fup_years/12,

    scale.age=scale(age),
    scale.lncreat=scale(lncreat),

    state =
      case_when(
        event %in% "clinic" ~ 1,
        str_detect(event, "^CV") ~ 2,
        str_detect(event, "^RRT") ~ 3,
        event=="death" ~ 4,
        event=="censored" ~ 99
      ),
    obstype =
      case_when(
        event %in% "clinic" ~ 1,
        str_detect(event, "^CV") ~ 3,
        str_detect(event, "^RRT") ~ 3,
        event %in% "death" ~ 3,
        event=="censored" ~ 1
      ),
    obstrue= 1,


    lncreat =
      case_when(
        state == 1 ~ 0,
        TRUE ~ lncreat
      ),

    SBP = case_when(
        state == 1 ~ 0,
        TRUE ~ SBP
      ),

    DBP = case_when(
        state == 1 ~ 0,
        TRUE ~ DBP
      ),

    BMI = case_when(
        state == 1 ~ 0,
        TRUE ~ BMI
      ),

  ) %>%
  filter(!is.na(state)) %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)




Pm1 <- rbind(
  c(0.9, 0.05, 0.025, 0.025),
  c(0.85, 0.05,  0.05,  0.05),
  c(0, 0, 1, 0),
  c(0, 0, 0, 1)
)
rowSums(Pm1)
Qm1 <- expm::logm(Pm1, method="Eigen")*2
rowSums(Qm1)

covariates.in1 <- list(
  "1-2" = ~ age + sex + BMI + lncreat + SBP,
  "1-3" = ~ age + sex + BMI + lncreat + SBP,
  "1-4" = ~ age + sex + BMI + lncreat + SBP,
  "2-1" = ~ age + sex,
  "2-3" = ~ age + sex,
  "2-4" = ~ age + sex
)

constraint.in1 <- list(
    age = c(1,2,3,4,2,3),
    sexfemale = c(1,2,3,4,2,3)
    #BMI=c(1,2,3,4,2,3),
    #lncreat=c(1,2,3,4,2,3),
    #SBP=c(1,2,3,4,2,3)
  )

hmodel.in1 = list(
  hmmIdent(1),
  hmmIdent(2),
  hmmIdent(3),
  hmmIdent(4)
)


# standard MSM with exact observed trnasitions to hosp, RRT, and death
msmFit1a <- msm::msm(
  data = data_msm1,
  formula = state ~ months,
  subject = ptid,
  qmatrix = Qm1,
  initprobs = c(1, 0, 0, 0),
  gen.inits = TRUE,
  #fixedpars = TRUE,
  censor.states = 99, censor.states = c(1)
  obstype = obstype,
  covariates = covariates.in1,
  constraint = constraint.in1,
  opt.method = "bobyqa" # default = "optim"

  #control=list(fnscale=20000, maxit=23000, reltol = 1e-10)
  #method="Nelder-Mead" #default="BFGS"
  )

write_rds(msmFit1a, path=here::here("data","models", "censored 1 year", "model 1a - age sex BMI lncreat SBP.rds"))

# dressed up as a HMM (but exactly the same as MSM above)
msmFit1b <- msm::msm(
  data = data_msm1,
  formula = state ~ months,
  subject = ptid,
  qmatrix = Qm1,
  initprobs = c(1, 0, 0, 0),
  gen.inits = TRUE,
  #fixedpars = TRUE,
  censor.states = 99, censor.states = c(1)
  obstype = obstype,
  covariates = covariates.in1,
  constraint = constraint.in1,
  hmodel = hmodel.in1,
  opt.method = "bobyqa" # default = "optim"

  #control=list(fnscale=20000, maxit=23000, reltol = 1e-10)
  #method="Nelder-Mead" #default="BFGS"
  )

write_rds(msmFit1a, path=here::here("data","models","censored 1 year", "model 1b - age sex BMI lncreat SBP.rds"))

summary(msmFit1a)
summary(msmFit1b)


# MODEL 2 -----------------------------------------------------------------



# misclassification model --------


data_msm2 %>% filter(event=="clinic" & !is.na(BMI) & !is.na(SBP) & !is.na(lncreat) & !is.na(obstype)) %$% ptid %>% unique() %>% `[`(821)

data_msm2 <- data_ts %>%
  group_by(ptid) %>%
  fill(weight, BMI, SBP, DBP, lncreat) %>% #this is done to stop msm removing rows where these values are missing LOCF is used for panel states, and values are not used for event-states (see value passed to covariates argument)
  ungroup() %>%
  mutate(
    months = fup_years/12,

    scale.age=scale(age),
    scale.lncreat=scale(lncreat),

    state =
      case_when(
        event %in% "clinic" & (SBP<=140 | is.na(SBP)) ~ 1,
        event %in% "clinic" & SBP>140 ~ 2,
        str_detect(event, "^CV") ~ 3,
        str_detect(event, "^RRT") ~ 4,
        event=="death" ~ 5,
        event=="censored" ~ 99,
      ),
    obstype =
      case_when(
        event %in% "clinic" ~ 1,
        str_detect(event, "^CV") ~ 3,
        str_detect(event, "^RRT") ~ 3,
        event == "death" ~ 3,
        event == "censored" ~ 1,
      ),

  ) %>%
  filter(!is.na(state)) %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)


Pm2 <- rbind(
      c(0.7, 0.2, 0.05, 0.025, 0.025),
      c(0.1, 0.6, 0.2,  0.05,  0.05),
      c(0.1, 0.5, 0,  0.15,  0.25),
      c(0, 0, 0, 1, 0),
      c(0, 0, 0, 0, 1)
    )
rowSums(Pm2)

Qm2 <- expm::logm(Pm2, method="Eigen")
rowSums(Qm2)

Em2 <- rbind(
  c(0, 0.3, 0, 0, 0),
  c(0.3, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0),
  c(0, 0, 0, 0, 0)
)

hmodel.in2 = list(
  hmmCat(c(0.7,0.3,0,0,0)),
  hmmCat(c(0.3,0.7,0,0,0)),
  hmmIdent(3),
  hmmIdent(4),
  hmmIdent(5)
)


covariates.in2 = list(
  "1-2" = ~ age + sex + BMI + lncreat + SBP,
  "1-3" = ~ age + sex + BMI + lncreat + SBP,
  "1-4" = ~ age + sex + BMI + lncreat + SBP,
  "1-5" = ~ age + sex + BMI + lncreat + SBP,
  "2-1" = ~ age + sex + BMI + lncreat + SBP,
  "2-3" = ~ age + sex + BMI + lncreat + SBP,
  "2-4" = ~ age + sex + BMI + lncreat + SBP,
  "2-5" = ~ age + sex + BMI + lncreat + SBP,
  "3-1" = ~ age + sex,
  "3-2" = ~ age + sex,
  "3-4" = ~ age + sex,
  "3-5" = ~ age + sex
)


constraint.in2 <- list(
    age = c(1,1,2,3,4,1,2,3, 5,6,7,8),
    sexfemale = c(1,1,2,3,4,1,2,3, 5,6,7,8)
    #BMI=c(1,2,3,4,2,3),
    #lncreat=c(1,2,3,4,2,3),
    #SBP=c(1,2,3,4,2,3)
  )




msmFit2 <- msm::msm(
  data = data_msm2,
  formula = state ~ fup_years,
  subject = ptid,
  qmatrix =Qm2,
  ematrix = Em2,
  #hmodel = hmodel.in2,
  initprobs = c(0.6, 0.4, 0, 0, 0),
  est.initprobs = TRUE,
  #gen.inits = TRUE,
  #fixedpars = TRUE,
  #hcovariates = list(~age+sex, ~age+sex, NULL, NULL, NULL),
  covariates = covariates.in2,
  constraint = constraint.in2,
  obstype = obstype,
  censor = 99, censor.states=c(1, 2), # does this work for HHMs?

  #opt.method = "bobyqa" # default = "optim"

  #method="Nelder-Mead" #default="BFGS"
  control = list(fnscale=2000,
                  maxit=1000, #default=500 for Nelder mead
                  reltol=1e-10#, #default~=1e-08
                 # ndeps = rep(1e-05, n.pars) #default=1e-3
                 )
  )

#length((msmFit1$data$mf$`(subject)`))
#((msmFit2$data$mf$`(subject)`))

write_rds(msmFit2, path=here::here("data","models", "censored 1 year", "model 2 - age sex BMI lncreat SBP.rds"))







# log creatinine and SBP are emission distributions, with CV hospitalisation events included as a transient state  -----------------------------------------------------------------------
# RRT and death are absorbing states


data_msm3 <-
  data_ts %>%
  group_by(ptid) %>%
  fill(weight, BMI, SBP, DBP, lncreat) %>% #this is done to stop msm removing rows where these values are missing LOCF is used for panel states, and values are not used for event-states (see value passed to covariates argument)
  ungroup() %>%
  mutate(
    state =
      case_when(
        event %in% "clinic" ~ -1,
        str_detect(event, "^CV") ~ -3,
        str_detect(event, "^RRT") ~ -4,
        event=="death" ~ -5,
        event == "censored" ~ -99
      ),

    obstype =
      case_when(
        event %in% c("clinic", "censored") ~ 1,
        event %ni% "clinic" ~ 3
      ),

    lncreat =
      case_when(
        state == -1 ~ lncreat,
        TRUE ~ state
      ),

    SBP = case_when(
        state == -1 ~ SBP,
        TRUE ~ state
      ),

    DBP = case_when(
        state == -1 ~ DBP,
        TRUE ~ state
      ),

    BMI = case_when(
        state == -1 ~ BMI,
        TRUE ~ state
      ),


  ) %>%
  filter(
    !(is.na(SBP) & is.na(lncreat) & is.na(BMI) & state=='clinic')
  ) %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)
data_msm3$emission <- cbind(data_msm3$lncreat, data_msm3$SBP, data_msm3$BMI)

print(msm::statetable.msm(-state, ptid, data=data_msm3))



Pm3 <- rbind(
      c(0.7, 0.2, 0.05, 0.025, 0.025),
      c(0.1, 0.6, 0.2,  0.05,  0.05),
      c(0.1, 0.5, 0,  0.15,  0.25),
      c(0, 0, 0, 1, 0),
      c(0, 0, 0, 0, 1)
    )
rowSums(Pm3)

Qm3 <- expm::logm(Pm3, method="Eigen")
rowSums(Qm3)



hmodel.in3 <- list(

  # contorlled
  hmmMV(
    hmmNorm(6, 0.5),
    hmmNorm(150, 20),
    hmmNorm(25, 5)
  ),

  # uncontrolled
  hmmMV(
    hmmNorm(5, 0.5),
    hmmNorm(120, 15),
    hmmNorm(25, 5)
  ),

  # hospitalisation
  hmmMV(
    hmmIdent(-3),
    hmmIdent(-3),
    hmmIdent(-3)
  ),

  # RRT
  hmmMV(
    hmmIdent(-4),
    hmmIdent(-4),
    hmmIdent(-4)
  ),

  # death
  hmmMV(
    hmmIdent(-5),
    hmmIdent(-5),
    hmmIdent(-5)
  )
)


covariates.in3 = list(
  "1-2" = ~ age + sex + BMI + lncreat + SBP,
  "1-3" = ~ age + sex + BMI + lncreat + SBP,
  "1-4" = ~ age + sex + BMI + lncreat + SBP,
  "1-5" = ~ age + sex + BMI + lncreat + SBP,
  "2-1" = ~ age + sex + BMI + lncreat + SBP,
  "2-3" = ~ age + sex + BMI + lncreat + SBP,
  "2-4" = ~ age + sex + BMI + lncreat + SBP,
  "2-5" = ~ age + sex + BMI + lncreat + SBP,
  "3-1" = ~ age + sex,
  "3-2" = ~ age + sex,
  "3-4" = ~ age + sex,
  "3-5" = ~ age + sex
)


constraint.in3 <- list(
    age = c(1,1,2,3,4,1,2,3, 5,6,7,8),
    sexfemale = c(1,1,2,3,4,1,2,3, 5,6,7,8)
    #BMI=c(1,2,3,4,2,3),
    #lncreat=c(1,2,3,4,2,3),
    #SBP=c(1,2,3,4,2,3)
  )


msmFit3 <- msm::msm(
  data = data_msm3,
  formula = emission ~ fup_years,
  subject = ptid,
  qmatrix =Qm3,
  initprobs = c(0.8, 0.2, 0, 0, 0),
  est.initprobs = TRUE,
  #gen.inits = TRUE,
  #fixedpars = TRUE,
  hmodel = hmodel.in3,
  #hcovariates = list(~age+sex, ~age+sex, NULL, NULL, NULL),
  covariates = covariates.in3,
  constraint = constraint.in3,
  obstype = obstype,
  censor = -99, censor.states=c(1, 2), #this doesn't currently work for HMMs

  opt.method = "bobyqa" # default = "optim"

  #method="Nelder-Mead" #default="BFGS"
  # control = list(fnscale=2000,
  #                 maxit=1000, #default=500 for Nelder mead
  #                 reltol=1e-10#, #default~=1e-08
  #                # ndeps = rep(1e-05, n.pars) #default=1e-3
  #                )
  )


write_rds(msmFit3, path=here::here("data","models", "censored 1 year", "model 3 - age sex BMI lncreat SBP.rds"))







# log creatinine and SBP are emission distributions, with CV hospitalisation events included as a transient state  -----------------------------------------------------------------------
# RRT is a censoring state


Pm4 <- rbind(
      c(0.7, 0.2, 0.05, 0.05),
      c(0.1, 0.6, 0.2,  0.15),
      c(0.1, 0.5, 0,  0.4),
      c(0, 0, 0, 1)
    )
rowSums(Pm4)

Qm4 <- expm::logm(Pm4, method="Eigen")
rowSums(Qm4)

data_msm4 <-
  data_ts %>%
  mutate(
    state =
      case_when(
        event %in% "clinic" ~ -1,
        str_detect(event, "^CV") ~ -2,
        str_detect(event, "^RRT") ~ -99,
        event=="death" ~ -4
      ),

    obstype =
      case_when(
        event %in% "clinic" ~ 1,
        event %ni% "clinic" ~ 3
      ),

    lmcreat =
      case_when(
        state == -1 ~ lncreat,
        TRUE ~ state
      ),

    SBP = case_when(
        state == -1 ~ SBP,
        TRUE ~ state
      ),

    DBP = case_when(
        state == -1 ~ DBP,
        TRUE ~ state
      ),

  ) %>%
  filter(
    !(is.na(SBP) & is.na(lncreat) & state=='clinic')
  ) %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)

data_msm4$emission <- cbind(data_msm4$lncreat, data_msm4$SBP)

print(msm::statetable.msm(-state, ptid, data=data_msm4))

hmodel.in4 <- list(

  # contorlled
  hmmMV(
    hmmNorm(6, 0.5),
    hmmNorm(150, 20)
  ),

  # uncontrolled
  hmmMV(
    hmmNorm(5, 0.5),
    hmmNorm(120, 15)
  ),

  # hospitalisation
  hmmMV(
    hmmIdent(-2),
    hmmIdent(-2)
  ),

  # death
  hmmMV(
    hmmIdent(-4),
    hmmIdent(-4)
  )
)


msmFit4 <- msm::msm(
  data = data_msm4,
  formula = emission ~ fup_years,
  subject = ptid,
  qmatrix =Qm4,
  initprobs = c(0.9, 0.1, 0, 0),
  est.initprobs = TRUE,
  #gen.inits = TRUE,
  #fixedpars = TRUE,
  hmodel = hmodel.in4,
  #hcovariates = list(~age+sex, ~age+sex, NULL, NULL, NULL),
  covariates = ~ age + sex, #age + sex affects all transition intensities
  obstype = obstype,
  censor = -99, #does this work for hmms?
  #censor.states=c(1, 2),

  #opt.method = "bobyqa" # default = "optim"

  #method="Nelder-Mead" #default="BFGS"
  control = list(fnscale=2000,
                  maxit=1000, #default=500 for Nelder mead
                  reltol=1e-10#, #default~=1e-08
                 # ndeps = rep(1e-05, n.pars) #default=1e-3
                 )
  )
