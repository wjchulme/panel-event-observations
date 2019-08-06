library('frailtypack')

data_ts <- readRDS(file = here::here("data", "data_longitudinal.rds")) %>%
  arrange(ptid, date) %>%
  group_by(ptid) %>%
  filter(!(date==first(date) & (str_detect(event,"^RRT")|str_detect(event, "death")))) %>%
  ungroup() %>%
  add_count(ptid, name="n_obs_pt") %>%
  filter(n_obs_pt!=1)




frailtyPenal(Surv(time, event) ~ cluster(ptid) + age + sex + lncreat, Frailty = TRUE, n.knots = 10, kappa1 = 1,
            data = data_ts, cross.validation = TRUE)