
data_canon <- read_rds(path = here::here("data", "data_canon.rds")) %>%
  filter(!is.na(age), !is.na(sex))

library('GGally')
library('cowplot')
# explore baseline variables -----------------------------------------------------------------

#age

plot_age <-
  data_canon %>%
  filter(is_enrolment) %>%
  ggplot() +
    #geom_blank() +
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=10)+
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=20)+
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=30)+
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=40)+
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=50)+
    geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=60)+
    geom_density(aes(age, stat(density), group=sex, colour=sex), adjust=0.85)+
    #geom_histogram(aes(age, stat(density), group=sex, fill=sex), position='identity', alpha=1/7, bins=70)+
    facet_grid(rows=vars(sex))+
    labs(x="Age at SKS recruitment, years")+
    theme_bw()+
    theme(legend.position='none')+
    NULL
plot_age
ggsave(here::here("figures", "plot_age.png"), width=150, height=150, units="mm")


#ethnicity
data_canon %>%
  filter(is_enrolment) %>%
  ggplot() +
    geom_bar(aes(ethnicity, fill=sex), position='dodge', size=0.4)+
    theme_bw()







# pairs -------------------------------------------------------------------



data_canon %>%
  as_tibble() %>%
  filter(is_enrolment) %>%
  select(age, sex, SBP, DBP, MAP, height, weight, BMI, creatinine, lncreat) %>%
  GGally::ggscatmat(columns = c('age', 'SBP', 'DBP', 'MAP', 'creatinine', 'lncreat', 'height', 'weight', 'BMI'), color='sex', alpha=0.2)


plot_pairs <-
data_canon %>%
  as_tibble() %>%
  filter(is_enrolment) %>%
  select(sex, age, SBP, DBP, MAP, height, weight, BMI, creatinine, lncreat) %>%
  GGally::ggpairs(
    columns = names(.)[-1],
    #columns = c('age', 'SBP', 'DBP', 'MAP', 'creatinine', 'lncreat', 'height', 'weight', 'BMI'),
    mapping=aes(color=sex, alpha=0.2),

    upper = list(
      continuous = GGally::wrap("points", alpha=0.2, size=0.1),
      combo = GGally::wrap("facethist")
    ),
    diag = list(
      continuous = GGally::wrap("densityDiag", alpha=0.3, colour="transparent"),
      discrete = GGally::wrap("barDiag")
    ),
    lower = list(
      continuous = GGally::wrap("points", alpha=0.2, size=0.1),
      combo = GGally::wrap("dot", alpha=0.2, size=0.1)
    ),
    title = "Correlation of measurements at basline"
  ) +
  theme_bw()


ggsave(here::here("figures", "plot_pairs.png"), width=150, height=150, units="mm")




# height versus weight ----------------------------------------------------


plot_height_weight <-
data_canon %>%
  filter(is_enrolment) %>%
  ggplot() +
  naniar::geom_miss_point(aes(height, weight, colour=sex), size=1, alpha=0.3)+
  labs(x="Height, m", y="Weight, kg")+
  theme_bw()+
  theme(legend.justification=c(0,1), legend.position = c(0.01, 0.99), legend.title = element_blank())

ggsave(here::here("figures","plot_height_weight.png"), plot_height_weight, width=150, height=150, units="mm")

# SBP versus DBP ----------------------------------------------------------


plot_DBP_SBP <-
data_canon %>%
  filter(is_enrolment) %>%
  ggplot() +
  naniar::geom_miss_point(aes(DBP, SBP, colour=sex), size=1, alpha=0.3)+
  labs(x="DBP, mmHg", y="SBP, mmHg")+
  theme_bw()+
  theme(legend.justification=c(1,0), legend.position = c(0.99, 0.01), legend.title = element_blank())

ggsave(here::here("figures","plot_BP.png"), plot_DBP_SBP, width=150, height=150, units="mm")



# BP / weight history -----------------------------------------------------------

plot_DBPagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, DBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="DBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_SBPagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, SBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="SBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_MAPagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, MAP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="MAP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_diffBPagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, SBP-DBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="SBP-DBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_creatinineagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, creatinine, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    scale_y_continuous(trans='log', breaks=c(50, 100, 200, 500, 1000, 3000))+
    labs(x=NULL, y=expression("Creatinine, "*mu*"mol"/"L"))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_weightagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, weight, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x="Patient age, years", y="Weight, kg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_BMIagetraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(age, BMI, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x="Patient age, years", y=expression("BMI, "*kg/m^2))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_DBPfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, DBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="DBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_SBPfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, SBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="SBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_MAPfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, MAP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="MAP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_diffBPfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, SBP-DBP, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x=NULL, y="SBP-DBP, mmHg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_creatininefuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, creatinine, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    scale_y_continuous(trans='log', breaks=c(50, 100, 200, 500, 1000, 3000))+
    labs(x=NULL, y=expression("Creatinine, "*mu*"mol"/"L"))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL


plot_weightfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, weight, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x="Follow-up, years", y="Weight, kg")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL

plot_BMIfuptraj <-
  data_canon %>%
  ggplot() +
    geom_line(aes(fup_years, BMI, group=ptid, colour=as.character(ptid)), size=0.5, alpha=0.2)+
    scale_colour_discrete(guide='none')+
    labs(x="Follow-up, years", y=expression("BMI, "*kg/m^2))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border=element_blank(),
          axis.line = element_line(colour = "black")
          ) +
    NULL



plot_BPtrajectories <-
  cowplot::plot_grid(
    plot_SBPagetraj, plot_SBPfuptraj,
    plot_DBPagetraj, plot_DBPfuptraj,
    plot_diffBPagetraj, plot_diffBPfuptraj,
    plot_MAPagetraj, plot_MAPfuptraj,
    ncol=2, align="hv"
  )

ggsave(here::here("figures", "plot_BPtrajectories.png"), plot_BPtrajectories, width=150, height=200, units="mm")


plot_trajectories <-
  cowplot::plot_grid(
    plot_creatinineagetraj, plot_creatininefuptraj,
    plot_weightagetraj, plot_weightfuptraj,
    plot_BMIagetraj, plot_BMIfuptraj,
    ncol=2, align="hv"
  )

ggsave(here::here("figures", "plot_trajectories.png"), plot_trajectories, width=150, height=200, units="mm")

# event history -----------------------------------------------------------

plot_eventhist <-
data_canon %>%
  group_by(ptid) %>%
  filter(event!="birth",  (event!="clinic" | is_enrolment | row_number()==n()), n()>1) %>%
  mutate(
    maxfup=max(fup_years),
    event=case_when(
      str_detect(event, "^RRT_TX") ~ "transplant",
      str_detect(event, "^RRT") ~ "dialysis",
      str_detect(event, "death") ~ "death",
      str_detect(event, "clinic") & is_enrolment ~ "enrolment",
      str_detect(event, "clinic") & !is_enrolment ~ "last follow-up",
      str_detect(event, "^CV") ~ "CV hospitalisation",
      TRUE ~ NA_character_
    ) %>% factor(levels=c("enrolment", "CV hospitalisation", "last follow-up", "transplant", "dialysis", "death"))
  ) %>%
  ungroup() %>%
  mutate(
    #maxfupbins = ntile(maxfup+ptid/1000000, 4),
    binid = dense_rank(maxfup+ptid/1000000),
    binid4 = binid%%((max(binid)+1)/4),
    grpid4 = floor(4*binid/(max(binid)+1))
  ) %>%
  arrange(binid) %>%
  # group_by(maxfupbins) %>%
  # mutate(binid=dense_rank(maxfup+ptid/1000000)) %>%
  # ungroup() %>%
  ggplot() +
    geom_line(aes(fup_years, -binid4, group=ptid), size=0.01, alpha=0.1)+
    geom_point(aes(fup_years, -binid4, colour=event
                   #,shape=case_when(is_egfrlt10 %in% TRUE ~ "eGFR<10", TRUE~"eGFR>10")
                   ), size=0.4)+
    scale_color_manual(
      values=c("enrolment"="lightgray", "CV hospitalisation"="#D55E00", "last follow-up"="transparent", "transplant"="#009E73", "dialysis"="#0072B2", "death"="#000000"),
      guide = guide_legend(override.aes = list(size=2))
    )+
    scale_shape_manual(values=c(4, 20),  guide = guide_legend(override.aes = list(size=2)))+
    scale_x_continuous(breaks=seq(0,20,1))+
    #facet_wrap(~item, ncol=1)+
    #scale_colour_discrete(guide='none')+
    labs(x='Follow-up year', y=NULL, shape=NULL)+
    facet_grid(cols=vars(grpid4), scales='free', space="free")+
    theme_bw() +
    #guides(colour = guide_legend(override.aes = list(size=2)))+
    #guide_legend(size=2)+
    theme(
      legend.position="bottom",
      #legend.key.size = unit(3,"line"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_blank(), strip.background = element_blank(),
      panel.grid = element_blank(), panel.border = element_blank()
      ) +
    NULL


ggsave(here::here("figures", "plot_eventhist.png"), plot_eventhist, width=150, height=220, units="mm", dpi=600)




# check CV history and event rates ----------------------------------------------------

# compare CV history with data from this paper:
# https://onlinelibrary.wiley.com/doi/full/10.1046/j.1525-139X.2003.16025.x?sid=nlm%3Apubmed


data_canon <- read_rds(file = here::here("data", "data_canon.rds")) %>% as_tibble()

write_rds(data_canon, path = here::here("data", "data_canon.rds"))

  arrange(ptid, date)

data_enrol <- data_canon %>%
  group_by(ptid) %>%
  filter(event=='clinic') %>%
  filter(min(fup_month)==fup_month)


data_preRRT <- data_canon %>%
  group_by(ptid) %>%
  filter(event=="clinic") %>%
  filter(max(fup_month)==fup_month)


"At least 35% of patients with CKD have evidence of an ischemic event (myocardial infarction or angina) at the time of presentation to a nephrologist."
table(data_enrol$had_MI) %>% prop.table()
table(data_preRRT$had_MI) %>% prop.table()


data_canon %>%
  group_by(ptid) %>%
  summarise(
    MI = any(had_MI>0 | str_detect(event, "^RRT"))
  ) %>%
  ungroup() %>%
  with(table(MI)) %>%
  prop.table()

# this suggests significant under-reporting of previous history of MI



# compare CV event rates with data from this paper:
# https://bjgp.org/content/68/673/e512/tab-figures-data

data_canon %>%
  as_tibble() %>%
  group_by(ptid) %>%
  summarise(
    #each follow-up clinic covers approximately one year so assume one year per event=="clinic" row
    rateMI = sum(str_detect(event, "CV_MI"))/sum(event=="clinic"),
    rateCVA = sum(str_detect(event, "CV_CVA"))/sum(event=="clinic"),
    rateCCF = sum(str_detect(event, "CV_CCF"))/sum(event=="clinic")
  ) %>%
  ungroup() %>%
  summarise(
    avgrateMI = mean(rateMI)*1000,
    avgrateCVA = mean(rateCVA)*1000,
    avgrateCCF = mean(rateCCF)*1000,
    )

# slightly underestimates incidence per 1000 person-years
# but not bad


