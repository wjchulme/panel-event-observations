msmFit1 <- read_rds(here::here("data","models", "uncensored", "model 1b - age sex BMI lncreat SBP.rds"))
msmFit2 <- read_rds(here::here("data","models", "uncensored", "model 2 - age sex BMI lncreat SBP.rds"))
msmFit3 <- read_rds(here::here("data","models", "uncensored", "model 3 - age sex BMI lncreat SBP.rds"))


msmFit1 <- read_rds(here::here("data","models", "censored 1 year", "model 1b - age sex BMI lncreat SBP.rds"))
msmFit2 <- read_rds(here::here("data","models", "censored 1 year", "model 2 - age sex BMI lncreat SBP.rds"))
msmFit3 <- read_rds(here::here("data","models", "censored 1 year", "model 3 - age sex BMI lncreat SBP.rds"))


nrow(msmFit3$data$mf)

.msm.LOOKUP <- data.frame(
        label = msm:::.msm.HMODELS,
        hmmname = c("hmmCat", "hmmIdent", "hmmUnif", "hmmNorm", "hmmLNorm", "hmmExp", "hmmGamma", "hmmWeibull", "hmmPois", "hmmBinom", "hmmBetaBinom",
                    "hmmTNorm", "hmmMETNorm", "hmmMEUnif", "hmmNBinom", "hmmBeta", "hmmT",
                    "hmmClmTNorm", "hmmClmTNorm7"),
        distname = c("multinom", "", "unif", "norm", "lnorm", "exp", "gamma", "weibull", "pois", "binom", "",
                     "", "", "", "nbinom", "beta", "t",
                     "", ""),
        distrname = c("", "", "Unif", "Norm", "LNorm", "Exp", "", "Weibull", "Pois", "Binom", "",
                      "", "", "", "Nbinom", "Beta", "",
                      "",""),
        stringsAsFactors = FALSE
      )

# calculate passage probability for hospital state at every t_ij ----------



mixturepredCat.msm <- function(obj){
      #This is a mixture model since we only have one time-point (one observation)

  hmodel.in <- hmodel2list(obj$hmodel)
  mf <- obj$data$mf
  K <- length(hmodel.in)
#  stateprobMtx <- matrix(data=NA_real_, nrow=nrow(mf), ncol=K)

  stateprobMtx <- t(sapply(mf$`(state)`, FUN=function(state){
    hmodel.in[[state]]$pars[3:7]
  }))

}

mixtureprednormMV.msm <- function(obj){
  hmodel.in <- hmodel2list(obj$hmodel)
  mf <- obj$data$mf
  K <- length(hmodel.in)

  if(is.null(nrow(obj$hmodel$labels))) stop("needs multivariate input")

  D <- nrow(obj$hmodel$labels)

  stateprobMtx <- t(sapply(array_tree(mf$`(state)`,1), FUN=function(state){



    labels <- rep(NA_character_, K)
    densMV <- rep(NA_real_, K)

    for(k in 1:K){

      labels[k] <- hmodel.in[[k]][[1]]$label

      if(labels[k] == "identity"){

        densMV[k] <- (state[1]==hmodel.in[[k]][[1]]$pars)*1

      } else if(labels[k] == "normal"){
        pars <- matrix(NA_real_, nrow=D, ncol=2)
        dens <- rep(NA_real_, nrow=D)
        for(d in 1:D){
            pars[d,] <- hmodel.in[[k]][[d]]$pars
            dens[d] <- dnorm(state[d], pars[d,1], pars[d,2])
        }
        densMV[k] <- mean(dens)
      } else {stop("response distribution neither 'normal' nor 'identity'")}
    }


    if(sum(densMV*(labels=="identity"))>0){
      densMV[labels!="identity"] <- 0
    }

    densMV/sum(densMV)
  }))

  stateprobMtx
}

nrow(msmFit1$data$mf)
nrow(msmFit2$data$mf)
nrow(msmFit3$data$mf)




eval1 <-  msmFit1$data$mf %>%
  group_by(`(subject)`) %>%
  mutate(
    leadint = lead(`(time)`)-`(time)`,
    state=`(state)`,
    leadhosp = (lead(state)==2)*1,
  ) %>%
  ungroup() %>%
  mutate(
    row.num=row_number(),

    ppass1 = pmap_dbl(
              list(leadint, age, sex, BMI, lncreat, SBP),
              function(leadint, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(state!=1 | is.na(leadint)|is.na(BMI)|is.na(SBP)|is.na(lncreat)) {
                      NA
                    }else{
                      ppass.msm(msmFit1, tot=leadint, start=c(1,0,0,0), covariates=covlist)[2]
                    }
                }
            ),

  )


eval2 <-  msmFit2$data$mf %>%
  {.$stateprobFB2 = viterbi.msm(msmFit2)$pstate %>% array_tree(1); .} %>%
  {.$stateprobMtx2 = mixturepredCat.msm(msmFit2) %>% array_tree(1); .} %>%
  group_by(`(subject)`) %>%
  mutate(
    leadint = lead(`(time)`)-`(time)`,
    state=`(state)`,
    leadhosp = (lead(state==3)*1)

  ) %>%
  ungroup() %>%
  mutate(
    row.num=row_number(),


    ppass2 = pmap(
              list(leadint, age, sex, BMI, lncreat, SBP),
              function(leadint, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint)|is.na(BMI)|is.na(SBP)|is.na(lncreat)) {
                      NA
                    }else{
                      ppass.msm(msmFit2, tot=leadint, start="all", covariates=covlist)[c(1,2),3]
                    }
                }
            ),

    ppass2_1=map_dbl(ppass2, ~.[1]),
    ppass2_2=map_dbl(ppass2, ~.[2]),

    ppass2Mtx = pmap_dbl(
              list(leadint, state, stateprobMtx2, age, sex, BMI, lncreat, SBP),
              function(leadint, state, stateprobMtx2, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint) | is.na(BMI) | is.na(SBP) | is.na(lncreat) | state==3) {
                      NA
                    }else{
                      ppass.msm(msmFit2, tot=leadint, start=stateprobMtx2, covariates=covlist)[3]
                    }
                }
            ),

    ppass2FB = pmap_dbl(
              list(leadint, state, stateprobFB2, age, sex, BMI, lncreat, SBP),
              function(leadint, state, stateprobFB2, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint) | is.na(BMI) | is.na(SBP) | is.na(lncreat) | state==3) {
                      NA
                    }else{
                      ppass.msm(msmFit2, tot=leadint, start=stateprobFB2, covariates=covlist)[3]
                    }
                }
            ),


  )


eval3 <- msmFit3$data$mf %>%
  {.$stateprobFB3 = viterbi.msm(msmFit3)$pstate %>% array_tree(1); .} %>%
  {.$stateprobMtx3 = mixtureprednormMV.msm(msmFit3) %>% array_tree(1); .} %>%
  group_by(`(subject)`) %>%
  mutate(
    leadint = lead(`(time)`)-`(time)`,
    state = ifelse(`(state)`[,1]<0, -`(state)`[,1], 1),
    leadhosp = (lead(state==3)*1)

  ) %>%
  ungroup() %>%
  mutate(
    row.num=row_number(),


    ppass3 = pmap(
              list(leadint, age, sex, BMI, lncreat, SBP),
              function(leadint, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint)|is.na(BMI)|is.na(SBP)|is.na(lncreat)) {
                      NA
                    }else{
                      ppass.msm(msmFit3, tot=leadint, start="all", covariates=covlist)[c(1,2),3]
                    }
                }
            ),

    ppass3_1=map_dbl(ppass3, ~.[1]),
    ppass3_2=map_dbl(ppass3, ~.[2]),

    ppass3Mtx = pmap_dbl(
              list(leadint, state, stateprobMtx3, age, sex, BMI, lncreat, SBP),
              function(leadint, state, stateprobMtx3, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint) | is.na(BMI) | is.na(SBP) | is.na(lncreat) | state==3) {
                      NA
                    }else{
                      ppass.msm(msmFit3, tot=leadint, start=stateprobMtx3, covariates=covlist)[3]
                    }
                }
            ),

    ppass3FB = pmap_dbl(
              list(leadint, state, stateprobFB3, age, sex, BMI, lncreat, SBP),
              function(leadint, state, stateprobFB3, age, sex, BMI, lncreat, SBP){
                    covlist <- list(age, sex, BMI, lncreat, SBP)

                    if(is.na(leadint) | is.na(BMI) | is.na(SBP) | is.na(lncreat) | state==3) {
                      NA
                    }else{
                      ppass.msm(msmFit3, tot=leadint, start=stateprobFB3, covariates=covlist)[3]
                    }
                }
            ),


  )

plot.prevalence.msm(msmFit1)



prevalence.msm(msmFit1)
pearson.msm(msmFit2)


# pseudo time-to-event GoF ------------------------------------------------








# pseudo-binary GoF -------------------------------------------------------


ggplot(eval2) +
  geom_point(aes(x=leadint, y=ppass2Mtx, colour=as.factor(leadhosp)), alpha=0.2, size=1)+
  theme_bw()

ggplot(eval2) +
  geom_point(aes(x=leadint, y=ppass2FB, colour=as.factor(leadhosp)), alpha=0.2, size=1)+
  theme_bw()



library('pROC')
roc_obj1 <- roc(eval1$leadhosp, eval1$ppass1)

roc_obj2FB <- roc(eval2$leadhosp, eval2$ppass2FB)
roc_obj2Mtx <- roc(eval2$leadhosp, eval2$ppass2Mtx)
roc_obj2_1 <- roc(eval2$leadhosp, eval2$ppass2_1)
roc_obj2_2 <- roc(eval2$leadhosp, eval2$ppass2_2)


roc_obj3FB <- roc(eval3$leadhosp, eval3$ppass3FB)
roc_obj3Mtx <- roc(eval3$leadhosp, eval3$ppass3Mtx)
roc_obj3_1 <- roc(eval3$leadhosp, eval3$ppass3_1)
roc_obj3_2 <- roc(eval3$leadhosp, eval3$ppass3_2)


auc(roc_obj1)
auc(roc_obj2FB)
auc(roc_obj2Mtx)
auc(roc_obj3FB)
auc(roc_obj3Mtx)

plot(roc_obj1)
plot(roc_obj2FB)
plot(roc_obj2Mtx)
plot(roc_obj3FB)
plot(roc_obj3Mtx)


# lumpability -------------------------------------------------------------



# set of states cannot be lumped in general, because the resulting transition matrix may not be memoryless - see "lumpability"
Q1 <- qmatrix.msm(msmFit1)$estimates
Q2 <- qmatrix.msm(msmFit2)$estimates
Q3 <- qmatrix.msm(msmFit3)$estimates


sum(Q2[1,3:5])
sum(Q2[2,3:5])

sum(Q2[3:5,1])
sum(Q2[3:5,2])


Q2lumped <- rbind(
  c(sum(Q2[1:2,1:2]), colSums(Q2[1:2,3:5])),
  c(sum(Q2[3,1:2]), Q2[3,3:5]),
  c(sum(Q2[4,1:2]), Q2[4,3:5]),
  c(sum(Q2[5,1:2]), Q2[5,3:5])
)

Q3lumped <- rbind(
  c(sum(Q3[1:2,1:2]), colSums(Q3[1:2,3:5])),
  c(sum(Q3[3,1:2]), Q3[3,3:5]),
  c(sum(Q3[4,1:2]), Q3[4,3:5]),
  c(sum(Q3[5,1:2]), Q3[5,3:5])
)


rowSums(Q2lumped)
rowSums(Q3lumped)



