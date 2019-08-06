# load libraries and create some project-wide functions that will sit in the "project_functions" environment

library('tidyverse')
library('lubridate')
library('willsutils')


import::from(magrittr, "%$%")


if("project_functions" %in% search()) detach(project_functions)

with((project_functions <- new.env(parent=as.environment("package:stats"))),
  {
    #list project-specific functions here:

    lag.skipna <- function(x)
    {
      index <- dplyr::na_if(cummax(seq_along(x)*!is.na(x)), 0)
      y <- x[dplyr::lag(index)*(x==x)]
      y
    }

    distday <- function(date1, date2, maxday=21) abs(lubridate::interval(date1, date2))/lubridate::days(maxday) <= 1




    hmodel2list <- function(hmodel, hmmdist = TRUE){

      if(!hmodel$hidden) stop("hmodel.object is not a Hidden Markov Model")

      .msm.LOOKUP <- data.frame(
        label = msm:::.msm.HMODELS,
        hmmname = c("hmmCat", "hmmIdent", "hmmUnif", "hmmNorm", "hmmLNorm", "hmmExp", "hmmGamma", "hmmWeibull", "hmmPois", "hmmBinom", "hmmBetaBinom",
                    "hmmTNorm", "hmmMETNorm", "hmmMEUnif", "hmmNBinom", "hmmBeta", "hmmT",
                    "hmmClmTNorm", "hmmClmTNorm7"),
        stringsAsFactors = FALSE
      )

      # makes a state-specific vector of parameters extracted from hmodel into a list of parameters (treating hmmCat as a special case)
      makeargslist <- function(params, label){
        # params = named vector of parameters for the distribution function
        #label = label (character) of the distribution function
        if(!(label %in% .msm.LOOKUP$label)) stop("Distribution ", label, " not currently supported for hmodel2list")
        if(label=="categorical")
          list(prob = params[names(params) %in% c("p", "p0", "pbase")], basecat = params[names(params)=="basecat"])
        else if(label=="identity")
          list(x = params[names(params) == "which"])
        else
          as.list(params)
      }

      labellist <- purrr::array_branch(hmodel$labels)
      paramlist <- split(hmodel$pars, list(hmodel$parout, hmodel$parstate))
      paramnestedlist <- mapply(makeargslist, paramlist, labellist, SIMPLIFY=FALSE, USE.NAMES=FALSE)
      distlist <- lapply(labellist, function(label){match.fun(.msm.LOOKUP$hmmname[.msm.LOOKUP$label==label])})

      if(hmodel$mv){
        hmmdistlist <- purrr::invoke_map(distlist, paramnestedlist)
        hmmdistnestedlist <- split(hmmdistlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
        msmlist <- lapply(hmmdistnestedlist, function(hmmdist){purrr::lift_dl(msm::hmmMV)(hmmdist)})

        if(hmmdist)
          msmlist
        else
          split(paramnestedlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
      } else {
        if(hmmdist)
          purrr::invoke_map(distlist, paramnestedlist)
        else
          split(paramnestedlist, rep(seq_len(hmodel$nstates), times=hmodel$nout))
      }

    }





  }
)
attach(project_functions);rm(project_functions)

search()




## alternatively can use local() :
# local(
#   {
#     #list project-specific functions here:
#
#   },
#   envir = (project_functions <- new.env(parent=as.environment("package:stats")))
# )

