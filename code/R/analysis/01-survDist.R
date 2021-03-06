# restore session
session::restore.session("code/R/analysis/00-initialization.rda")

#load parameters
surv_parameters <- 
  "code/parameters/parametricDist.toml" %>%
  RcppTOML::parseTOML() %>%
  '[['(MODE)
  
#load functions
source("code/R/functions/fit_survival_distribution.R")

#load data
surv_data <- readxl::read_excel("data/Appendix_S1_MasterDataset.xlsx",
                                sheet = "Final Dataset")

#fit KM curve to parametric distributions
survival_distribution <- fit_survDistr(surv_data,
                                       surv_parameters$base_distr,
                                       surv_parameters$custom_distr,
                                       surv_parameters$t0,
                                       surv_parameters$t1,
                                       surv_parameters$t2)

# save session
session::save.session("code/R/analysis/01-survDist.rda")
