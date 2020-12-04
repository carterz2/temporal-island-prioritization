# restore session
session::restore.session("code/R/analysis/01-survDist.rda")

#load parameters
candidate_parameters <- 
  "code/parameters/candidateModel.toml" %>%
  RcppTOML::parseTOML() %>%
  '[['(MODE)

#load functions
source("code/R/functions/surv_candidate_model.R")

#load data
surv_2010 <- readxl::read_excel("data/Appendix_S1_MasterDataset.xlsx",
                                sheet = "Erad 2010")
candidates_txt <- read.table("data/Candidate_Models_SLURM_Output.txt", 
                             sep = ",", header = TRUE) 
selected_distribution <- survival_distribution[4]$ig_dist

#fit candidate models and make predictions
survival_predictions <- suppressMessages(surv_candidate_model(surv_data,
                                                              selected_distibution,
                                                              candidates_txt,
                                                              candidate_parameters$number_models))
# save session
session::save.session("code/R/analysis/02-candidateModel.rda")

