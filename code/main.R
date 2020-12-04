message("starting 00")
source("code/R/analysis/00-initialization.R")

message("\nstarting 01")
source("code/R/analysis/01-survDist.R")
message("use 'survival_distribution' to view results.")

message("starting 02")
source("code/R/analysis/02-candidateModel.R")
message("use 'survival_predictions' to view results.")

#source("code/R/analysis/03-exports.R")


