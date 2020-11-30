# load packages
library(sf)
library(dplyr)
library(survival)
library(flexsurv)
library(statmod)
library(invgamma)
library(actuar)
library(invgamma)

# run mode
MODE <- "debug" # for debugging the code
# MODE <- "release" # for running the code

# set RNG
#set.seed(500)

# load parameters
general_parameters <-
  "code/parameters/general.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)


# save session
session::save.session("code/R/analysis/00-initialization.rda")
