#' fit parametric survival distributions to kaplan-meier survival curve
#'
#' @param surv_data 'data.frame' excel table from carter et al.
#' 
#' @param base_distr 'character' relevant distributions provided in R packages
#' 
#' @param custom_distr 'character' custom distributions created for analysis
#' 
#' @return `list` of fitted dist results. This list contains the 
#'    following elements:
#' \describe{
#'   \item{"loglik_table"}{`data.frame` with summary statistics for model.}
#' }
#'
#' @export
#' 
fit_survDistr <- function(
  surv_data,
  base_distr,
  custom_distr) {
  
  #assert that arguments are valid
  assertthat::assert_that(
    inherits(surv_data, "data.frame"),
    is.character(base_distr),
    is.character(custom_distr))
  #manipulate excel version of dataset for analysis
  suppressMessages(attach(surv_data)) %>%
    na.omit(surv_data)
  surv_data <- mutate(surv_data, all = Event != 0)
  
  #------------------------------------------------------------------
  
  #Fit non-parametric KM survival curve to dataset
  surv.obj <- Surv(surv_data$Erad_Time, surv_data$all)
  fit_km <- survfit(surv.obj ~ 1,
                    data = surv_data)
  km_curve <-fit_km$surv
  
  #------------------------------------------------------------------
  
  #fit basic parametric distributions to KM curve
  base_dists <- plyr::llply(
    seq_len(length(surv_parameters$base_distr)), function(i) {
      flexsurvreg(surv.obj ~ 1,
                  data = surv_data,
                  dist = base_distr[i])
    })
  
  #Extract loglikelihoods
  base_logvals <- plyr::llply(
    seq_len(length(surv_parameters$base_distr)), function(i) {
      base_dists[[i]]$loglik
    })
  
  #Extract shape and dispersion values
  base_params <- plyr::llply(
    seq_len(length(surv_parameters$base_distr)), function(i) {
      c(base_dists[[i]]$cov[1:2])
    })
  
  #------------------------------------------------------------------
  
  #fit custom parametric distributions to km curve
  #inverse gaussian distribution
  ig <- list(name=custom_distr[1], pars=c("mean", "dispersion"), 
             location="mean",
             transforms=c(I,log), 
             inv.transforms=c(I, exp),
             inits=function(t){ c((median(t)), log(median(t))) })
  invgaus_surv <- suppressMessages(flexsurvreg(surv.obj ~ 1, 
                                               data = surv_data,
                                               ci= TRUE,
                                               dist=ig))
  invgaus_param1 <- invgaus_surv$res[1]
  invgaus_param2 <- invgaus_surv$res[2]
  #custom: inverse gamma distribution
  invgam <- list(name=custom_distr[2], pars = c("shape", "scale"),
                 location = "scale",
                 transforms =c(log, log),
                 inv.transforms = c(exp, exp),
                 inits = function(t){c(1, (median(t))) })
  invgam_surv <- suppressMessages(flexsurvreg(surv.obj ~ 1,
                                              data = surv_data,
                                              dist=invgam))
  invgam_param1 <- invgam_surv$res[1]
  invgam_param2 <- invgam_surv$res[2]
  #custom: burr distribution
  burr <- list(name = custom_distr[3], 
               pars = c("shape1", "shape2", "scale"),
               location = "scale",
               transforms = c(log, log, log),
               inv.transforms = c(exp, exp, exp),
               inits = function(t){c(log(median(t)),
                                     log(median(t)),
                                     median(t)) })
  burr_surv <- suppressMessages(flexsurvreg(surv.obj ~ 1,
                                            data = surv_data,
                                            dist=burr))
  burr_param1 <- burr_surv$res[1]
  burr_param2 <- burr_surv$res[2]
  burr_param3 <- burr_surv$res[3]
  # compile loglikelihoods
  loglik_table <- data.frame(
    distribution = c(base_distr,custom_distr),
    loglikelihood = c(base_logvals[[1]],
                      base_logvals[[2]],
                      base_logvals[[3]],
                      base_logvals[[4]],
                      base_logvals[[5]],
                      base_logvals[[6]],
                      invgaus_surv$loglik,
                      invgam_surv$loglik,
                      burr_surv$loglik),
    parameter1 = c(base_params[[1]][1],
                   base_params[[2]][1],
                   base_params[[3]][1],
                   base_params[[4]][1],
                   base_params[[5]][1],
                   base_params[[6]][1],
                   invgaus_param1,
                   invgam_param1,
                   burr_param1),
    parameter2 = c(base_params[[1]][2],
                   base_params[[2]][2],
                   base_params[[3]][2],
                   base_params[[4]][2],
                   base_params[[5]][2],
                   base_params[[6]][2],
                   invgaus_param2,
                   invgam_param2,
                   burr_param2),
    parameter3 = c(NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   NA,
                   burr_param3)
    )
  #print output for posterity
 knitr::kable(loglik_table,
              caption = "summary table of parametric distributions fitted to KM curve.")

}
