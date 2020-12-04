#' fit parametric survival distributions to kaplan-meier survival curve
#'
#' @param surv_data 'data.frame' excel table from carter et al.
#' 
#' @param base_distr 'character' relevant distributions provided in R packages
#' 
#' @param custom_distr 'character' custom distributions created for analysis
#' 
#' @param t0 'integer' start time (yrs), 0 = 1980
#' 
#' @param t1 'integer' mid-point time of interest, 40 = 2020 (current day)
#' 
#' @param t2 'integer' end time (yrs), 70 = 2050
#' 
#' @return `list` of fitted dist results. This list contains the 
#'    following elements:
#' \describe{
#'   \item{"loglik_table"}{`data.frame` with summary statistics for model.}
#'   \item{"ig_confint_plot"}{`ggplot` with fitted inverse gaussian curves.}
#' }
#'
#' @export
#' 
fit_survDistr <- function(
  surv_data,
  base_distr,
  custom_distr,
  t0,
  t1,
  t2) {
  
  #assert that arguments are valid
  assertthat::assert_that(
    inherits(surv_data, "data.frame"),
    is.character(base_distr),
    is.character(custom_distr))
  #manipulate excel version of dataset for analysis
  suppressMessages(attach(surv_data)) %>%
    na.omit(surv_data)
  surv_data <- mutate(surv_data, all = Event != 0)
  
  # step 1-----------------------------------------------------------
  
  #Fit non-parametric KM survival curve to dataset
  surv.obj <- Surv(surv_data$Erad_Time, surv_data$all)
  fit_km <- survfit(surv.obj ~ 1,
                    data = surv_data)
  km_curve <-fit_km$surv
  
  #Step 2-----------------------------------------------------------
  
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
  
  #Step 3-----------------------------------------------------------
  
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
  invgaus_curve <- curve(1-pinvgauss(x, invgaus_param1,
                                     dispersion=invgaus_param2),
                         from=t0,to=t2,ylim=c(0,1))
  dev.off()
  invgaus_curve_df <- data.frame(invgaus_curve)
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
    parameter3 = c(NA, NA,NA,NA,
                   NA,NA,NA,NA,
                   burr_param3)
    )

  #Step 4-----------------------------------------------------------
  
  #parametric bootstrap of selected inverse gaussian curve for CI's
  #Translated in to code from Efran and Tibshirani (1993)
  
  #v-cov matrix and coefficient estimates of invgaus curve
  vcov_invgaus <- vcov(invgaus_surv)
  surv_sigma <- matrix(c(vcov_invgaus), 2)
  surv_mu <- c(invgaus_param1, invgaus_param2)
  
  #create 250 iterations of different curves based on v-cov
  # & mean/dispersion
  set.seed(1)
  seed <- .Random.seed
  
  #keep the randomly generated numbers between instances
  rand_restore <- function() {
    if (exists(".Random.seed", .GlobalEnv))
      seed <- .GlobalEnv$.Random.seed
    else
      seed <- NULL
    
    if (!is.null(seed)) 
      .GlobalEnv$.Random.seed <- seed
    else
      rm(".Random.seed", envir = .GlobalEnv)
  }
  confint_preds<- MASS::mvrnorm(n = 250, mu = surv_mu,
                                Sigma = surv_sigma)
  colnames(confint_preds) <- c("mean", "dispersion")
  
  #Remove negative dispersion values by adding positive values 
  # to a separate dataframe
  surv_curv_df <- data.frame(mean = numeric(0),
                             dispersion = numeric(0))
  rows <- 1:nrow(confint_preds)
  for (row in rows){
    surv_disp = confint_preds[row,2]
    if(surv_disp > 0){
      bootstrap_val <- confint_preds[row,2]
      surv_curv_df <- rbind(surv_curv_df,
                            confint_preds[row, ])
    }
    
    colnames(surv_curv_df) <- c("mean", "dispersion")
  }
  
    #Create different lines as confidence intervals &
    #add to invgaus curve
    combined_curves <- data.frame()
    curv_numb <- 1:100 #100 lines
    invisible(for(numb in curv_numb){
      curv_numb<- numb
      surv_mean <- surv_curv_df[numb,1]
      surv_disp <- surv_curv_df[numb,2]
      bootstrap_curv <- curve(1-pinvgauss(x, 
                                          mean = surv_mean,
                                          dispersion = surv_disp),
                              from=0,to=100,ylim=c(0,1),
                              type = "l")
      dev.off() #dont plot 100 F*cking curves R.
      bootstrap_curv_df <- data.frame(curv_numb, bootstrap_curv)
      combined_curves <- rbind(combined_curves, bootstrap_curv_df)
    })
    
    #sort bootstrap curves and select 95% confidence intervals
    # at t = 40 (current day)
    cie95_df <- combined_curves[combined_curves$x == t1,]
    cie95_df <- cie95_df[order(cie95_df$y),]
    percentile95 <- quantile(cie95_df$y, c(0.05, 0.95))
    
    #find curves closest to 95% CIE values
    cie95_low_curve <- cie95_df[findInterval(percentile95[1], cie95_df$y),]
    cie95_high_curve <- cie95_df[findInterval(percentile95[2], cie95_df$y),]
    
    #Isolate curves from df
    conf_bound_curve1 <- combined_curves[combined_curves$curv_numb == 
                                           cie95_low_curve$curv_numb,]
    conf_curve1_df <- data.frame(conf_bound_curve1$x, conf_bound_curve1$y)
    conf_bound_curve2 <- combined_curves[combined_curves$curv_numb == 
                                           cie95_high_curve$curv_numb,]
    conf_curve2_df <- data.frame(conf_bound_curve2$x, conf_bound_curve2$y)
    
    #Confidence curves from publication
    fit.ig.CIE5 <- curve(1-pinvgauss(x, 73.84748,dispersion=0.2875484688),
                         from=t0,to=t2,ylim=c(0,1))
    dev.off()
    fit.ig.CIE5.df <- data.frame(fit.ig.CIE5)
    fit.ig.CIE95 <- curve(1-pinvgauss(x, 80.20728,dispersion=0.0101536839),
                          from=t0,to=t2,ylim=c(0,1))
    dev.off()
    fit.ig.CIE95.df <- data.frame(fit.ig.CIE95)
    
    cie_table <- data.frame(
      conf_Bounds = c("5%","95%"),
      parameter1 = c(73.84748,80.20728),
      parameter2 = c(0.2875, 0.0101)
    )
    
  #Output----------------------------------------------------------
  
    ig_confint_plot<- ggplot() +
      geom_line(data = fit.ig.CIE5.df, aes(x = x, y = y), colour = 'red') +
      geom_line(data = fit.ig.CIE95.df, aes(x = x, y = y), colour = 'red') +
      geom_line(data = invgaus_curve_df, aes(x = x, y = y), colour = 'black')+
      labs(x = "Year", y = "Survival probability")+
      scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70),
                         labels = c(1980,1990,2000,2010,2020,2030,2040,2050)) +
      theme_classic()  
    
  #loglikelihood statistics table
    
    list(knitr::kable(loglik_table,
                      caption = "summary table of parametric distributions fitted to KM curve."),
         knitr::kable(cie_table,
                      caption = "95-5% CIE for inverse gaussian distribution\nbased on bootstrapping."),
         inverse_gaussian_plot = ig_confint_plot,
         ig_dist = ig)
}

