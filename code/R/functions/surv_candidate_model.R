#' Select an appropriate candidate model for survival parameters with fitted
#' inverse gaussian distribution
#'
#' @param surv_data 'data.frame' excel table from carter et al.
#' 
#' @param surv_distribution 'list' output from 'fit_survival_distribution.R'
#' 
#' @param candidates_txt 'list' txt file with all candidate model combinations
#' 
#' @param nmbr_models 'integer' number of top performing models to retrieve
#' 
#' @return `list` of fitted dist results. This list contains the 
#'    following elements:
#' \describe{
#'   \item{"loglik_table"}{`data.frame` with summary statistics for model.}
#'   \item{"ig_confint_plot"}{`ggplot` with fitted inverse gaussian curves.}
#' }
#'
#' @export
surv_candidate_model <- function(
  surv_data,
  surv_distribution,
  candidates_txt,
  nmbr_models) {
  
  #assert that arguments are valid
  assertthat::assert_that(
    inherits(surv_data, "data.frame"),
    inherits(nmbr_models, "integer"))
  
  #manipulate excel version of dataset for analysis
  suppressMessages(attach(surv_data)) %>%
    na.omit(surv_data)
  surv_data <- mutate(surv_data, all = Event != 0)
  #convert variables
  surv_data <- mutate(surv_data, Mnld_Dist = log10(Mnld_Dist))
  surv_data <-mutate(surv_data, Area = log10(Area))
  surv_data <-mutate(surv_data, Buffer = (Buffer)^(1/3)) 
  surv_data <-mutate(surv_data, Landing = factor(Landing)) 
  surv_data <-mutate(surv_data, Tenure_Public = as.factor(Tenure_Public)) 
  surv_data <-mutate(surv_data, Tenure_Private = as.factor(Tenure_Private))
  surv_data <-mutate(surv_data, Tenure_Maori = as.factor(Tenure_Maori))
  surv_data <-mutate(surv_data, Nmbr_PrivOwn = 
             factor(Nmbr_PrivOwn, levels = 
                      c("none", "low", "medium", "high")))
  surv_data <-mutate(surv_data, Nmbr_MaoriOwn = 
             factor(Nmbr_MaoriOwn, levels = 
                      c("none", "low", "medium", "high")))
  
  #variable list of shortened names for convention below
  
  surv_data <- surv_data %>%
    rename(dist = Mnld_Dist,
           area = Area,
           SS = Stp_Stn,
           Buffer = Buffer,
           landing = Landing,
           public = Tenure_Public,
           private = Tenure_Private,
           maori = Tenure_Maori,
           privNumb = Nmbr_PrivOwn,
           maoriNumb = Nmbr_MaoriOwn)
  
  #step1-----------------------------------------------------------
  #Candidate model statistics
  
  #Gather top performing models from Candidate_Models_SLURM_Output.txt
  candidates <- candidates_txt
  
  #put desired number of top performing model (& parameters) in a table
  nmbr_models <- seq(1,nmbr_models)
  model_table <- data.frame()
  for (number in nmbr_models) {
    variables <- sub(".*~","",candidates[[2]][number])
    AICc <- round(candidates[[3]][number], digits = 2)
    delta_AICc <- round(candidates[[6]][number], digits = 2)
    k <- candidates[[4]][number]
    logL <- round(candidates[[5]][number], digits = 2)
    
    model_table<- rbind(model_table, c(number, variables, AICc, delta_AICc, k, logL))
  }
  colnames(model_table) <- c("candidate", 
                             "Variables",
                             "AICc", 
                             "delta_AICc",
                             "k", 
                             "logL")
    
  #Extract the top model
  surv.obj <- Surv(surv_data$Erad_Time, surv_data$all)
  candidate_model <- flexsurvreg(surv.obj ~
                                  area + SS + Buffer + landing + privNumb + maoriNumb, 
                                  data = surv_data,
                                  dist=selected_distribution)
  
  #step2-----------------------------------------------------------
  #predictions made from candidate model
  
  #Year 2050 -- 70 years
  invaded_islands <- surv_data %>%
    filter(data.frame(surv_data$Rat_Status == "Present"))
  invaded_islands <- subset(invaded_islands, select = c("Island",
                                                       "Lat",
                                                       "Lon",
                                                       "Rat_Status",
                                                       "Erad_Time",
                                                       "area",
                                                       "SS",
                                                       "Buffer",
                                                       "landing",
                                                       "privNumb",
                                                       "maoriNumb"))
  erad_prob_table <- data.frame(matrix(NA, ncol = 8))
  colnames(erad_prob_table) <- c("Island",
                                 "Lat",
                                 "Lon",
                                 "Survival_2050",
                                 "Erad_Prob_2050",
                                 "Expected_Erad_Yr",
                                 "UCL_Erad_Yr",
                                 "LCL_Erad_Yr")

  for (i in  1:nrow(invaded_islands)) {
    island_of_interest <- invaded_islands[i, ]
    surv_prediction = data.frame(summary(candidate_model,
                                         newdata = island_of_interest,
                                         t = 70))
    erad_prediction = 1-(surv_prediction[,2])
    
    #expected eradication year preds (sufficient time for success)
    island_erad_yr <- data.frame(summary(candidate_model,
                                         newdata = island_of_interest,
                                         t = c(0:450)))
    names(island_erad_yr)[1] <- "time"
    names(island_erad_yr)[2] <- "est"
    names(island_erad_yr)[3] <- "lcltime"
    names(island_erad_yr)[4] <- "ucltime"
    expected_erad_yr <- island_erad_yr[which.min(abs(0.20 -island_erad_yr$est)),][1:2]
    ucl_erad_yr <- island_erad_yr[which.min(abs(0.20 -island_erad_yr$ucltime)),][1]
    lcl_erad_yr <- island_erad_yr[which.min(abs(0.20 -island_erad_yr$lcltime)),][1]
    
    #Add relevant data from both dataframes to a single table
    erad_prob_table <- rbind(erad_prob_table,
                             list(island_of_interest$Island,
                                  island_of_interest$Lat,
                                  island_of_interest$Lon,
                                  round(surv_prediction[,2], digits = 5)*100,
                                  round(erad_prediction, digits = 5)*100,
                                  expected_erad_yr$time + 1980,
                                  ucl_erad_yr$time + 1980,
                                  lcl_erad_yr$time + 1980))
  }
  
  erad_prob_table <- erad_prob_table[order(-(erad_prob_table$Erad_Prob_2050)),]
  erad_prob_table <- slice(erad_prob_table, 
                           1:(n()-1)) #removes na col
  
  #Step3-----------------------------------------------------------
  #Determining differences in islands experiencing reinvasion and those
  # not experiencing reinvasion
  
  #Cox proportional hazards model
  cox_data <- surv_data
  cox_event <- cox_data %>%
    filter(cox_data$Event == "1")
  
  cox_model <- survfit(coxph(Surv(cox_event$Erad_Time,
                                  cox_event$Event) ~ 1 +
                         strata(cox_event$Reinvade)))
  
  cox_plot <- survminer::ggsurvplot(cox_model, 
                                    data = cox_event, 
                                    conf.int = TRUE, 
                                    conf.int.alpha = 0.3, 
                                    xlab = "Year", 
                                    ylab = "Survival Probability",
                                    legend.labs = c("not-reinvaded",
                                                    "reinvaded"), 
                                    legend = c(0.25,0.25))
  
  
  #Step4-----------------------------------------------------------
  #Assessing candidate model predictive power with 2010 dataset
  
  #2010 eradication dataset
  suppressMessages(attach(surv_2010)) %>%
    na.omit(surv_2010)
  surv_2010 <- mutate(surv_2010, all = Event != 0)
  #convert variables
  surv_2010 <-mutate(surv_2010, Mnld_Dist = log10(Mnld_Dist))
  surv_2010 <-mutate(surv_2010, Area = log10(Area))
  surv_2010 <-mutate(surv_2010, Buffer = (Buffer)^(1/3)) 
  surv_2010 <-mutate(surv_2010, Landing = factor(Landing)) 
  surv_2010 <-mutate(surv_2010, Tenure_Public = as.factor(Tenure_Public)) 
  surv_2010 <-mutate(surv_2010, Tenure_Private = as.factor(Tenure_Private))
  surv_2010 <-mutate(surv_2010, Tenure_Maori = as.factor(Tenure_Maori))
  surv_2010 <-mutate(surv_2010, Nmbr_PrivOwn = 
                       factor(Nmbr_PrivOwn, levels = 
                                c("none", "low", "medium", "high")))
  surv_2010 <-mutate(surv_2010, Nmbr_MaoriOwn = 
                       factor(Nmbr_MaoriOwn, levels = 
                                c("none", "low", "medium", "high")))
  #Rename
  surv_2010 <- surv_2010 %>%
    rename(dist_2010 = Mnld_Dist,
           area_2010 = Area,
           SS_2010 = Stp_Stn,
           Buffer_2010 = Buffer,
           landing_2010 = Landing,
           public_2010 = Tenure_Public,
           private_2010 = Tenure_Private,
           maori_2010 = Tenure_Maori,
           privNumb_2010 = Nmbr_PrivOwn,
           maoriNumb_2010 = Nmbr_MaoriOwn)
  
  #Run model
  candidate_model_2010 <- flexsurvreg(formula = Surv(Erad_Time_2010,
                                                     Event_2010) ~ 
                                        area_2010 + 
                                        SS_2010 + 
                                        Buffer_2010 + 
                                        landing_2010 + 
                                        privNumb_2010 + 
                                        maoriNumb_2010, 
                                      data = surv_2010, 
                                      dist = survival_distribution[4]$ig_dist)
  
  #prioritise for 2020
  invaded_islands_2010 <- surv_2010 %>%
    filter(surv_2010$Status_2010 == "Present")
  invaded_islands_2010 <- data.frame(invaded_islands_2010)
  
  surv_2010 <- subset(invaded_islands_2010, select = c("Island",
                                                       "Latitude",
                                                       "Longitude",
                                                       "Status_2010",
                                                       "Erad_Time_2010",
                                                       "area_2010",
                                                       "SS_2010",
                                                       "Buffer_2010",
                                                       "landing_2010",
                                                       "privNumb_2010",
                                                       "maoriNumb_2010"))
  erad_prob_2010_table <- data.frame(matrix(NA, ncol = 5))
  colnames(erad_prob_2010_table) <- c("Island",
                                      "Lat",
                                      "Lon",
                                      "Survival_2020",
                                      "Erad_prob_2020")
  
  for (i in 1:nrow(invaded_islands_2010)){
    focal_island_2010 = invaded_islands_2010[i,]
    surv_prediction_2010 <- data.frame(summary(candidate_model_2010,
                                               newdata = focal_island_2010,
                                               t = 41))
    erad_prediction_2010 <- 1-(surv_prediction_2010[,2])
    
    erad_prob_2010_table <- rbind(erad_prob_2010_table,
                                  list(focal_island_2010$Island,
                                       round(focal_island_2010$Latitude, digits = 3),
                                       round(focal_island_2010$Longitude, digits = 3),
                                       round(surv_prediction_2010[,2], digits = 3)*100,
                                       round(erad_prediction_2010, digits = 5)*100))
  }
  erad_prob_2010_table <- erad_prob_2010_table[order(rev(erad_prob_2010_table$
                                                       Erad_prob_2020)),]
  erad_prob_2010_table <- slice(erad_prob_2010_table, 1:(n()-1)) #removes na col
  
  erad_islands_2010_df <- erad_prob_2010_table %>% filter(
    Island %in% c("Taranga (Hen)",
                  "Rakitu",
                  "Wharepuaitaha",
                  "Kaihuka",
                  "Great Mercury",
                  "Harakeke (Galakek)",
                  "Horomamae (Owen)",
                  "Pakatoa",
                  "Tia Island (Entrance)",
                  "Rotoroa",
                  "Rukawahakura Island",
                  "Ulva"))
  
  #Output----------------------------------------------------------
  
  list(knitr::kable(model_table,
                    caption = "Summary table of selected candidate models."),
       knitr::kable(erad_prob_table,
                    caption = "Eradication predictions for
                    New Zealand rat-invaded islands."),
       cox_plot,
       knitr::kable(erad_islands_2010_df,
                    caption = "2010 Eradication predictions for islands
                    actually rat-eradicated 2010 - 2010."))
}

