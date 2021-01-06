# Dataset and code for Carter et al. 2020 
##  The clock is ticking: Temporally prioritizing eradications on islands

### 1. To view annotated outputs, use the `main_markDown` file:
       - open the `HTML` version for tidy viewing, or
       - open the `Rmd` version to run the code in chunks

### 2. To view the inner workings of the code:
       - navigate to the `code/R/functions` folder
	   
       - open `fit_survival_distribution.R` to see information related to:
         - fitting the non-parametric kaplan-meier survival curve to the survival data
         - parametric distribution selection based on the fitted KM curve
		 
       - open `surv_candidate_model.R` to see information related to:
         - candidate model selection
         - temporal eradication predictions
         - model validation