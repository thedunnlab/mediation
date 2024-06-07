# Updated fn that removes non variable probes (as of 10/07/20 update) and removes non-whites and twins (as of 1/25/21 update)
# updated again on 2/5/21 to reincorporate non-whites into the analysis
# updated on 2/24/21 to use newest df that removes young moms under 20, uses a binary variable for mom age, and removes twins
# updated on 3/25/21 to change maternal age to a continuous variable - mom_age - and reincorporate young moms
# Note: Have broken fn into two parts because script with entire function will not run all the way through
# Second part of this fn is in script: "5. Schaid_function_part2_2021-03-26.R"

# load data to use
load("/data/js95/DunnLab/Brooke/Mediation/data/mediation_data_2021-06-03.Rdata")

## defining adversities
## adversity measures before 7
parc <- grep("^parc", names(df), value=T)
#parc <- parc[parc != "parcruelty_9y"] 
onead <- grep("^onead", names(df), value=T)
#onead <- onead[onead != "oneadult_8y"] # remove year 8 estimate
abus <- grep("^abus", names(df), value=T)
#abus <- abus[abus != "abuse_8y"] # remove year 8 estimate 
r_faminst <- grep("^r_faminst", names(df), value=T)
Fscor <- grep("^Fscor", names(df), value=T)
mompsy <- grep("^mompsy", names(df), value=T)
nbhq <- grep("^nbhqual", names(df), value=T)

# to run fn for all 7 adversities
schaid_mediation_reduc_probes_p1(parc, "parcruelty_8m")
schaid_mediation_reduc_probes_p1(onead, "oneadult_47m")
schaid_mediation_reduc_probes_p1(abus, "abuse_18m")
schaid_mediation_reduc_probes_p1(r_faminst, "r_faminstability_57m")
schaid_mediation_reduc_probes_p1(Fscor, "Fscore_61m")
schaid_mediation_reduc_probes_p1(mompsy, "mompsych_33m")
schaid_mediation_reduc_probes_p1(nbhq, "nbhqual_21mD")


schaid_mediation_reduc_probes_p1 <- function(adver, timepoint){
  
  print(adver)
  print(timepoint)
  
  # loading libraries ----------------
  cat("Loading packages \n")
  
  suppressMessages(library(lumi))
  suppressMessages(library(knitr))
  suppressMessages(library(dplyr))
  suppressMessages(library(regmed))
  library(lavaan)
  
  # regmed will get you the selected mediation sites
  # lavaan will get you the SEs for those estimates
  
  # preparing outcome and exposure data for analysis ---------------------------------------------------
  
  options(digits = 7)
  
  # defining df ----------------------------------------------------------------------------------
  # data is loaded in the file where this function runs.
  
  # load data
  load("/data/js95/DunnLab/Brooke/Mediation/data/mediation_data_2021-06-03.Rdata")
  
  covars <- c("WHITE", "Female", "mom_age", "ppregnum", "birthweight","sustained.smoke","ed_momgest")
  
  ## of the observations available to use for adversity, narrow down to columns of interest
  df <- df[, colnames(df) %in% c("ID", covars, adver, "SMFQ_10y")]
  df <- df[complete.cases(df),]
  n <- nrow(df)
  
  # sort df here to make sure in the same order as adj.m when it was created
  df <- df %>% arrange(ID)

  # load appropriate adj.m previously created in "3. create_regressed_mvalues_2021-03-25" function
  load(paste0("/data/js95/DunnLab/Brooke/Mediation/data/Adjusted_M-value_dataframes/adjusted_mvalues_w_covs_",timepoint,"_2021-03-25.Rdata"))

  
  # removing nonvariable and sex probes -----------------------------------------------------------
  cat("Removing nonvariable and sex probes \n")
  
  # check rows align with ID col in df
  identical(df$ID, rownames(adj.m))
  
  # convert matrix to data frame
  adj.m <- as.data.frame(adj.m)
  
  # create list of probes to remove
  
  # load sex probes list
  load("/data/js95/ALSPAC/ARIES/DNAm_2020/F7/sex.probes.list.Rdata")
  
  # load nonvariable probes
  load("/data/js95/ALSPAC/ARIES/DNAm_2020/F7/nonvariable.probes.list.Rdata")
  
  probes_to_remove <- c(nonvariable.probes, sex.probes)
  probes_to_remove <- unique(probes_to_remove) # find unique probes among combined list
  # remove probes from adj.m
  adj.m <- adj.m %>% dplyr::select(!all_of(probes_to_remove))
  
  # convert back to a matrix
  adj.m <- as.matrix(adj.m)
  
  # conditioning out confounders from analysis -----------------------------------------------------
  # all model inputs must be regressed on confounders prior to analysis
  # all model inputs must be standardized prior to analysis to help with later interpretation of results
  cat("conditioning out confounders from outcome and exposures \n")
  
  ## let x be the exposure variable and med be the potential mediators
  x <- df[,timepoint]
  y.response <- df$SMFQ_10y
  # med matrix is M-values after regressing on technical and confounding covs
  
  ## regess y.response on covariates to get residuals as y-adjusted outcome
  fit <- lm(y.response ~ WHITE+Female+mom_age+ppregnum+birthweight+sustained.smoke+ed_momgest, data = df)
  y.resid <- fit$residuals
  
  ## center and scale y - standardization
  y.std <- (y.resid - mean(y.resid))/sd(y.resid)
  
  ## PURIST APPROACH
  # regress exposures on covariates
  x <- as.numeric(as.character(x))
  fit.x <- lm(x ~ WHITE+Female+mom_age+ppregnum+birthweight+sustained.smoke+ed_momgest, data = df)
  x.resid <- fit.x$residuals
  
  # center and scale x
  x.std <- (x.resid - mean(x.resid))/sd(x.resid)
  
  # standardize med before calculating correlation - this helps values converge - same matrix as before except reduced by number of sex/nonvar probes
  med <- scale(adj.m)
  
  # check if successful med standardization
  # should get mean of 0 and sd of 1 - this is how you know data is standardized. These values were once normalized,
  # but that's not the same thing. normalized means it was scaled to a value between 0 and 1.
  colMeans(med)  # faster version of apply(scaled.dat, 2, mean)
  apply(med, 2, sd) 
  
  # Perform SIS -------------------------------------------------------------------
  cat("Perform SIS \n")
  print(dim(x.std))
  print(dim(med))
  # finding correlations between x and meds and meds and y to identify CpGs mediation is most likely to come from
  r.x <- cor(x.std, med)
  r.y <- cor(y.std, med)
  r.xy <- abs(r.x*r.y)
  
  ## want to choose those with higher value for r.xy, which is what rank does
  rr.xy <- rank(r.xy) # with rank, 1 = lowest
  
  ## choose how many top potential mediators to analyze
  nmed.test <- floor((n-3-1)/2) # -2 for x and y variances
  
  #test, looking at top 10,000
  #top.10000 <- rr.xy > (length(rr.xy) - 10000)
  
  ## logical, chooses “top” mediators by highest r.xy.
  is.top <- rr.xy > (length(rr.xy) - nmed.test)

  # mediators to keep
  medkeep <- med[, is.top]
  
  
  ## Implementing regmed package ------------------------------------------------------------
  cat("Implementing regmed package \n")
  
  # renaming inputs to make commands/typing more concise
  meds <- medkeep
  x <- x.std
  y <- y.std
  
  # reformatting before saving
  med.df <- as.data.frame(medkeep)
  # saving inputs for later manipulation
  inputs <- cbind(x, y, med.df)
  save(inputs, file = paste0("/data/js95/DunnLab/Brooke/Mediation/data/med_step2_inputs/inputs_",timepoint,"_", 
                             Sys.Date(),".Rdata"))
  

  # load inputs data
  #load(paste0("/data/js95/DunnLab/Brooke/Mediation/data/med_step2_inputs/inputs_",timepoint,"_",Sys.Date(),".Rdata")) ### THIS IS ALWAYS WHERE THE FUNCTION BREAKS DOWN. GLOBAL ENVIR?
  
  ## ideally would like to be able to remove the next four lines since these are just defined in lines 124:126
 # load(paste0("/data/js95/DunnLab/Brooke/Mediation/data/med_step2_inputs/inputs_",timepoint,"_", 
 #             Sys.Date(),".Rdata"))
 # meds <- as.matrix(inputs[,-c(1:2)]) 
 # dim(meds)
 # x <- inputs[,1]
 # y <- inputs[,2]
  
  ##### Due to changes in df structure have to look at ENTIRE grid again ###  ------------
  ## We will use a lambda = 0.2 penalty for all adversities going forward.
  # Fit grid of values to find ideal lambda.vec for this data. frac.lasso =.8 based on simulation studies id-ing this as ideal value
 cat("Fitting entire grid \n")
  
 fit.grid <- regmed.grid(x, meds, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
  
 save(fit.grid, file = paste0("/data/js95/DunnLab/Brooke/Mediation/fit.grids_SMFQ_10y/fit.grid_",timepoint,"_", 
                              Sys.Date(),".Rdata"))
  # using this grid can see if 0.2 still the ideal lambda value. Repeat this process for all 7 adversities
  
  # fit grid with lambda = 0.2 -------------------
  # fitting grid of values to find ideal lambda.vec for this data. lasso =.8 based on simulation studies id-ing this as ideal value
  cat("Fitting grid with 0.2 penalty \n")
  
  # as every grid showed mediation signal at .2, creating a fit grid for .2
  fit.grid2 <- regmed.grid(x, meds, y, lambda.vec= 0.2, frac.lasso=.8)
  
  save(fit.grid2, file = paste0("/data/js95/DunnLab/Brooke/Mediation/data/fit.grid2s_SMFQ_10y/fit.grid2_",timepoint,"_", 
                               Sys.Date(),".Rdata"))
  
  cat("DONE with part one \n")
}


### Can't find any way to automate the last part of the function so must do manually for each adversity
# This last part of the fn is contained in a new script named "5. Schaid_function_part2_2021-03-26.R"
