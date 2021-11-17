
# Note: Schaid function is broken in two parts
# Second part of this fn is in script: "Schaid_function_part2.R"

# load Rdata on exposures and outcome 
load()

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

# to run fn for all 7 adversities - exposure point previously identified using the SLCMA
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

  # load Rdata on exposures and outcome again
  load()
  
  # define covariates
  covars <- c("WHITE", "Female", "mom_age", "ppregnum", "birthweight","sustained.smoke","ed_momgest")
  
  ## of the observations available to use for adversity, narrow down to columns of interest
  df <- df[, colnames(df) %in% c("ID", covars, adver, "SMFQ_10y")]
  df <- df[complete.cases(df),]
  n <- nrow(df)
  
  # sort df here to make sure in the same order as adj.m when it was created
  df <- df %>% arrange(ID)

 # load exposure-specific adj.m previously created in "create_regressed_mvalues_2021-03-25" function
  load()
  
  
  # removing nonvariable and sex probes -----------------------------------------------------------
  cat("Removing nonvariable and sex probes \n")
  
  # check rows align with ID col in df
  identical(df$ID, rownames(adj.m))
  
  # convert matrix to data frame
  adj.m <- as.data.frame(adj.m)
  
  # create list of probes to remove
  
  # load sex probes list
  load()
  
  # load nonvariable probes list
  load()
  
  probes_to_remove <- c(nonvariable.probes, sex.probes)
  probes_to_remove <- unique(probes_to_remove)
  
  adj.m <- adj.m %>% dplyr::select(!all_of(probes_to_remove))
  
  # convert back to a matrix
  adj.m <- as.matrix(adj.m)
  
  # conditioning out confounders from analysis -----------------------------------------------------
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
  
  # regress exposures on covariates
  x <- as.numeric(as.character(x))
  fit.x <- lm(x ~ WHITE+Female+mom_age+ppregnum+birthweight+sustained.smoke+ed_momgest, data = df)
  x.resid <- fit.x$residuals
  
  # center and scale x
  x.std <- (x.resid - mean(x.resid))/sd(x.resid)
  
  # standardize med before calculating correlation - this helps values converge - same matrix as before except reduced by number of sex/nonvar probes
  med <- scale(adj.m)
  
  # check of successful med standardization
  # should get mean of 0 and sd of 1 - this is how you know data is standardized. 
  colMeans(med)  # faster version of apply(scaled.dat, 2, mean)
  apply(med, 2, sd) 
  
  # Perform SIS -------------------------------------------------------------------
  cat("Perform SIS \n")
  print(dim(x.std))
  print(dim(med))
  # finding correlations between x and meds and meds and y to identify group mediators most likely to come from
  r.x <- cor(x.std, med)
  r.y <- cor(y.std, med)
  r.xy <- abs(r.x*r.y)
  
  ## want to choose those with higher value for r.xy, which is what rank does
  rr.xy <- rank(r.xy) # with rank, 1 = lowest
  
  ## choose how many top potential mediators to analyze
  nmed.test <- floor((n-3-1)/2) # substract 2 for x and y variances
  
  ## logical, chooses “top” mediators by highest r.xy.
  is.top <- rr.xy > (length(rr.xy) - nmed.test)

    #mediators to keep
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
  filepath= "insert path here"
  save(inputs, file = paste0(filepath,timepoint,"_",Sys.Date(),".Rdata"))
  
 # Fit grid of values to find ideal lambda.vec for this data. frac.lasso =.8 based on simulation studies identifying this as ideal value
 # cat("Fitting grid \n")
  
  fit.grid <- regmed.grid(x, meds, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
  
  filepath= "insert path here"
  save(fit.grid, file = paste0(filepath,timepoint,"_",Sys.Date(),".Rdata"))

  # fit grid with lambda = 0.2 -------------------
  cat("Fitting grid with 0.2 penalty \n")
  
  # creating a fit grid for .2 to balance false positives and false negatives
  fit.grid2 <- regmed.grid(x, meds, y, lambda.vec= 0.2, frac.lasso=.8)
  
  filepath= "insert path here"
  save(fit.grid2, file = paste0(filepath,timepoint,"_",Sys.Date(),".Rdata"))
  
  cat("DONE with part one \n")
}


# This last part of the fn is contained in a second script named "Schaid_function_part2.R"
