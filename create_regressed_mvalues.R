# Function to create M-value datasets for each adversity to use the in Schaid_mediation_function script

# load Rdata containing exposures and outcome
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

# to create for all 7 adversities (function code is below)
m_values(parc, "parcruelty_8m")
m_values(onead, "oneadult_47m")
m_values(abus, "abuse_18m")
m_values(r_faminst, "r_faminstability_57m")
m_values(Fscor, "Fscore_61m")
m_values(mompsy, "mompsych_33m")
m_values(nbhq, "nbhqual_21mD")

m_values <- function(adver, timepoint){
  
  print(all_of(adver))
  print(timepoint)
  
  # loading libraries ----------------
  cat("Loading packages \n")
  
  suppressMessages(library(lumi))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyselect))
  
  # preparing outcome and exposure data for analysis ---------------------------------------------------
  
  options(digits = 7)
  
  # creating M values dataset ------------------------------------------------------------
  
  # load the betas df containing cell type and confounder data
  cat("Loading betas \n")
  filepath="define filepath here"
  load(filepath, verbose = T)
  
  ## regress out technical variables and confounders from df subsetted to each adversity's complete cases -----------------
  cat("Regress out tech vars and confounders from df SUBSET \n")
  cat(date(), "\n")
  
  # load Rdata on exposures and outcome
  load()
  
  # define covariates
  covars <- c("WHITE", "Female", "mom_age", "ppregnum", "birthweight","sustained.smoke","ed_momgest")
  
  ## of the observations available to use for adversity, narrow down to columns of interest
  df <- df[, colnames(df) %in% c("ID", covars, adver, "SMFQ_10y")]
  print(colnames(df))
  df <- df[complete.cases(df),]
  n <- nrow(df)
  print(n)
  
  # first, reduce tech_vars and betas.pheno.F7 to those IDs with complete observations for adversity of interest
  # betas.pheno.F7 contains DNAm data + tech vars + covariates w/rows by ID #
  betas.pheno.F7.adver <- betas.pheno.F7 %>% filter(ID %in% df$ID) 
  print(nrow(betas.pheno.F7.adver))
  
  # sort them both here to make sure in the same order for later regressions
  df <- df %>% arrange(ID)
  betas.pheno.F7.adver<- betas.pheno.F7.adver %>% arrange(ID)
  if (identical(df$ID, betas.pheno.F7.adver$ID)==FALSE){
    stop("df and betas.pheno.F7.adver matrices do not align by ID number.")
  }
  
  # Regression specifically for adversity's complete cases
  cat("Regress out tech and confounders from adversity-specific df \n")
  
  # creating a matrix of ALL covariates including sex, race, etc. (tech_vars + covs)
  regress_vars <- betas.pheno.F7.adver[, c("ID","Bcell","CD4T","CD8T","Gran","Mono","NK", "sample_type2", 
                                     all_of(covars))]
  str(regress_vars)
  for(i in  c("sample_type2","WHITE","Female",'sustained.smoke',"ed_momgest", "ppregnum")){
    regress_vars[,i] <- as.factor(regress_vars[,i])
  }
  print(str(regress_vars))
  print(summary(regress_vars$mom_age))
  
  betas.new <- betas.pheno.F7.adver[, c(1, 32:ncol(betas.pheno.F7))]
  print(dim(betas.new)) #652 samples, 450745 CpGs + 1 ID column
  if(identical(regress_vars$ID, betas.new$ID)==FALSE){ 
    stop("regression var and betas.pheno.F7.adver matrices do not align by ID number.")
  }
  rownames(betas.new) <- betas.new$ID
  betas.new <- betas.new[,-1]
  
  
  #regressing covariates
  beta.corr.new <- do.call(cbind, lapply(1:ncol(betas.new), function(x){
    if(x %% 1000 ==0){print(x)}
    dnam <- betas.new[,x]
    res <- residuals(lm(dnam ~ Bcell+CD4T+CD8T+Gran+Mono+NK+sample_type2+
                          WHITE+Female+mom_age+ppregnum+birthweight+
                          sustained.smoke+ed_momgest, data = regress_vars))
    adj.betas <- res + mean(dnam)
    adj.betas
  }))
  
  dim(beta.corr.new)
  colnames(beta.corr.new) <- colnames(betas.new)
  rownames(beta.corr.new) <- rownames(betas.new)
  options(digits = 22)
  beta.corr.new[which(beta.corr.new>=1)] <- max(beta.corr.new[which(beta.corr.new<1)])
  beta.corr.new[which(beta.corr.new<=0)] <- min(beta.corr.new[which(beta.corr.new>0)])
  
  # conversion of beta residuals to M values
  cat("Conversion of beta residuals to M values \n")
  adj.m <- lumi::beta2m(beta.corr.new) #convert to m-values
  
  # check no infinite values
  if (sum(is.infinite(adj.m))>0){
    stop("Error in beta to M-value conversion, infinite values found")
  }
  cat("Save df \n")
  filepath= "define filepath and filename here"
  save(adj.m, file = paste0(filepath,timepoint,"_", Sys.Date(), 
                            ".Rdata"))
  
  cat("DONE \n")
}
