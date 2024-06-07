## Second stage of the mediation function
## 2021-03-25

# Summary of where I am in the analysis -------------------------------------------------
# At this stage, I have rerun the creation of the M-value matrices to incorporate continuous maternal age including young moms < 20 and have removed twins
# I then ran, for each of the 7 adversities, the Schaid function part 1, which conducts the S&S method through the creation of the grids
# Because previously the "ideal" lambda value for seeing mediation signal was at a lamdba = 0.2, I also ran a grid limited to 0-0.2 
# so that if I choose that value again, the fit grid is already completed. 
# Once lambda value is settled on, can rerun the grid for the smaller range if necessary and if not, can use the fit trim command to calculate effect estimates
# Today (3/25), just want to find out how many mediators are selected for each adversity and what they total. The week I'm back, I can dig into the details.

# NOTE: This script does not contain a function, it is line-by-line code to run and will change based on which lambda penalty selected.


# load each input & fit.grid2 & define timepoint saved from part 1
timepoint = "mompsych_33m"
timepoint = "nbhqual_21mD"
timepoint = "Fscore_61m"
timepoint = "parcruelty_8m"
timepoint = "abuse_18m"
timepoint = "oneadult_47m"
timepoint = "r_faminstability_57m"

## Code to check if 0.2 still the best lambda penalty using data from 3/25/21 update --------------
cat("Fit trim full \n")
fit.trim <- trim.best(fit.grid, mediator.epsilon = 1e-07)
summary.regmed(fit.trim) ## HERE check and see what it looks like best lambda value is. Choose penalty based on when betas aren't all 0s.

test8 <- cbind(fit.grid[["fit.list"]][[8]][["alpha"]],fit.grid[["fit.list"]][[9]][["beta"]])
test9 <- cbind(fit.grid[["fit.list"]][[9]][["alpha"]],fit.grid[["fit.list"]][[9]][["beta"]])
test10 <- cbind(fit.grid[["fit.list"]][[10]][["alpha"]],fit.grid[["fit.list"]][[9]][["beta"]])

# Again, on the whole, 0.2 best lambda value

# PART 2 of Schaid function -----------------------

# upload inputs for fit.trim2 created in part 1 of function
load(paste0("/data/js95/DunnLab/Brooke/Mediation/data/med_step2_inputs/inputs_",timepoint,"_", 
            Sys.Date()-1,".Rdata"))
meds <- as.matrix(inputs[,-c(1:2)]) 
dim(meds)
x <- inputs[,1]
y <- inputs[,2]

cat("Fit trim .2 \n")
# upload fit.grid2 created in part 1 of function
load(paste0("/data/js95/DunnLab/Brooke/Mediation/data/fit.grid2s_SMFQ_10y/fit.grid2_",timepoint,"_2021-03-26",".Rdata"))

fit.trim2 <- trim.best(fit.grid2, mediator.epsilon = 1e-07)
summary.regmed(fit.trim2) # 9 mediators, this checks out.

# code to use for fit.trim2:
which.med <- colnames(meds) %in% dimnames(fit.trim2$alpha)[[1]]
med.selected <- meds[, which.med]
length(dimnames(fit.trim2$alpha)[[1]])

# relaxed lasso fit - model is refit without the penalization term to get unpenalized effect estimates 
fit.lam0 <- regmed.fit(x, med.selected, y, lambda = 0, frac.lasso = 0.8)
sum.fit.lam0 <- summary.regmed(fit.lam0)
write.csv(sum.fit.lam0, paste0("/data/js95/DunnLab/Brooke/Mediation/output/summ.fit.lam0.",timepoint,"_",Sys.Date(),".csv"))

# plot mediator signal, if desired
plot(fit.lam0, lty = 2, lwd = 3, cex = 0.7)

## Code to use to estimate parameter SEs ---------------------------------------------------------
cat("Fit SEs \n") 
## choose subset of mediators that are in fit.trim
mediator.names <- dimnames(med.selected)[[2]]

## setup lavaan model
med.model <- lavaan.model(y.name = "y", x.name = "x", med.name = mediator.names, 
                          medcov = fit.trim2$MedCov)
## note: fit.trim2 is the same as fit.lam0

## setup data for lavaan
dat <- data.frame(cbind(scale(x, center = TRUE, scale = TRUE), scale(y, center = TRUE, 
                                                                     scale = TRUE), scale(med.selected, center = TRUE, scale = TRUE)))
names(dat) <- c("x", "y", dimnames(med.selected)[[2]])

## fit sem with lavaan - read lavaan documentation to interpret the output
fit.lavaan <- sem(model = med.model, data = dat)
sum.fit.lavaan <- summary(fit.lavaan)
write.csv(sum.fit.lavaan, paste0("/data/js95/DunnLab/Brooke/Mediation/output/summ.fit.lavaan.",timepoint,"_",Sys.Date(),".csv"))

cat("DONE \n")
