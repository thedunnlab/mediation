## Second stage of the mediation function


# load each input & fit.grid2 & define timepoint saved from part 1
timepoint = "mompsych_33m"
timepoint = "nbhqual_21mD"
timepoint = "Fscore_61m"
timepoint = "parcruelty_8m"
timepoint = "abuse_18m"
timepoint = "oneadult_47m"
timepoint = "r_faminstability_57m"

# PART 2 of Schaid function -----------------------

# upload inputs for fit.trim2 created in part 1 of function
load()
meds <- as.matrix(inputs[,-c(1:2)]) 
dim(meds)
x <- inputs[,1]
y <- inputs[,2]

cat("Fit trim .2 \n")
# upload fit.grid2 created in part 1 of function - could also load fit.grid if using the entire grid of values
load()

fit.trim2 <- trim.best(fit.grid2, mediator.epsilon = 1e-07)
summary.regmed(fit.trim2)

# code to use for fit.trim2:
which.med <- colnames(meds) %in% dimnames(fit.trim2$alpha)[[1]]
med.selected <- meds[, which.med]
length(dimnames(fit.trim2$alpha)[[1]])

# relaxed lasso fit - model is refit without the penalization term to get unpenalized effect estimates 
fit.lam0 <- regmed.fit(x, med.selected, y, lambda = 0, frac.lasso = 0.8)
sum.fit.lam0 <- summary.regmed(fit.lam0)
write.csv(sum.fit.lam0, paste0("filepath/filename",timepoint,"_",Sys.Date(),".csv"))

# plot mediator signal
plot(fit.lam0, lty = 2, lwd = 3, cex = 0.7)

## Estimate parameter SEs ---------------------------------------------------------
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

## fit sem with lavaan - read lavaan package documentation to interpret the output
fit.lavaan <- sem(model = med.model, data = dat)
sum.fit.lavaan <- summary(fit.lavaan)
write.csv(sum.fit.lavaan, paste0("filepath/filename",timepoint,"_",Sys.Date(),".csv"))

cat("DONE \n")
