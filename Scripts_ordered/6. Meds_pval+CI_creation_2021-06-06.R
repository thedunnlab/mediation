## Mediation Results
## Calculating p-value, CI, and finding closest gene for each site

# import med data - data from Schaid fn compiled into a single Excel doc
library(readxl)
med.data <- read_excel("/data/js95/DunnLab/Brooke/Mediation/data/Mediation_sites_all_2021-06-04.xls")
med.data <- med.data[-71,] # blank row

## Code from quantpsy.org to calculate a p-value -------
require(MASS)
a=-0.09671372
b=-0.1274389
rep=20000
conf=95
pest=c(a,b)
acov <- matrix(c(
  0.001496979, 0,
  0, 0.001361938
),2,2)
mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
ab <- mcmc[,1]*mcmc[,2]
low=(1-conf/100)/2
upp=((1-conf/100)/2)+(conf/100)
LL=quantile(ab,low)
UL=quantile(ab,upp)
LL4=format(LL,digits=4)
UL4=format(UL,digits=4)
p0 <- ecdf(ab)(0)
P=1-2*abs(p0-0.5)
P

# Converting p-value code to a function and adding pvalue to med.data df

pval.calc <- function(a, b, var.a, var.b) {
  require(MASS)
  set.seed(123)
  rep=20000
  conf=95
  pest=c(a,b)
  acov <- matrix(c(
    var.a, 0,
    0, var.b
  ),2,2)
  mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
  ab <- mcmc[,1]*mcmc[,2]
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(ab,low)
  UL=quantile(ab,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)
  p0 <- ecdf(ab)(0)
  P=1-2*abs(p0-0.5)
  P
}

pval.calc(a, b, 0.001496979, 0.001361938) # test, works!

# converted to a fn with help from Alex
pval <- unlist(lapply(1:nrow(med.data), function(i){
  dat <- med.data[i,]
  pval.calc(a=dat$alpha, b=dat$beta, var.a = dat$var.a, var.b = dat$var.b)
}))
med.data$pvalue <- pval

# calculating confidence intervals
LL.calc <- function(a, b, var.a, var.b) {
  require(MASS)
  set.seed(123)
  rep=20000
  conf=95
  pest=c(a,b)
  acov <- matrix(c(
    var.a, 0,
    0, var.b
  ),2,2)
  mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
  ab <- mcmc[,1]*mcmc[,2]
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(ab,low)
  UL=quantile(ab,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)
  p0 <- ecdf(ab)(0)
  P=1-2*abs(p0-0.5)
  LL4
}

UL.calc <- function(a, b, var.a, var.b) {
  require(MASS)
  set.seed(123)
  rep=20000
  conf=95
  pest=c(a,b)
  acov <- matrix(c(
    var.a, 0,
    0, var.b
  ),2,2)
  mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
  ab <- mcmc[,1]*mcmc[,2]
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(ab,low)
  UL=quantile(ab,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)
  p0 <- ecdf(ab)(0)
  P=1-2*abs(p0-0.5)
  UL4
}
 
## add CI to data frame
LL4 <- unlist(lapply(1:nrow(med.data), function(i){
  dat <- med.data[i,]
  LL.calc(a=dat$alpha, b=dat$beta, var.a = dat$var.a, var.b = dat$var.b)
}))
med.data$LL <- LL4

UL4 <- unlist(lapply(1:nrow(med.data), function(i){
  dat <- med.data[i,]
  UL.calc(a=dat$alpha, b=dat$beta, var.a = dat$var.a, var.b = dat$var.b)
}))
med.data$UL <- UL4

# add nearest gene data --------
price.annot <- read.table("/data/js95/DunnLab/alussier/price_annotation_450k.txt", 
                          header=T, sep='\t')

annot.hits <- price.annot[which(price.annot$ID %in% med.data$CpG),]
dim(annot.hits)
annot.hits <- annot.hits[match(med.data$CpG, as.character(annot.hits$ID)),]
identical(as.character(annot.hits$ID), as.character(med.data$CpG))
#med.data[,c("CpG","Adversity")]

med.data$Gene <- annot.hits$Closest_TSS_gene_name

# save changes
#save(med.data, file = "/data/js95/DunnLab/Brooke/Mediation/data/Mediation_sites_2021-06-06.Rdata")
#write.csv(med.data, file = "/data/js95/DunnLab/Brooke/Mediation/data/Mediation_sites_2021-06-06.csv")