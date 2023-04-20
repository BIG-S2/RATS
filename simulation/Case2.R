library(tidyverse)
library(MASS)
library(mvtnorm)
library(Ball)
library(lobstr)

args <- commandArgs(trailingOnly = TRUE)
signal_dim=as.numeric(args[1])
ri = as.numeric(args[2])
power_rep_n=as.numeric(args[3])

print(paste0("rep number: ", power_rep_n))

r= seq(1, 9, length.out = 9)


source("cal_test_statistic.R")


##########
### high dimensional sparse signal setting: normal location shift
### case 14 in summary
### generate null distribution
##########
n = m = 50
dim = 100
var_X = r[ri]
print(paste0("ss: ", n))
print(paste0("Signal dimension: ", signal_dim))
print(paste0("var: ", var_X))

### get the sample
set.seed(power_rep_n **2 + 3*power_rep_n)

if (!file.exists(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))){
  X = mvrnorm(n, rep(0, dim), diag(dim))
  X[,1:signal_dim] = mvrnorm(n, rep(0, signal_dim), diag(x = var_X, nrow = signal_dim, ncol = signal_dim))
  Y = mvrnorm(m, rep(0, dim), diag(dim))
  save(X, file=paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  save(Y, file=paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
} else {
  load(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  load(paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
}

pnull = minP.getnull(X, Y, perm = 1000)
write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull.csv"))

###########



## load the null distribution files
nperm = dim(pnull)[1]

#### sample test statistic and p values
########
val_mmd = mmd(X, Y, distance = "euclidean", median = TRUE)
val_ed = energy(X, Y)
val_ed1 = energymod(X, Y, gamma=1)
val_ed2 = energymod(X, Y, gamma=2)
val_ed3 = energymod(X, Y, gamma=3)
val_ed4 = energymod(X, Y, gamma=4)
val_ed5 = energymod(X, Y, gamma=5)
val_ed6 = energymod(X, Y, gamma=6)
val_edmax = energymod(X, Y, gamma="max")

val_bd = bd(X, Y)
val_bd1 = bdmod(X, Y, gamma=1)
val_bd2 = bdmod(X, Y, gamma=2)
val_bd3 = bdmod(X, Y, gamma=3)
val_bd4 = bdmod(X, Y, gamma=4)
val_bd5 = bdmod(X, Y, gamma=5)
val_bd6 = bdmod(X, Y, gamma=6)
val_bdmax = bdmod(X, Y, gamma="max")

pval_mmd = (sum(pnull$MMD >= val_mmd) + 1) / (nperm + 1)
pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
pval_ed1 = (sum(pnull$ED1 >= val_ed1) + 1) / (nperm + 1)
pval_ed2 = (sum(pnull$ED2 >= val_ed2) + 1) / (nperm + 1)
pval_ed3 = (sum(pnull$ED3 >= val_ed3) + 1) / (nperm + 1)
pval_ed4 = (sum(pnull$ED4 >= val_ed4) + 1) / (nperm + 1)
pval_ed5 = (sum(pnull$ED5 >= val_ed5) + 1) / (nperm + 1)
pval_ed6 = (sum(pnull$ED6 >= val_ed6) + 1) / (nperm + 1)
pval_edmax = (sum(pnull$EDMAX >= val_edmax) + 1) / (nperm + 1)
minpval_ed = min(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax)

pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
pval_bd1 = (sum(pnull$BD1 >= val_bd1) + 1) / (nperm + 1)
pval_bd2 = (sum(pnull$BD2 >= val_bd2) + 1) / (nperm + 1)
pval_bd3 = (sum(pnull$BD3 >= val_bd3) + 1) / (nperm + 1)
pval_bd4 = (sum(pnull$BD4 >= val_bd4) + 1) / (nperm + 1)
pval_bd5 = (sum(pnull$BD5 >= val_bd5) + 1) / (nperm + 1)
pval_bd6 = (sum(pnull$BD6 >= val_bd6) + 1) / (nperm + 1)
pval_bdmax = (sum(pnull$BDMAX >= val_bdmax) + 1) / (nperm + 1)
minpval_bd = min(pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax)
minpval = min(minpval_ed, minpval_bd)
##########




## permute the sample to get null distribution of min p value
minp_null = numeric(nperm)
minped_null = numeric(nperm)
minpbd_null = numeric(nperm)

for (t in 1:nperm){
  #if (t %% 100 == 0) {print(paste0("rep: ", t)); gc(verbose = TRUE)}
  
  pval_rep_ed = mean(pnull$ED >= pnull$ED[t])
  pval_rep_ed1 = mean(pnull$ED1 >= pnull$ED1[t])
  pval_rep_ed2 = mean(pnull$ED2 >= pnull$ED2[t])
  pval_rep_ed3 = mean(pnull$ED3 >= pnull$ED3[t])
  pval_rep_ed4 = mean(pnull$ED4 >= pnull$ED4[t])
  pval_rep_ed5 = mean(pnull$ED5 >= pnull$ED5[t])
  pval_rep_ed6 = mean(pnull$ED6 >= pnull$ED6[t])
  pval_rep_edmax = mean(pnull$EDMAX >= pnull$EDMAX[t])
  
  pval_rep_bd = mean(pnull$BD >= pnull$BD[t])
  pval_rep_bd1 = mean(pnull$BD1 >= pnull$BD1[t])
  pval_rep_bd2 = mean(pnull$BD2 >= pnull$BD2[t])
  pval_rep_bd3 = mean(pnull$BD3 >= pnull$BD3[t])
  pval_rep_bd4 = mean(pnull$BD4 >= pnull$BD4[t])
  pval_rep_bd5 = mean(pnull$BD5 >= pnull$BD5[t])
  pval_rep_bd6 = mean(pnull$BD6 >= pnull$BD6[t])
  pval_rep_bdmax = mean(pnull$BDMAX >= pnull$BDMAX[t])
  
  minped_null[t] = min(pval_rep_ed, pval_rep_ed1, pval_rep_ed2, pval_rep_ed3, pval_rep_ed4, 
                       pval_rep_ed5, pval_rep_ed6, pval_rep_edmax)
  minpbd_null[t] = min(pval_rep_bd, pval_rep_bd1, pval_rep_bd2, pval_rep_bd3, pval_rep_bd4,
                       pval_rep_bd5, pval_rep_bd6, pval_rep_bdmax)
  minp_null[t] = min(minped_null[t], minpbd_null[t])
}

minp_rep = tibble(minped_null, minpbd_null, minp_null)
colnames(minp_rep) = c("minP_ED", "minP_BD", "minP")
write_csv(minp_rep, paste0("./r", ri, "/n", power_rep_n, "/minprep.csv"))


## get the replicates to calculate p value of min p
minped_pval = (sum(minp_rep$minP_ED <= minpval_ed) + 1) / (nperm + 1)
minpbd_pval = (sum(minp_rep$minP_BD <= minpval_bd) + 1) / (nperm + 1)
minp_pval = (sum(minp_rep$minP <= minpval) + 1) / (nperm + 1)

var_X = r[ri]

pvec = c(var_X, pval_mmd, pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax,
         pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax,
         minped_pval, minpbd_pval, minp_pval)
names(pvec) = c("var", "MMD", "ED", "ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX",
                "BD", "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX",
                "RATS_ED", "RATS_BD", "RATS")

pvalue_out = as_tibble(t(pvec))

write_csv(pvalue_out, paste0("./r", ri, "/n", power_rep_n, "/pvalue.csv"))

###########





