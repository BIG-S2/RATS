library(tidyverse)
library(MASS)
library(mvtnorm)
library(Ball)
library(lobstr)
library(pracma)

args <- commandArgs(trailingOnly = TRUE)
ri = as.numeric(args[2])
power_rep_n=as.numeric(args[3])
rot = args[4]

r = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

print(paste0("rep number: ", power_rep_n))


source("get_pnull_orig.R")
source("get_pnull_rotation.R")
source("cal_RATS_null.R")
source("cal_RATS.R")


generate.rotation = function(d, n = 10){
  rot = vector("list", length = n)
  for (i in 1:n){
    rot[[i]] = randortho(d, type = "orthonormal")
  }
  return(rot)
}



##########
### Bernoulli distribution case from Zhu et al
### generate null distribution
##########
n = m = 50
dim = 100
dim_X = r[ri]
print(paste0("ss: ", n))
print(paste0("Signal dimension: ", dim_X))

### get the sample
set.seed(power_rep_n **2 + 3*power_rep_n)

if (!file.exists(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))){
  orig.X = matrix(rbinom(n * dim, 1, 0.5), ncol = dim, nrow = n)
  orig.Y = matrix(rbinom(n * dim, 1, 0.5), ncol = dim, nrow = n)
  if (dim_X >= 2){
    orig.X[, seq(2, dim_X, by = 2)] = orig.X[, seq(1, dim_X-1, by = 2)]
  }
  save(orig.X, file=paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  save(orig.Y, file=paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
} else {
  load(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  load(paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
}

if (rot == "rotation"){
  print("Performed rotation")
  A.list = generate.rotation(dim, 25)
  save(A.list, file=paste0("./r", ri, "/n", power_rep_n, "/rotation_list_", power_rep_n, ".RData"))
  #rotated = find.max.rotation.from.list(orig.X, orig.Y, A.list)
  #X = rotated[[1]]
  #Y = rotated[[2]]
  pnull = lgamma.getnull.rot.from.list(orig.X, orig.Y, A.list, perm = 500)
} else {
  print("Original version")
  X = orig.X
  Y = orig.Y
  pnull = lgamma.getnull(X, Y, perm = 500)
}

write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull_", rot, ".csv"))

###########


#### sample test statistic and p values
nperm = dim(pnull)[1]
if (rot == "rotation"){
  pval_lst = cal.RATS.stat.rot(X = orig.X, Y = orig.Y, A.list, pnull)
} else {
  pval_lst = cal.RATS.stat(orig.X, orig.Y, pnull)
}
minpval_ed = pval_lst[[2]][1]
minpval_bd = pval_lst[[2]][2]
minpval = pval_lst[[2]][3]

minpval_ed_simp = pval_lst[[3]][1]
minpval_bd_simp = pval_lst[[3]][2]
minpval_simp = pval_lst[[3]][3]

## permute the sample to get null distribution of min p value
minp_rep = cal.RATS.null(pnull = pnull)
write_csv(minp_rep, paste0("./r", ri, "/n", power_rep_n, "/minprep_", rot, ".csv"))


## get the replicates to calculate p value of min p
minped_pval = (sum(minp_rep$minP_ED <= minpval_ed) + 1) / (nperm + 1)
minpbd_pval = (sum(minp_rep$minP_BD <= minpval_bd) + 1) / (nperm + 1)
minp_pval = (sum(minp_rep$minP <= minpval) + 1) / (nperm + 1)

minped_simp_pval = (sum(minp_rep$minP_ED_simp <= minpval_ed_simp) + 1) / (nperm + 1)
minpbd_simp_pval = (sum(minp_rep$minP_BD_simp <= minpval_bd_simp) + 1) / (nperm + 1)
minp_simp_pval = (sum(minp_rep$minP_simp <= minpval_simp) + 1) / (nperm + 1)

pvec = c(dim_X, pval_lst[[1]],
         minped_pval, minpbd_pval, minp_pval,
         minped_simp_pval, minpbd_simp_pval, minp_simp_pval)
names(pvec) = c("dim", "ED", "ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX",
                "BD", "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX",
                "RATS_ED", "RATS_BD", "RATS", "RATS_ED_simp", "RATS_BD_simp", "RATS_simp")

pvalue_out = as_tibble(t(pvec))

write_csv(pvalue_out, paste0("./r", ri, "/n", power_rep_n, "/pvalue_", rot, ".csv"))

###########





