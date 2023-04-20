## R version 4.1.3
library(tidyverse)
library(MASS)
library(mvtnorm)
library(Ball)
library(lobstr)
library(pracma)
library(data.table)
library(energy)

source("get_pnull_orig_parallel.R")
source("cal_RATS_null.R")
source("cal_RATS.R")

args <- commandArgs(trailingOnly = TRUE)
seed = as.numeric(args[1])
nperm = as.numeric(args[2])
print(paste0("Seed: ", seed))

dir = file.path("/work/users/t/i/tianyou/two_sample/UKB_Yalin/processed/776/2500_rerun/", seed)
savedir = file.path("/work/users/t/i/tianyou/two_sample/UKB_Yalin/RATS_776_2500/", seed, "perms")


### Test on the original dataset ###
# dat = fread(file.path(dir, "UKB_p12_select_hippo_resid.tsv"))
# AD = dat %>%
#   filter(AD_proxy == "Yes")
# AD_mat = as.matrix(AD[,3:dim(AD)[2]])
# CN = dat %>%
#   filter(AD_proxy == "No")
# CN_mat = as.matrix(CN[,3:dim(CN)[2]])
# test = rbind(AD_mat, CN_mat)
# eqdist.etest(test, sizes = c(dim(AD_mat)[1], dim(CN_mat)[1]), R=101)


### Test on the encoded dataset ###
dat = fread(file.path(dir, "UKB_p12_select_hippo_encoded.tsv"))
AD = dat %>%
  filter(AD_proxy == "Yes")
AD_mat = as.matrix(AD[,3:dim(AD)[2]])
CN = dat %>%
  filter(AD_proxy == "No")
CN_mat = as.matrix(CN[,3:dim(CN)[2]])
test = rbind(AD_mat, CN_mat)
eqdist.etest(test, sizes = c(dim(AD_mat)[1], dim(CN_mat)[1]), R=201)
bd.test(AD_mat, CN_mat, num.permutations = 201)

#set.seed(1024)
# pnull = lgamma.getnull(AD_mat, CN_mat, perm = 200)


X = AD_mat
Y = CN_mat

pvals_df = list()
for (perm in 1:nperm){
  pvals_df[[perm]] = read_csv(file.path(savedir, paste0("null_stats_perm",perm,".csv")))
}
pnull = bind_rows(pvals_df)
pnull = as.data.frame(pnull)

#### sample test statistic and p values
########
# pval_lst = cal.RATS.stat(X, Y, pnull)
dist = readRDS(file.path(savedir, "dist_array.rds"))
dist_ed = readRDS(file.path(savedir, "dist_ed.rds"))
gamma_seq = c(1:6, "MAX")
pval_lst = cal.RATS.stat.mat.seq(X, Y, dist, dist_ed, pnull, gamma_seq)
minpval_ed = pval_lst[[2]][1]
minpval_bd = pval_lst[[2]][2]
minpval = pval_lst[[2]][3]
##########


## permute the sample to get null distribution of min p value
minp_rep = cal.RATS.null(pnull = pnull)


## get the replicates to calculate p value of min p
minped_pval = (sum(minp_rep$minP_ED <= minpval_ed) + 1) / (nperm + 1)
minpbd_pval = (sum(minp_rep$minP_BD <= minpval_bd) + 1) / (nperm + 1)
minp_pval = (sum(minp_rep$minP <= minpval) + 1) / (nperm + 1)


pvec = c(pval_lst[[1]],
         minped_pval, minpbd_pval, minp_pval)
names(pvec) = c("ED", "ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX",
                "BD", "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX",
                "RATS_ED", "RATS_BD", "RATS")

pvalue_out = as.matrix(t(pvec))
write_csv(as_tibble(pvalue_out), file.path(savedir, "UKB_p12_select_hippo_encoded_pvals.csv"))
###########





