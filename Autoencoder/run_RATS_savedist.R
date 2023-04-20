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
print(paste0("Seed: ", seed))

dir = file.path("/work/users/t/i/tianyou/two_sample/UKB_Yalin/processed/776/2500_rerun/", seed)
savedir = file.path("/work/users/t/i/tianyou/two_sample/UKB_Yalin/RATS_776_2500/", seed)
if (!dir.exists(savedir)){
  dir.create(savedir)
}
if (!dir.exists(file.path(savedir, "perms"))){
  dir.create(file.path(savedir, "perms"))
}

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
# eqdist.etest(test, sizes = c(dim(AD_mat)[1], dim(CN_mat)[1]), R=101)

#set.seed(1024)
# pnull = lgamma.getnull(AD_mat, CN_mat, perm = 200)

X = AD_mat
Y = CN_mat
dat = rbind(X, Y)
dist_ed = l2.dist(X, Y)
saveRDS(dist_ed, file.path(savedir, "perms", "dist_ed.rds"))
dist_rbf = rbf(X, Y, distance="euclidean", median=TRUE)
saveRDS(dist_rbf, file.path(savedir, "perms", "dist_rbf.rds"))
dist = dist.array(X,Y)
saveRDS(dist, file.path(savedir, "perms", "dist_array.rds"))




