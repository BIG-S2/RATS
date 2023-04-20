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
perm = as.numeric(args[2])
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
# eqdist.etest(test, sizes = c(dim(AD_mat)[1], dim(CN_mat)[1]), R=101)

#set.seed(1024)
# pnull = lgamma.getnull(AD_mat, CN_mat, perm = 200)

X = AD_mat
Y = CN_mat
dat = rbind(X, Y)



# nx = dim(X)[1]
# ny = dim(Y)[1]
# dat = rbind(X, Y)
# dist = readRDS(file.path(savedir, "dist_array.rds"))
# dist_ed = readRDS(file.path(savedir, "dist_ed.rds"))
# dist_rbf = readRDS(file.path(savedir, "dist_rbf.rds"))
# 
# set.seed(NULL)
# 
# sample = sample(1:(nx+ny), (nx+ny))
# X_sample = dat[sample[1:nx],]
# Y_sample = dat[sample[(nx+1):(nx+ny)],]
# mat_sample = dist[sample, sample, ]
# value_mmd_perm <- mmd.mat(dist_rbf[sample, sample], nx, ny)
# value_ed_perm <- energy.mat(dist_ed[sample, sample], nx, ny)
# value_perm = lgamma.cal(mat_sample, nx, ny)
# value_bd_perm <- bd(dist_ed[sample, sample], distance = TRUE, size = c(nx,ny))
# 
# names(value_perm) = c("ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX", 
#                          "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX")
# value_perm = as.data.frame(t(value_perm))
# value_perm$MMD = value_mmd_perm
# value_perm$ED = value_ed_perm
# value_perm$BD = value_bd_perm

dist = readRDS(file.path(savedir, "dist_array.rds"))
dist_ed = readRDS(file.path(savedir, "dist_ed.rds"))
# dist_rbf = readRDS(file.path(savedir, "dist_rbf.rds"))

gamma_seq = c(1:6, "MAX")

stats_null = lgamma.getnull.parallel(X=X, Y=Y, dist=dist, dist_ed=dist_ed, gamma_seq=gamma_seq)
write_csv(stats_null, file.path(savedir, paste0("null_stats_perm", perm, ".csv")))
