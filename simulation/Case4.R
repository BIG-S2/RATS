library(tidyverse)
library(MASS)
library(mvtnorm)
library(Ball)
library(lobstr)

args <- commandArgs(trailingOnly = TRUE)
signal_dim=as.numeric(args[1])
ri = as.numeric(args[2])
power_rep_n=as.numeric(args[3])
nullperm = as.numeric(args[4])

print(paste0("rep number: ", power_rep_n))

### distance function to calculate (X_ki - Y_li)^2
### when gamma = 1, this is equivalent to L1 modification when summed up
dist.array <- function(X, Y){
  Z = rbind(X,Y)
  if (dim(X)[2]!=dim(Y)[2]){
    print("Dimensions do not match!")
  }
  
  d = dim(Z)[2]
  l2_dist = array(0, dim = c(dim(Z)[1], dim(Z)[1], d))
  for (i in 1:d){
    l2_dist[,,i] = unname(as.matrix(dist(Z[,i])))
  }
  return(l2_dist)
}



##### further modification: take dimension-wise ED as input and no need to calculate it for every gamma
## Energy distance modified
energymod <- function(X, Y, gamma=1){
  n = dim(X)[1]
  m = dim(Y)[1]
  dist = dist.array(X, Y)
  
  dist_XX = dist[1:n,1:n,]
  dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  dist_XY = dist[1:n,(n+1):(n+m),]
  
  #### sum for ith dimensions
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    mean_XX = apply(dist_XX, 3, sum)/(n*(n-1))
    mean_YY = apply(dist_YY, 3, sum)/(m*(m-1))
    mean_XY = apply(dist_XY, 3, sum)/(m*n)
    energy.stat = max(2*mean_XY - mean_XX - mean_YY)
  } else {
    mean_XX = apply(dist_XX, 3, sum)/(n*(n-1))
    mean_YY = apply(dist_YY, 3, sum)/(m*(m-1))
    mean_XY = apply(dist_XY, 3, sum)/(m*n)
    energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
    #energy.stat = sum(((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)^(gamma))
  }
  
  return(energy.stat)
}


## energy statistic with matrix as input
energymod.mat <- function(dist, nx, ny, gamma=1){
  n = nx
  m = ny
  dist_XX = dist[1:n,1:n,]
  dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  dist_XY = dist[1:n,(n+1):(n+m),]
  
  #### sum for ith dimensions
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    mean_XX = apply(dist_XX, 3, sum)/(n*(n-1))
    mean_YY = apply(dist_YY, 3, sum)/(m*(m-1))
    mean_XY = apply(dist_XY, 3, sum)/(m*n)
    #energy.stat = max((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)
    energy.stat = max(2*mean_XY - mean_XX - mean_YY)
  } else {
    mean_XX = apply(dist_XX, 3, sum)/(n*(n-1))
    mean_YY = apply(dist_YY, 3, sum)/(m*(m-1))
    mean_XY = apply(dist_XY, 3, sum)/(m*n)
    #energy.stat = sum(((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)^(gamma))
    energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
  }
  
  return(energy.stat)
}


### energy distance test
energymod.test <- function(X, Y, gamma = 1, perm=101, lst=FALSE){
  edmod_stat = energymod(X, Y, gamma = gamma)
  
  dist = dist.array(X,Y)
  dat = rbind(X,Y)
  nx = dim(X)[1]
  ny = dim(Y)[1]
  energymod_mat_perm = numeric(perm)
  for (i in 1:perm){
    sample = sample(1:(nx+ny), (nx+ny))
    X_sample = dat[sample[1:nx],]
    Y_sample = dat[sample[(nx+1):(nx+ny)],]
    mat_sample = dist[sample, sample, ]
    energymod_mat_perm[i] <- energymod.mat(mat_sample, nx, ny, gamma = gamma)
  }
  pval_mat = mean(energymod_mat_perm >= edmod_stat)
  if (lst==FALSE) {
    return(pval_mat)
  }
  else {
    return(data.frame(pval_mat, ed_stat))
  }
}



##### ball divergence modified
## ball distance modified
bdmod <- function(X, Y, gamma=1){
  n = dim(X)[1]
  m = dim(Y)[1]
  p = dim(X)[2]
  dist = dist.array(X, Y)
  
  bdm = apply(dist, 3, bd, distance = TRUE, size = c(n,m))
  
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    bd.stat = max(bdm)
  } else {
    bd.stat = sum(bdm^(gamma))
    #bd.stat = sum(((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)^(gamma))
  }
  
  return(bd.stat)
}


## bd statistic with matrix as input
bdmod.mat <- function(dist, nx, ny, gamma=1){
  bdm = apply(dist, 3, bd, distance = TRUE, size = c(nx,ny))
  
  #### sum for ith dimensions
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    bd.stat = max(bdm)
  } else {
    bd.stat = sum(bdm^(gamma))
  }
  
  return(bd.stat)
}


### classical energy distances
##L1 distance without taking square root
l1.dist <- function(X, Y){
  Z = rbind(X,Y)
  d = dim(Z)[2]
  l1_dist = matrix(0, nrow = dim(Z)[1], ncol = dim(Z)[1])
  for (i in 1:dim(Z)[1]){
    for (j in 1:dim(Z)[1]){
      if (j < i) {l1_dist[i,j] = l1_dist[j,i]}
      else {
        l1_dist[i,j] = sum(abs(Z[i,]-Z[j,]))
      }
    }
  }
  return(l1_dist)
}

##L2 distance matrix
l2.dist <- function(X, Y){
  Z = rbind(X,Y)
  ZZ = Z %*% t(Z)
  diag_ZZ = diag(ZZ)
  diag_ZZ_F = matrix(diag_ZZ, nrow = dim(Z)[1], ncol = dim(Z)[1])
  diag_ZZ_F_T = matrix(diag_ZZ, nrow = dim(Z)[1], ncol=dim(Z)[1], byrow = TRUE)
  sq_dist = diag_ZZ_F + diag_ZZ_F_T - 2*ZZ
  dist = sqrt(abs(sq_dist))
  return(dist)
}


## Energy distance
energy <- function(X, Y, distance="euclidean"){
  if (distance == "euclidean") {
    dist = l2.dist(X, Y)
  } else if (distance == "l1") {
    dist = l1.dist(X, Y)
  } else if (distance == "l1_sqrt") {
    dist = sqrt(l1.dist(X, Y))
  } else if (distance == "l3") {
    dist = l3.dist(X, Y)^(1/3)
  } else if (distance == "l4") {
    dist = l4.dist(X, Y)^(1/4)
  } else if (distance == "lmax") {
    dist = sqrt(lmax.dist(X, Y))
  }
  
  n = dim(X)[1]
  m = dim(Y)[1]
  energy_XX = dist[1:n,1:n]
  energy_YY = dist[(n+1):(n+m),(n+1):(n+m)]
  energy_XY = dist[1:n,(n+1):(n+m)]
  mean_XX = (sum(energy_XX))/(n*(n-1))  #the -n is for minusing the i=j terms
  mean_YY = (sum(energy_YY))/(m*(m-1))
  mean_XY = (sum(energy_XY))/(m*n)
  #return((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)
  return(2*mean_XY - mean_XX - mean_YY)
}

## energy statistic with matrix as input
energy.mat <- function(mat, nx, ny){
  n = nx
  m = ny
  energy_XX = mat[1:n,1:n]
  energy_YY = mat[(n+1):(n+m),(n+1):(n+m)]
  energy_XY = mat[1:n,(n+1):(n+m)]
  mean_XX = (sum(energy_XX))/(n*(n-1))  #the -n is for minusing the i=j terms
  mean_YY = (sum(energy_YY))/(m*(m-1))
  mean_XY = (sum(energy_XY))/(m*n)
  #return((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)
  return(2*mean_XY - mean_XX - mean_YY)
}

### energy distance test
energy.test <- function(X, Y, distance="euclidean", perm=101, lst=FALSE){
  ed_stat = energy(X, Y, distance = distance)
  if (distance == "euclidean") {
    dist = l2.dist(X, Y)
  } else if (distance == "l1") {
    dist = l1.dist(X, Y)
  } else if (distance == "l1_sqrt") {
    dist = sqrt(l1.dist(X, Y))
  } else if (distance == "l3") {
    dist = l3.dist(X, Y)^(1/3)
  } else if (distance == "l4") {
    dist = l4.dist(X, Y)^(1/4)
  } else if (distance == "lmax") {
    dist = sqrt(lmax.dist(X, Y))
  }
  
  dat = rbind(X,Y)
  nx = dim(X)[1]
  ny = dim(Y)[1]
  energy_mat_perm = numeric(perm)
  for (i in 1:perm){
    sample = sample(1:(nx+ny), (nx+ny))
    X_sample = dat[sample[1:nx],]
    Y_sample = dat[sample[(nx+1):(nx+ny)],]
    mat_sample = dist[sample, sample]
    energy_mat_perm[i] <- energy.mat(mat_sample, nx, ny)
  }
  pval_mat = mean(energy_mat_perm >= ed_stat)
  if (lst==FALSE) {
    return(pval_mat)
  }
  else {
    return(data.frame(pval_mat, ed_stat))
  }
}



### function to calculate all statistics for a permuted sample
minP.cal = function(mat, nx, ny){
  value_edl1_perm <- energymod.mat(mat, nx, ny, gamma = 1)
  value_edl2_perm <- energymod.mat(mat, nx, ny, gamma = 2)
  value_edl3_perm <- energymod.mat(mat, nx, ny, gamma = 3)
  value_edl4_perm <- energymod.mat(mat, nx, ny, gamma = 4)
  value_edl5_perm <- energymod.mat(mat, nx, ny, gamma = 5)
  value_edl6_perm <- energymod.mat(mat, nx, ny, gamma = 6)
  value_edlmax_perm <- energymod.mat(mat, nx, ny, gamma = "max")
  
  value_bdl1_perm <- bdmod.mat(mat, nx, ny, gamma = 1)
  value_bdl2_perm <- bdmod.mat(mat, nx, ny, gamma = 2)
  value_bdl3_perm <- bdmod.mat(mat, nx, ny, gamma = 3)
  value_bdl4_perm <- bdmod.mat(mat, nx, ny, gamma = 4)
  value_bdl5_perm <- bdmod.mat(mat, nx, ny, gamma = 5)
  value_bdl6_perm <- bdmod.mat(mat, nx, ny, gamma = 6)
  value_bdlmax_perm <- bdmod.mat(mat, nx, ny, gamma = "max")
  
  value = c(value_edl1_perm, value_edl2_perm, value_edl3_perm, value_edl4_perm,
            value_edl5_perm, value_edl6_perm, value_edlmax_perm, value_bdl1_perm,
            value_bdl2_perm, value_bdl3_perm, value_bdl4_perm, value_bdl5_perm,
            value_bdl6_perm, value_bdlmax_perm)
  return(value)
}




############ 
#### function to calculate the min p test
minP.getnull = function(X, Y, perm = 100){
  nx = dim(X)[1]
  ny = dim(Y)[1]
  dat = rbind(X, Y)
  dist = dist.array(X,Y)
  dist_ed = l2.dist(X, Y)
  
  ## get the null distribution
  value_ed_perm = numeric(perm)
  value_bd_perm = numeric(perm)
  
  value_perm = matrix(NA, nrow = perm, ncol = 14)
  
  set.seed(NULL)
  for (i in 1:perm){
    if (i %% 10 == 0) {print(paste0("perm: ", i)); print(mem_used())}
    sample = sample(1:(nx+ny), (nx+ny))
    X_sample = dat[sample[1:nx],]
    Y_sample = dat[sample[(nx+1):(nx+ny)],]
    mat_sample = dist[sample, sample, ]
    value_ed_perm[i] <- energy.mat(dist_ed[sample, sample], nx, ny)
    value_perm[i,] = minP.cal(mat_sample, nx, ny)
    value_bd_perm[i] <- bd(dist_ed[sample, sample], distance = TRUE, size = c(nx,ny))
  }
  
  value_perm = as.data.frame(value_perm)
  colnames(value_perm) = c("ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX", 
                           "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX")
  value_perm$ED = value_ed_perm
  value_perm$BD = value_bd_perm
  return(value_perm)
}



##########
### high dimensional sparse signal setting: normal location shift
### case 14 in summary
##########
n = m = 50
dim = 100
r= c(0, 0.175, 0.35, 0.525, 0.7, 0.875, 1.05, 1.225, 1.4)
mu_X = r[ri]
print(paste0("ss: ", n))
print(paste0("dimension: ", dim))
print(paste0("mu: ", mu_X))
print(paste0("nullperm:", nullperm))

### get the sample
set.seed(power_rep_n **2 + 3*power_rep_n)

if (!file.exists(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))){
  X = mvrnorm(n, rep(0, dim), diag(dim))
  X[,1:signal_dim] = mvrnorm(n, rep(mu_X, signal_dim), diag(signal_dim))
  Y = mvrnorm(m, rep(0, dim), diag(dim))
  save(X, file=paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  save(Y, file=paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
} else {
  load(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  load(paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
}

pnull = minP.getnull(X, Y, perm = 100)
write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull_", nullperm, ".csv"))

###########





##########
### high dimensional sparse signal setting: normal scale shift
### case 15 in summary
##########
n = m = 50
dim = 100
r = c(1, 1.75, 2.5, 3.25, 4, 4.75, 5.5, 6.25, 7)
var_X = r[ri]
print(paste0("ss: ", n))
print(paste0("dimension: ", dim))
print(paste0("var: ", var_X))
print(paste0("nullperm:", nullperm))

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

pnull = minP.getnull(X, Y, perm = 100)
write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull_", nullperm, ".csv"))

###########


##########
### high dimensional sparse signal setting: normal and poisson
### case 16 in summary
##########
n = m = 50
dim = 100
r = c(1,2,5,10,15,20,25,30)
dim_X = r[ri]
print(paste0("ss: ", n))
print(paste0("dimension: ", dim))
print(paste0("signal_dim: ", dim_X))
print(paste0("nullperm:", nullperm))

### get the sample
set.seed(power_rep_n **2 + 3*power_rep_n)

if (!file.exists(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))){
  X = mvrnorm(n, rep(1, dim), diag(dim))
  X[,1:dim_X] = matrix(rpois(n * dim_X, 1), ncol = dim_X, nrow = n)
  Y = mvrnorm(m, rep(1, dim), diag(dim))
  save(X, file=paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  save(Y, file=paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
} else {
  load(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  load(paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
}

pnull = minP.getnull(X, Y, perm = 100)
write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull_", nullperm, ".csv"))

###########




##########
### high dimensional sparse signal setting: normal and t
### case 18 in summary
##########
n = m = 50
dim = 100
r = c(1,10,20,30, 40,50, 60, 70)
dim_X = r[ri]
print(paste0("ss: ", n))
print(paste0("dimension: ", dim))
print(paste0("signal_dim: ", dim_X))
print(paste0("nullperm:", nullperm))

### get the sample
set.seed(power_rep_n **2 + 3*power_rep_n)

if (!file.exists(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))){
  X = mvrnorm(n, rep(0, dim), 5 * diag(dim))
  X[,1:dim_X] = as.matrix(rmvt(n, sigma = 3 * diag(dim_X), df = 5))
  Y = mvrnorm(m, rep(0, dim), 5 * diag(dim))
  save(X, file=paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  save(Y, file=paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
} else {
  load(paste0("./r", ri, "/n", power_rep_n, "/X_", power_rep_n, ".RData"))
  load(paste0("./r", ri, "/n", power_rep_n, "/Y_", power_rep_n, ".RData"))
}

pnull = minP.getnull(X, Y, perm = 100)
write_csv(pnull, paste0("./r", ri, "/n", power_rep_n, "/pnull_", nullperm, ".csv"))

###########

