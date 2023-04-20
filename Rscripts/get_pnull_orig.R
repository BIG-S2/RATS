source("cal_ED.R")
source("cal_BD.R")
source("cal_MMD.R")


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





### function to calculate all statistics for a permuted sample
lgamma.cal = function(mat, nx, ny){
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
lgamma.getnull = function(X, Y, perm = 100){
  nx = dim(X)[1]
  ny = dim(Y)[1]
  dat = rbind(X, Y)
  dist = dist.array(X,Y)
  dist_ed = l2.dist(X, Y)
  dist_rbf = rbf(X, Y, distance="euclidean", median=TRUE)
  
  ## get the null distribution
  value_mmd_perm = numeric(perm)
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
    value_mmd_perm[i] <- mmd.mat(dist_rbf[sample, sample], nx, ny)
    value_ed_perm[i] <- energy.mat(dist_ed[sample, sample], nx, ny)
    value_perm[i,] = lgamma.cal(mat_sample, nx, ny)
    value_bd_perm[i] <- bd(dist_ed[sample, sample], distance = TRUE, size = c(nx,ny))
  }
  
  value_perm = as.data.frame(value_perm)
  colnames(value_perm) = c("ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX", 
                           "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX")
  value_perm$MMD = value_mmd_perm
  value_perm$ED = value_ed_perm
  value_perm$BD = value_bd_perm
  return(value_perm)
}


###############
#### generate AR(1) covariance matrix
ar1_cor <- function(size, rho) {
  exponent <- abs(matrix(1:size - 1, nrow = size, ncol = size, byrow = TRUE) - 
                    (1:size - 1))
  return(rho^exponent)
}
###############

