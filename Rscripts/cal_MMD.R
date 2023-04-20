##rbf kernel
## X,Y are n/m-by-p matrix
rbf <- function(X, Y, sigma=1, distance="euclidean", median = TRUE){
  if (distance == "euclidean") {
    dist = l2.dist(X,Y)
  } else if (distance == "l1") {
    dist = l1.dist(X, Y)
  }
  sq_dist = dist ^ 2
  if (median==TRUE) {
    sigma = median(dist)
  }
  rbf_mat = exp(-dist^2/(2*sigma^2))
  return(rbf_mat)
}


## MMD
mmd <- function(X, Y, sigma=1, distance="euclidean", median = TRUE){
  rbf_mat = rbf(X=X, Y=Y, sigma=sigma, distance=distance, median = median)
  n = dim(X)[1]
  m = dim(Y)[1]
  rbf_XX = rbf_mat[1:n,1:n]
  rbf_YY = rbf_mat[(n+1):(n+m),(n+1):(n+m)]
  rbf_XY = rbf_mat[1:n,(n+1):(n+m)]
  mean_XX = (sum(rbf_XX)-n)/(n*(n-1))  #the -n is for minusing the i=j terms
  mean_YY = (sum(rbf_YY)-m)/(m*(m-1))
  mean_XY = (sum(rbf_XY))/(m*n)
  return(mean_XX + mean_YY - 2*mean_XY)
}

## MMD with rbf matrix as input
mmd.mat <- function(rbf_mat, nx, ny){
  n = nx
  m = ny
  rbf_XX = rbf_mat[1:n,1:n]
  rbf_YY = rbf_mat[(n+1):(n+m),(n+1):(n+m)]
  rbf_XY = rbf_mat[1:n,(n+1):(n+m)]
  mean_XX = (sum(rbf_XX)-n)/(n*(n-1))
  mean_YY = (sum(rbf_YY)-m)/(m*(m-1))
  mean_XY = (sum(rbf_XY))/(m*n)
  return(mean_XX + mean_YY - 2*mean_XY)
}


## MMD TEST
mmd.test <- function(X, Y, sigma=1, distance="euclidean", perm=101, lst=FALSE, median = TRUE){
  mmd_stat = mmd(X, Y, sigma = sigma, distance = distance, median=median)
  rbf_mat = rbf(X=X, Y=Y, sigma=sigma, distance=distance, median=median)
  dat = rbind(X,Y)
  nx = dim(X)[1]
  ny = dim(Y)[1]
  mmd_mat_perm = numeric(perm)
  for (i in 1:perm){
    sample = sample(1:(nx+ny), (nx+ny))
    X_sample = dat[sample[1:nx],]
    Y_sample = dat[sample[(nx+1):(nx+ny)],]
    rbf_mat_sample = rbf_mat[sample, sample]
    mmd_mat_perm[i] <- mmd.mat(rbf_mat_sample, nx, ny)
  }
  pval_mat = mean(mmd_mat_perm > mmd_stat)
  if (lst==FALSE) {
    return(pval_mat)
  }
  else {
    return(data.frame(pval_mat, mmd_stat))
  }
}


