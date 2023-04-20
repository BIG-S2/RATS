## Energy distance modified
energymod <- function(X, Y, gamma=1){
  n = dim(X)[1]
  m = dim(Y)[1]
  dist = dist.array(X, Y)
  
  # dist_XX = dist[1:n,1:n,]
  # dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  # dist_XY = dist[1:n,(n+1):(n+m),]
  
  #### sum for ith dimensions
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
    mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
    mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
    energy.stat = max(2*mean_XY - mean_XX - mean_YY)
  } else {
    mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
    mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
    mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
    energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
    #energy.stat = sum(((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)^(gamma))
  }
  
  return(energy.stat)
}


## energy statistic with matrix as input
energymod.mat <- function(dist, nx, ny, gamma=1){
  n = nx
  m = ny
  # dist_XX = dist[1:n,1:n,]
  # dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  # dist_XY = dist[1:n,(n+1):(n+m),]
  
  #### sum for ith dimensions
  if (!is.numeric(gamma) & gamma != "max") {
    print("Gamma must be numeric!")
  } else if (gamma=="max"){
    mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
    mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
    mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
    #energy.stat = max((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)
    energy.stat = max(2*mean_XY - mean_XX - mean_YY)
  } else {
    mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
    mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
    mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
    #energy.stat = sum(((mean_XY - mean_XX)^2 + (mean_XY - mean_YY)^2)^(gamma))
    energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
  }
  
  return(energy.stat)
}


## energy statistic with matrix as input, for a sequence of gamma values
energymod.mat.seq <- function(dist, nx, ny, gamma_seq){
  n = nx
  m = ny
  # dist_XX = dist[1:n,1:n,]
  # dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  # dist_XY = dist[1:n,(n+1):(n+m),]
  mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
  mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
  mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
  
  #### sum for ith dimensions
  ED_seq = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    gamma = ifelse(tolower(gamma_seq[i]) == "max", "max", as.numeric(gamma_seq[i]))
    if (!is.numeric(gamma) & gamma != "max") {
      print("Gamma must be numeric!")
    } else if (gamma=="max"){
      energy.stat = max(2*mean_XY - mean_XX - mean_YY)
    } else {
      energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
    }
    ED_seq[i] = energy.stat
  }
  names(ED_seq) = paste0("ED", gamma_seq)
  return(ED_seq)
}

## energy distance modified, for a sequence of gamma values
energymod.seq <- function(X, Y, gamma_seq){
  n = dim(X)[1]
  m = dim(Y)[1]
  dist = dist.array(X, Y)
  # dist_XX = dist[1:n,1:n,]
  # dist_YY = dist[(n+1):(n+m),(n+1):(n+m),]
  # dist_XY = dist[1:n,(n+1):(n+m),]
  mean_XX = apply(dist[1:n,1:n,], 3, sum)/(n*(n-1))
  mean_YY = apply(dist[(n+1):(n+m),(n+1):(n+m),], 3, sum)/(m*(m-1))
  mean_XY = apply(dist[1:n,(n+1):(n+m),], 3, sum)/(m*n)
  
  #### sum for ith dimensions
  ED_seq = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    gamma = ifelse(tolower(gamma_seq[i]) == "max", "max", as.numeric(gamma_seq[i]))
    if (!is.numeric(gamma) & gamma != "max") {
      print("Gamma must be numeric!")
    } else if (gamma=="max"){
      energy.stat = max(2*mean_XY - mean_XX - mean_YY)
    } else {
      energy.stat = sum((2*mean_XY - mean_XX - mean_YY)^(gamma))
    }
    ED_seq[i] = energy.stat
  }
  names(ED_seq) = paste0("ED", gamma_seq)
  return(ED_seq)
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


