##### ball divergence modified
## ball divergence modified
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


## bd statistic with matrix as input, for a sequence of gamma values
bdmod.mat.seq <- function(dist, nx, ny, gamma_seq){
  bdm = apply(dist, 3, bd, distance = TRUE, size = c(nx,ny))
  
  #### sum for ith dimensions
  BD_seq = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    gamma = ifelse(tolower(gamma_seq[i]) == "max", "max", as.numeric(gamma_seq[i]))
    if (!is.numeric(gamma) & gamma != "max") {
      print("Gamma must be numeric!")
    } else if (gamma=="max"){
      bd.stat = max(bdm)
    } else {
      bd.stat = sum(bdm^(gamma))
    }
    BD_seq[i] = bd.stat
  }
  names(BD_seq) = paste0("BD", gamma_seq)
  return(BD_seq)
}


## ball divergence modified, for a sequence of gamma values
bdmod.seq <- function(X, Y, gamma_seq){
  n = dim(X)[1]
  m = dim(Y)[1]
  p = dim(X)[2]
  dist = dist.array(X, Y)
  
  bdm = apply(dist, 3, bd, distance = TRUE, size = c(n,m))
  
  BD_seq = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    gamma = ifelse(tolower(gamma_seq[i]) == "max", "max", as.numeric(gamma_seq[i]))
    if (!is.numeric(gamma) & gamma != "max") {
      print("Gamma must be numeric!")
    } else if (gamma=="max"){
      bd.stat = max(bdm)
    } else {
      bd.stat = sum(bdm^(gamma))
    }
    BD_seq[i] = bd.stat
  }
  names(BD_seq) = paste0("BD", gamma_seq)
  return(BD_seq)
}

