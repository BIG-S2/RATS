############ 
#### function to calculate the min p test
lgamma.getnull.rot.from.list = function(X, Y, rotation_list, perm = 100){
  nx = dim(X)[1]
  ny = dim(Y)[1]
  
  #dist = dist.array(X,Y)
  #dist_ed = l2.dist(X, Y)
  #dist_rbf = rbf(X, Y, distance="euclidean", median=TRUE)
  
  ## get the null distribution
  value_mmd_perm = numeric(perm)
  value_ed_perm = numeric(perm)
  value_bd_perm = numeric(perm)
  
  value_perm = matrix(NA, nrow = perm, ncol = 14)
  
  set.seed(NULL)
  for (i in 1:perm){
    if (i %% 5 == 0) {print(paste0("perm: ", i)); print(mem_used())}

    dat = rbind(X, Y)
    sample = sample(1:(nx+ny), (nx+ny))
    new_X = dat[sample[1:nx],]
    new_Y = dat[sample[(nx+1):(nx+ny)],]
    
    for (g in 1:6){
      max_rotated = find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                                basefunc = "ED", gamma = g)
      value_perm[i,g] = max_rotated
    }
    max_rotated = find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                              basefunc = "ED", gamma = "max")
    value_perm[i,7] = max_rotated
    for (g in 8:13){
      max_rotated = find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                                basefunc = "BD", gamma = g-7)
      value_perm[i,g] = max_rotated
    }
    max_rotated = find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                              basefunc = "BD", gamma = "max")
    value_perm[i,14] = max_rotated

    #value_mmd_perm[i] <- mmd(X_sample, Y_sample, distance="euclidean", median = TRUE)
    value_ed_perm[i] <- find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                                    basefunc = "EDorig")
    value_bd_perm[i] <- find.max.rotation.from.list(new_X, new_Y, rotation_list = rotation_list,
                                                    basefunc = "BDorig")
  }
  
  value_perm = as.data.frame(value_perm)
  colnames(value_perm) = c("ED1", "ED2", "ED3", "ED4", "ED5", "ED6", "EDMAX", 
                           "BD1", "BD2", "BD3", "BD4", "BD5", "BD6", "BDMAX")
  #value_perm$MMD = value_mmd_perm
  value_perm$ED = value_ed_perm
  value_perm$BD = value_bd_perm
  return(value_perm)
}



## calculate the max rotation for given rotation list and test type
find.max.rotation.from.list = function(X, Y, rotation_list, basefunc, gamma=NULL){
  best.X = X
  best.Y = Y
  if (basefunc == "ED"){
    best.cur = energymod(X, Y, gamma = gamma)
  } else if (basefunc == "BD"){
    best.cur = bdmod(X, Y, gamma = gamma)
  } else if (basefunc == "EDorig"){
    best.cur = energy(X, Y)
  } else if (basefunc == "BDorig"){
    best.cur = bd(X, Y, distance = FALSE)
  }
  for (rot in 1:length(rotation_list)){
    newX = X %*% rotation_list[[rot]]
    newY = Y %*% rotation_list[[rot]]
    if (basefunc == "ED"){
      newdist = energymod(newX, newY, gamma=gamma)
    } else if (basefunc == "BD"){
      newdist = bdmod(newX, newY, gamma=gamma)
    } else if (basefunc == "EDorig"){
      newdist = energy(newX, newY)
    } else if (basefunc == "BDorig"){
      newdist = bd(newX, newY, distance = FALSE)
    }
    if (newdist > best.cur){
      best.cur = newdist
      best.X = newX
      best.Y = newY
    }
  }
  return(best.cur)
}

