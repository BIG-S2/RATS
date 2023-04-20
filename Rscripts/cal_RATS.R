cal.RATS.stat = function(X, Y, pnull){
  #### sample test statistic and p values
  ########
  nperm = dim(pnull)[1]
  
  # val_mmd = mmd(X, Y, distance = "euclidean", median = TRUE)
  val_ed = energy(X, Y)
  val_ed1 = energymod(X, Y, gamma=1)
  val_ed2 = energymod(X, Y, gamma=2)
  val_ed3 = energymod(X, Y, gamma=3)
  val_ed4 = energymod(X, Y, gamma=4)
  val_ed5 = energymod(X, Y, gamma=5)
  val_ed6 = energymod(X, Y, gamma=6)
  val_edmax = energymod(X, Y, gamma="max")
  
  val_bd = bd(X, Y)
  val_bd1 = bdmod(X, Y, gamma=1)
  val_bd2 = bdmod(X, Y, gamma=2)
  val_bd3 = bdmod(X, Y, gamma=3)
  val_bd4 = bdmod(X, Y, gamma=4)
  val_bd5 = bdmod(X, Y, gamma=5)
  val_bd6 = bdmod(X, Y, gamma=6)
  val_bdmax = bdmod(X, Y, gamma="max")
  
  # pval_mmd = (sum(pnull$MMD >= val_mmd) + 1) / (nperm + 1)
  pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
  pval_ed1 = (sum(pnull$ED1 >= val_ed1) + 1) / (nperm + 1)
  pval_ed2 = (sum(pnull$ED2 >= val_ed2) + 1) / (nperm + 1)
  pval_ed3 = (sum(pnull$ED3 >= val_ed3) + 1) / (nperm + 1)
  pval_ed4 = (sum(pnull$ED4 >= val_ed4) + 1) / (nperm + 1)
  pval_ed5 = (sum(pnull$ED5 >= val_ed5) + 1) / (nperm + 1)
  pval_ed6 = (sum(pnull$ED6 >= val_ed6) + 1) / (nperm + 1)
  pval_edmax = (sum(pnull$EDMAX >= val_edmax) + 1) / (nperm + 1)
  minpval_ed = min(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax)
  minpval_ed_simp = min(pval_ed, pval_ed1, pval_ed2, pval_ed4)
  
  pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
  pval_bd1 = (sum(pnull$BD1 >= val_bd1) + 1) / (nperm + 1)
  pval_bd2 = (sum(pnull$BD2 >= val_bd2) + 1) / (nperm + 1)
  pval_bd3 = (sum(pnull$BD3 >= val_bd3) + 1) / (nperm + 1)
  pval_bd4 = (sum(pnull$BD4 >= val_bd4) + 1) / (nperm + 1)
  pval_bd5 = (sum(pnull$BD5 >= val_bd5) + 1) / (nperm + 1)
  pval_bd6 = (sum(pnull$BD6 >= val_bd6) + 1) / (nperm + 1)
  pval_bdmax = (sum(pnull$BDMAX >= val_bdmax) + 1) / (nperm + 1)
  minpval_bd = min(pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax)
  minpval_bd_simp = min(pval_bd, pval_bd1, pval_bd2, pval_bd4)
  minpval = min(minpval_ed, minpval_bd)
  minpval_simp = min(minpval_ed_simp, minpval_bd_simp)
  return(list(c(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax,
                pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax),
              c(minpval_ed, minpval_bd, minpval),
              c(minpval_ed_simp, minpval_bd_simp, minpval_simp)))
}


#### input a sequence of gamma values
cal.RATS.stat.seq = function(X, Y, pnull, gamma_seq){
  #### sample test statistic and p values
  ########
  nperm = dim(pnull)[1]
  
  # val_mmd = mmd(X, Y, distance = "euclidean", median = TRUE)
  val_ed = energy(X, Y)
  val_ed_mod = energymod.seq(X, Y, gamma_seq)
  
  val_bd = bd(X, Y)
  val_bd_mod = bdmod.seq(X, Y, gamma_seq)
  
  # pval_mmd = (sum(pnull$MMD >= val_mmd) + 1) / (nperm + 1)
  pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
  pval_ed_mod = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    pval_ed_mod[i] = (sum(pnull[[paste0("ED", gamma_seq[i])]] >= val_ed_mod[paste0("ED", gamma_seq[i])]) + 1) / (nperm + 1)
  }
  pval_ed_mod = c(pval_ed, pval_ed_mod)
  names(pval_ed_mod) = c("ED", paste0("ED", gamma_seq))
  minpval_ed = min(pval_ed_mod)
  # if (!is.null(gamma_seq_simp)){
  #   minpval_ed_simp = min(pval_ed, pval_ed_mod[gamma_seq_simp])  ## gamma_seq_simp corresponds to the order in gamma_seq that is included
  # }
  
  pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
  pval_bd_mod = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    pval_bd_mod[i] = (sum(pnull[[paste0("BD", gamma_seq[i])]] >= val_bd_mod[paste0("BD", gamma_seq[i])]) + 1) / (nperm + 1)
  }
  pval_bd_mod = c(pval_bd, pval_bd_mod)
  names(pval_bd_mod) = c("BD", paste0("BD", gamma_seq))
  minpval_bd = min(pval_bd_mod)
  # if (!is.null(gamma_seq_simp)){
  #   minpval_bd_simp = min(pval_bd, pval_bd_mod[gamma_seq_simp])
  # }
  minpval = min(minpval_ed, minpval_bd)
  # minpval_simp = min(minpval_ed_simp, minpval_bd_simp)
  return(list(c(pval_ed_mod,
                pval_bd_mod),
              c(minpval_ed, minpval_bd, minpval)))
}


cal.RATS.stat.mat = function(X, Y, dist, dist_l2, pnull){
  #### sample test statistic and p values
  ########
  nperm = dim(pnull)[1]
  nx = dim(X)[1]
  ny = dim(Y)[1]
  
  # val_mmd = mmd(X, Y, distance = "euclidean", median = TRUE)
  # val_ed = energy(X, Y)
  val_ed = energy.mat(dist_l2, nx, ny)
  # val_ed1 = energymod(X, Y, gamma=1)
  # val_ed2 = energymod(X, Y, gamma=2)
  # val_ed3 = energymod(X, Y, gamma=3)
  # val_ed4 = energymod(X, Y, gamma=4)
  # val_ed5 = energymod(X, Y, gamma=5)
  # val_ed6 = energymod(X, Y, gamma=6)
  # val_edmax = energymod(X, Y, gamma="max")
  val_ed1 = energymod.mat(dist, nx, ny, gamma=1)
  val_ed2 = energymod.mat(dist, nx, ny, gamma=2)
  val_ed3 = energymod.mat(dist, nx, ny, gamma=3)
  val_ed4 = energymod.mat(dist, nx, ny, gamma=4)
  val_ed5 = energymod.mat(dist, nx, ny, gamma=5)
  val_ed6 = energymod.mat(dist, nx, ny, gamma=6)
  val_edmax = energymod.mat(dist, nx, ny, gamma="max")
  
  val_bd = bd(X, Y)
  # val_bd1 = bdmod(X, Y, gamma=1)
  # val_bd2 = bdmod(X, Y, gamma=2)
  # val_bd3 = bdmod(X, Y, gamma=3)
  # val_bd4 = bdmod(X, Y, gamma=4)
  # val_bd5 = bdmod(X, Y, gamma=5)
  # val_bd6 = bdmod(X, Y, gamma=6)
  # val_bdmax = bdmod(X, Y, gamma="max")
  val_bd1 = bdmod.mat(dist, nx, ny, gamma=1)
  val_bd2 = bdmod.mat(dist, nx, ny, gamma=2)
  val_bd3 = bdmod.mat(dist, nx, ny, gamma=3)
  val_bd4 = bdmod.mat(dist, nx, ny, gamma=4)
  val_bd5 = bdmod.mat(dist, nx, ny, gamma=5)
  val_bd6 = bdmod.mat(dist, nx, ny, gamma=6)
  val_bdmax = bdmod.mat(dist, nx, ny, gamma="max")
  
  # pval_mmd = (sum(pnull$MMD >= val_mmd) + 1) / (nperm + 1)
  pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
  pval_ed1 = (sum(pnull$ED1 >= val_ed1) + 1) / (nperm + 1)
  pval_ed2 = (sum(pnull$ED2 >= val_ed2) + 1) / (nperm + 1)
  pval_ed3 = (sum(pnull$ED3 >= val_ed3) + 1) / (nperm + 1)
  pval_ed4 = (sum(pnull$ED4 >= val_ed4) + 1) / (nperm + 1)
  pval_ed5 = (sum(pnull$ED5 >= val_ed5) + 1) / (nperm + 1)
  pval_ed6 = (sum(pnull$ED6 >= val_ed6) + 1) / (nperm + 1)
  pval_edmax = (sum(pnull$EDMAX >= val_edmax) + 1) / (nperm + 1)
  minpval_ed = min(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax)
  minpval_ed_simp = min(pval_ed, pval_ed1, pval_ed2, pval_ed4)
  
  pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
  pval_bd1 = (sum(pnull$BD1 >= val_bd1) + 1) / (nperm + 1)
  pval_bd2 = (sum(pnull$BD2 >= val_bd2) + 1) / (nperm + 1)
  pval_bd3 = (sum(pnull$BD3 >= val_bd3) + 1) / (nperm + 1)
  pval_bd4 = (sum(pnull$BD4 >= val_bd4) + 1) / (nperm + 1)
  pval_bd5 = (sum(pnull$BD5 >= val_bd5) + 1) / (nperm + 1)
  pval_bd6 = (sum(pnull$BD6 >= val_bd6) + 1) / (nperm + 1)
  pval_bdmax = (sum(pnull$BDMAX >= val_bdmax) + 1) / (nperm + 1)
  minpval_bd = min(pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax)
  minpval_bd_simp = min(pval_bd, pval_bd1, pval_bd2, pval_bd4)
  minpval = min(minpval_ed, minpval_bd)
  minpval_simp = min(minpval_ed_simp, minpval_bd_simp)
  return(list(c(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax,
                pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax),
              c(minpval_ed, minpval_bd, minpval),
              c(minpval_ed_simp, minpval_bd_simp, minpval_simp)))
}


#### input a sequence of gamma values
cal.RATS.stat.mat.seq = function(X, Y, dist, dist_l2, pnull, gamma_seq){
  #### sample test statistic and p values
  ########
  nperm = dim(pnull)[1]
  nx = dim(X)[1]
  ny = dim(Y)[1]
  
  # val_mmd = mmd(X, Y, distance = "euclidean", median = TRUE)
  # val_ed = energy(X, Y)
  val_ed = energy.mat(dist_l2, nx, ny)
  val_ed_mod = energymod.mat.seq(dist, nx, ny, gamma_seq)
  
  val_bd = bd(X, Y)
  val_bd_mod = bdmod.mat.seq(dist, nx, ny, gamma_seq)
  
  # pval_mmd = (sum(pnull$MMD >= val_mmd) + 1) / (nperm + 1)
  pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
  pval_ed_mod = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    pval_ed_mod[i] = (sum(pnull[[paste0("ED", gamma_seq[i])]] >= val_ed_mod[paste0("ED", gamma_seq[i])]) + 1) / (nperm + 1)
  }
  pval_ed_mod = c(pval_ed, pval_ed_mod)
  names(pval_ed_mod) = c("ED", paste0("ED", gamma_seq))
  minpval_ed = min(pval_ed_mod)
  # if (!is.null(gamma_seq_simp)){
  #   minpval_ed_simp = min(pval_ed, pval_ed_mod[gamma_seq_simp])  ## gamma_seq_simp corresponds to the order in gamma_seq that is included
  # }
  
  pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
  pval_bd_mod = numeric(length(gamma_seq))
  for (i in 1:length(gamma_seq)){
    pval_bd_mod[i] = (sum(pnull[[paste0("BD", gamma_seq[i])]] >= val_bd_mod[paste0("BD", gamma_seq[i])]) + 1) / (nperm + 1)
  }
  pval_bd_mod = c(pval_bd, pval_bd_mod)
  names(pval_bd_mod) = c("BD", paste0("BD", gamma_seq))
  minpval_bd = min(pval_bd_mod)
  # if (!is.null(gamma_seq_simp)){
  #   minpval_bd_simp = min(pval_bd, pval_bd_mod[gamma_seq_simp])
  # }
  minpval = min(minpval_ed, minpval_bd)
  # minpval_simp = min(minpval_ed_simp, minpval_bd_simp)
  return(list(c(pval_ed_mod,
                pval_bd_mod),
              c(minpval_ed, minpval_bd, minpval)))
}



cal.RATS.stat.rot = function(X, Y, rotation_list, pnull){
  #### sample test statistic and p values
  ########
  nperm = dim(pnull)[1]
  
  val_ed = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                       basefunc = "EDorig")
  val_ed1 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 1)
  val_ed2 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 2)
  val_ed3 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 3)
  val_ed4 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 4)
  val_ed5 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 5)
  val_ed6 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "ED", gamma = 6)
  val_edmax = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                          basefunc = "ED", gamma = "max")
  
  val_bd = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                       basefunc = "BDorig")
  val_bd1 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 1)
  val_bd2 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 2)
  val_bd3 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 3)
  val_bd4 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 4)
  val_bd5 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 5)
  val_bd6 = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                        basefunc = "BD", gamma = 6)
  val_bdmax = find.max.rotation.from.list(X, Y, rotation_list = rotation_list,
                                          basefunc = "BD", gamma = "max")
  
  pval_ed = (sum(pnull$ED >= val_ed) + 1) / (nperm + 1)
  pval_ed1 = (sum(pnull$ED1 >= val_ed1) + 1) / (nperm + 1)
  pval_ed2 = (sum(pnull$ED2 >= val_ed2) + 1) / (nperm + 1)
  pval_ed3 = (sum(pnull$ED3 >= val_ed3) + 1) / (nperm + 1)
  pval_ed4 = (sum(pnull$ED4 >= val_ed4) + 1) / (nperm + 1)
  pval_ed5 = (sum(pnull$ED5 >= val_ed5) + 1) / (nperm + 1)
  pval_ed6 = (sum(pnull$ED6 >= val_ed6) + 1) / (nperm + 1)
  pval_edmax = (sum(pnull$EDMAX >= val_edmax) + 1) / (nperm + 1)
  minpval_ed = min(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax)
  minpval_ed_simp = min(pval_ed, pval_ed1, pval_ed2, pval_ed4)
  
  pval_bd = (sum(pnull$BD >= val_bd) + 1) / (nperm + 1)
  pval_bd1 = (sum(pnull$BD1 >= val_bd1) + 1) / (nperm + 1)
  pval_bd2 = (sum(pnull$BD2 >= val_bd2) + 1) / (nperm + 1)
  pval_bd3 = (sum(pnull$BD3 >= val_bd3) + 1) / (nperm + 1)
  pval_bd4 = (sum(pnull$BD4 >= val_bd4) + 1) / (nperm + 1)
  pval_bd5 = (sum(pnull$BD5 >= val_bd5) + 1) / (nperm + 1)
  pval_bd6 = (sum(pnull$BD6 >= val_bd6) + 1) / (nperm + 1)
  pval_bdmax = (sum(pnull$BDMAX >= val_bdmax) + 1) / (nperm + 1)
  minpval_bd = min(pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax)
  minpval_bd_simp = min(pval_bd, pval_bd1, pval_bd2, pval_bd4)
  minpval = min(minpval_ed, minpval_bd)
  minpval_simp = min(minpval_ed_simp, minpval_bd_simp)
  return(list(c(pval_ed, pval_ed1, pval_ed2, pval_ed3, pval_ed4, pval_ed5, pval_ed6, pval_edmax,
                pval_bd, pval_bd1, pval_bd2, pval_bd3, pval_bd4, pval_bd5, pval_bd6, pval_bdmax),
              c(minpval_ed, minpval_bd, minpval),
              c(minpval_ed_simp, minpval_bd_simp, minpval_simp)))
}



