cal.RATS.null = function(pnull){
  nperm = dim(pnull)[1]
  ## permute the sample to get null distribution of min p value
  minp_null = numeric(nperm)
  minped_null = numeric(nperm)
  minpbd_null = numeric(nperm)
  minp_simp_null = numeric(nperm)
  minped_simp_null = numeric(nperm)
  minpbd_simp_null = numeric(nperm)
  
  for (t in 1:nperm){
    #if (t %% 100 == 0) {print(paste0("rep: ", t)); gc(verbose = TRUE)}
    
    pval_rep_ed = mean(pnull$ED >= pnull$ED[t])
    pval_rep_ed1 = mean(pnull$ED1 >= pnull$ED1[t])
    pval_rep_ed2 = mean(pnull$ED2 >= pnull$ED2[t])
    pval_rep_ed3 = mean(pnull$ED3 >= pnull$ED3[t])
    pval_rep_ed4 = mean(pnull$ED4 >= pnull$ED4[t])
    pval_rep_ed5 = mean(pnull$ED5 >= pnull$ED5[t])
    pval_rep_ed6 = mean(pnull$ED6 >= pnull$ED6[t])
    pval_rep_edmax = mean(pnull$EDMAX >= pnull$EDMAX[t])
    
    pval_rep_bd = mean(pnull$BD >= pnull$BD[t])
    pval_rep_bd1 = mean(pnull$BD1 >= pnull$BD1[t])
    pval_rep_bd2 = mean(pnull$BD2 >= pnull$BD2[t])
    pval_rep_bd3 = mean(pnull$BD3 >= pnull$BD3[t])
    pval_rep_bd4 = mean(pnull$BD4 >= pnull$BD4[t])
    pval_rep_bd5 = mean(pnull$BD5 >= pnull$BD5[t])
    pval_rep_bd6 = mean(pnull$BD6 >= pnull$BD6[t])
    pval_rep_bdmax = mean(pnull$BDMAX >= pnull$BDMAX[t])
    
    minped_null[t] = min(pval_rep_ed, pval_rep_ed1, pval_rep_ed2, pval_rep_ed3, pval_rep_ed4, 
                         pval_rep_ed5, pval_rep_ed6, pval_rep_edmax)
    minpbd_null[t] = min(pval_rep_bd, pval_rep_bd1, pval_rep_bd2, pval_rep_bd3, pval_rep_bd4,
                         pval_rep_bd5, pval_rep_bd6, pval_rep_bdmax)
    minp_null[t] = min(minped_null[t], minpbd_null[t])
    
    minped_simp_null[t] = min(pval_rep_ed, pval_rep_ed1, pval_rep_ed2, pval_rep_ed4)
    minpbd_simp_null[t] = min(pval_rep_bd, pval_rep_bd1, pval_rep_bd2, pval_rep_bd4)
    minp_simp_null[t] = min(minped_simp_null[t], minpbd_simp_null[t])
  }
  
  minp_rep = tibble(minped_null, minpbd_null, minp_null, minped_simp_null, minpbd_simp_null, minp_simp_null)
  colnames(minp_rep) = c("minP_ED", "minP_BD", "minP", "minP_ED_simp", "minP_BD_simp", "minP_simp")
  return(minp_rep)
}


