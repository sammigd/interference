
make_hajek <- function(yhat_grp, df){
  #yhat_group is output of potential outcome ipw functon
  #sdf lets you get number per cluster
  
  z <- yhat_grp$yhat_group 
  po = apply(z, c(2,3), mean) #unadj pot outcome
  
  oe_z = yhat_grp$oe_yhat_group
  
  #hajek adj
  #sum the weights for each cluster
  wts = yhat_grp$wt_list #the sum of this over cluster is n_hz
  wts = apply(wts,c(2,3),sum) #this is n_hz
  
  oe_wts = apply(yhat_grp$oe_wt_list, 2, sum)
  
  #multiply each cluster yhat by the number of units in that cluster, then sum over all clusters
  size = df %>% group_by(neigh) %>% summarise(cluster_pop = n())

  #get the sum of the potential outcomes in each cluster (rather than the mean over individual units)
  sum_po = apply(z, c(2,3), "*", size$cluster_pop) #this is yhat * num units
  sum_po = apply(sum_po, c(2,3), sum) #this is summed over all clusters
  
  oe_sum_po = apply(oe_z, 2, '*', size$cluster_pop)
  oe_sum_po = apply(oe_sum_po, 2, sum)

  #print out the test / wts (shoudl be the estiamdn with hajek weight)
  haj = sum_po / wts
  
  oe_haj = oe_sum_po / oe_wts
  
  return(list(haj = haj, wts = wts, oe_haj = oe_haj, oe_wts = oe_wts))
} 

