# ########################################## #
# some additional stuff that might be useful #
# ########################################## #


## calculate the prior the effects graph as described in the paper (according to the log-odss ratio matrix); see the paper for a more detailled description

calculateHiddenPrior = function(ratioMat){

  prior.hidden = t(exp(ratioMat)/(1+exp(ratioMat)))
  prior.hidden = rbind(prior.hidden,colMeans(prior.hidden,2))
  prior.hidden = apply(prior.hidden,2,FUN=function(c){c/sum(c)})
  
  return(prior.hidden)

}


## calculate the preferred attachement of effects to signals, according to the final effects graph prior (if Empirical Bayes has been applied)
## attachs the effect to the signal with the highest probability (prints a warning if more than one signal has the highest probability)

final_attachement = function(hidden_list){

  f_a = hidden_list[[length(hidden_list)]]
  f_a = apply(f_a,2,FUN=function(col){res = rep(0,length(col)); res[which(col==max(col))[1]] = 1;if(sum(res)>1){print("warning: more than one attachement");print(col)};return(res)})
  return(f_a)

}

## reduce the list of local maxima to a final one: draws an edge if it is present in at least 50% of the MCMC steps in the burnin phase (a burnin value has to be specified)
## Attention: This does NOT take into account dependencies between edges; for a more detailed discussion see the paper

getFinalTheta = function(theta_list, burnin){

  nrRuns = length(theta_list)
  signals = nrow(theta_list[[1]])

  start = burnin+1; end = nrRuns
  
  sEF_list = array(unlist(theta_list[start:end]),dim=c(signals,signals,(end-start+1)))
  sEF = rowSums(sEF_list,dims=2)
  theta_res = matrix(as.numeric(sEF>=(0.5*(end-start+1))),ncol=signals)

  return(theta_res)
  
}

prior.EgeneAttach.EB = function(ratioMat){
	Pe = t(calculateHiddenPrior(ratioMat))
	colnames(Pe) = c(colnames(ratioMat), "null")
	Pe
}
