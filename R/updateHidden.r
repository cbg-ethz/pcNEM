# #################################### #
# MCMC Sampling:  Empirical Bayes step #
# #################################### #

# Method description: 	this method updates the effects graph prior, depending on the signals graphs of the previous <changeHfreq> MCMC steps
#			for a detailed description of the algorithm, see the paper
# parameters: see runMCMC.m

updateHidden = function(ratio,prior.hidden,theta_list,signals,effects,i,changeHfreq){

  hidden_list = lapply(theta_list[(i-changeHfreq+1):i],FUN=function(theta){		# for each of the previous <changeHfreq> signals graphs...
		  f_m_k = log(prior.hidden) + rbind(t(ratio%*%theta),rep(0,effects))	# calculate f
		  f_m_k = exp(f_m_k)
		  S = colSums(f_m_k)							# for all effects: calculate sum over all signals+notAttached
		  res = t(t(f_m_k)/S)							# calculate the new prior, based on one signals graph
		  return(res)
		})
		
  hidden_mat = array(unlist(hidden_list),dim=c(signals+1,effects,changeHfreq))
  
  # calculate the mean over all <changeHfreq> signals graphs
  prior.hidden_new = apply(hidden_mat,c(1,2),sum)/changeHfreq
  
  return(prior.hidden_new)

}  

