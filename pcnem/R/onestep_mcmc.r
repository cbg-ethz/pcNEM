# ######################## #
# MCMC Sampling: Core part #
# ######################## #

# Method description: 	This method provides the core part of the MCMC sampling: suggestion of a new set of parameters and acceptance / rejection
#			For a more detailed description of the underlying formulas see the paper
# parameters: see runMCMC.m

onestep_mcmc = function(ratio,theta_i,graph2_i,ll_i,pP_i,nrSteps_i,changeList,logitprior.theta,prior.hidden,maxsteps_eminem,sd_val,probVal,ep){

    ## initialization

    # initially set acceptance variable to "rejected"
    acc_i = -1
    
    # initialize the new sampled graph with the old one -> will be changed in the following
    graph2_new = graph2_i
    
    # initially set transition probabilites (new to old / old to new) and sparseness priors (old and new) to zero (log-scale!)
    prop_n2o = 0
    prop_o2n = 0
    ep_o = 0
    ep_n = 0
    
    # simplify further calculation
    probVec = c(1-probVal,probVal)
    
    ## sampling of new signals graph + calculation of probabilites that are calculated for the sampled signals graph

    # randomly select the <sd_val> edges that could be changed (from the ones that may be changed)
    edgesToChange = sample(changeList,sd_val,replace=FALSE,prob=NULL)
    # randomly assing new values to them, according to <probVal> (i.e. not all of them will actually be changed)
    newValues = sample(c(0,1),length(edgesToChange),replace=TRUE,prob=probVec)
    graph2_new[edgesToChange] = newValues
    # for the calculation of the transition probabilites / sparseness prior, only actually changed edges are taken into account, since all other values will be the same and cancel each other out
    changes2one = length(which(graph2_new-graph2_i==1))
    changes2zero = length(which(graph2_new-graph2_i==(-1)))
    # calculate transition propabilities (when a new signals graph is sampled, the newly assigned values are independent of the old ones; number of changed edges is the same in the old as in the new signals graph)
    prop_o2n = changes2zero*log(probVec[1]) + changes2one*log(probVec[2]) # log( (probVec[1]^changes2zero) * (probVec[2]^changes2one) )
    prop_n2o = changes2one*log(probVec[1]) + changes2zero*log(probVec[2]) # log( (probVec[1]^changes2one) * (probVec[2]^changes2zero) )
    # calculate sparseness priors for the old and the new sampled signals graph; only take into account edges that actually may be changed
    ep_n = length(which(graph2_new[changeList]==0))*log(probVec[1]) + length(which(graph2_new[changeList]==1))*log(probVec[2])
    ep_o = length(which(graph2_i[changeList]==0))*log(probVec[1]) + length(which(graph2_i[changeList]==1))*log(probVec[2])

    ## calculation of the corresponding local maximum + calculation of probabilites that are calculated for the local maximum

    # apply EMiNEM
    res_EMiNEM = EMforMCMC(ratio,graph2_new,logitprior.theta,prior.hidden,maxsteps_eminem)
    theta_new = res_EMiNEM[["res"]]
    nrSteps_new = res_EMiNEM[["nrSteps"]]
    
    # calculate the log-likelihood for the local maximum
    ll_new = calcLikelihood(ratio,theta_new,prior.hidden)
    
   # calculate the prior for the local maximum -> see the paper for a derivation of this formula
    pP_new = sum(theta_new[changeList]*logitprior.theta[changeList])
    
    ## acceptance / rejection step
    
    # (adapted) Metropolis-Hastings - see paper for further explanation
    u = (ll_new+pP_new+prop_n2o)-(ll_i+pP_i+prop_o2n)+ep*(ep_n-ep_o)
    rand_u = log(runif(1,0,1))
    
    # if ... -> acception, otherwisse rejection (i.e. the old values will be returned, without beeing updated)
    if(rand_u < u){
    
      # distinguish between the acception of a new local maximum (1) and the acception of a new sampled signals graph in the range of the old local maximum (0)
      if(all(theta_i==theta_new)){
	acc_i = 0
      }else{
	acc_i = 1    
      }
    
      # if the new sampled signals graph has been accepted -> update the variables
      theta_i = theta_new
      graph2_i = graph2_new
      ll_i = ll_new
      pP_i = pP_new
      nrSteps_i = nrSteps_new

    }
    
    # return the variables - whether updated or not
    return(list(theta=theta_i,graph2=graph2_i,ll=ll_i,pP=pP_i,nrSteps=nrSteps_i,accepted=acc_i))

}

