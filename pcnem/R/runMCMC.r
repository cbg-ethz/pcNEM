# ############################# #
# MCMC Sampling: Initialization #
# ############################# #

# Method description: this method initializes the parameters for the MCMC sampling, reads in the data, organizes the results, etc

# Important: Here, the arrangement of the ratio matrix is not the same as in the paper, it is transposed! (i.e. here: R_kj, in the paper: R_jk, for any signal j and effect k)

# Input parameters:

# ratio			numeric matrix: the data matrix (log-odds ratios), must include row- and column-names - rows = effects, columns = signals
# nrRuns		positive number: the number of MCMC steps
# theta_init		binary matrix 0/1, ncol(ratio) x ncol(ratio): the initial signals graph - if not provided by the user, it is randomly sampled based on the signals graph prior (includes all edges with a probability > 0.5)
# prior.theta		numeric matrix, values between 0 and 1, ncol(ratio) x ncol(ratio): the prior for the signals graph - if not provided by the user, it is set to 0.5 for all edges
# prior.hidden		numeric matrix, values between 0 and 1, ncol(ratio) x nrow(ratio)+1, columns have to sum up to one:the prior for the effects graph - if not provided by the user, attachement probabilities are assigned uniformily (1/(|signals|+1))
# maxsteps_eminem	positive number: the maximal number of iterations for the EM algorithm
# sd_val	  	positive number, between 1 and ncol(ratio)*(ncol(ratio)-1): number of edges to change in one MCMC step; see paper for the author's choice
# changeHfreq		positive number, nrRuns must be a multiple: the Empirical Bayes step is performed every <changeHfreq> steps (see paper for the author's choice); set >= nrRuns (or leave to default) to exclude this step
# probVal		positive number, between 0 and 1: expected fraction of edges present in the signals graph; see paper for the author's choice
# ep			positive number: weight of sparsitiy prior for calculation of acceptance rate; see paper for the author's choice

# Return values

# theta_list		the list of local maxima accepted during the sampling - length(theta_list)=nrRuns
# graph2_list		the list of sampled signals graphs, accepted during the sampling;graph2_list[[1]] is the initial signals graph - length(graph2_list)=nrRuns
# ll_list		the list of log-likelihoods, corresponding to the local maxima - length(ll_list)=nrRuns
# pP_list		the list of log-posteriors, corresponding to the local maxima - length(pP_list)=nrRuns
# nrSteps_list		the list of number of steps needed by the EM-algorithm to find the local maximum - length(nrSteps_list)=nrRuns
# acc_list		the list that indicates whether the corresponding sampled signals graph has been accepted (new local maximum (1), same local maximum (0)) or rejected(-1) in the MCMC sampling process - length(acc_list)=nrRuns
# hidden_list 		the list of effects graphs priors, adapted during the sampling process;hidden_list[[1]] is the initial prior  - length(hidden_list)=nrRuns/changeHfreq or 1 (if Empirical Bayes is not applied)
# prior.theta		the prior for the signals graph (for the sake of completeness)
# changeList		the list of edges that can be changed during the sampling process (for the sake of completeness)

runMCMC = function(ratio=NULL,nrRuns=100,theta_init=NULL,prior.theta=NULL,prior.hidden=NULL,maxsteps_eminem=1000,sd_val=1,changeHfreq=NULL,probVal=0.2,ep=0){

  # check if ratio matrix is valid
  if(is.null(ratio)){print("Error: no valid ratio matrix provided");return(list())} ## check on non-numeric values (Inf,...), too?
  
  # extract number of signal- and effect-node
  signals = ncol(ratio)
  effects = nrow(ratio)

  # set changeHfreq to a valid value for further calculation
  if(is.null(changeHfreq)){changeHfreq=nrRuns+1}

  # check if ... are valid
  if(nrRuns<=0){print("Error: no valid value for nrRuns provided");return(list())}
  if(maxsteps_eminem<=0){print("Error: no valid value for maxsteps_eminem provided");return(list())}
  if(changeHfreq<=0 || (changeHfreq<nrRuns && nrRuns%%changeHfreq!=0)){print("Error: no valid value for changeHfreq provided");return(list())}
  if(probVal<0 || probVal>1){print("Error: no valid value for probVal provided");return(list())}
  if(sd_val<=0 || sd_val > signals*(signals-1)){print("Error: no valid value for sd_val provided");return(list())}
  if(ep<0){print("Error: no valid value for ep provided");return(list())}

  # set initial signals graph and priors, if they are not provided - otherwise check for consistence
  if (is.null(prior.theta)){
	  prior.theta = matrix(1/2,nrow=signals,ncol=signals)
	  diag(prior.theta)=1
	  rownames(prior.theta) = colnames(prior.theta) = colnames(ratio)
  }else{
	  if(nrow(prior.theta)!=signals || ncol(prior.theta)!=signals){print("Error: the dimensions of prior.theta are not consistent with the ratio matrix");return(list())}
	  if(min(prior.theta)<0 || max(prior.theta)>1){print("Error: prior.theta contains invalid values");return(list())}
  }
  if (is.null(theta_init)){
	  theta_init = matrix(as.numeric(prior.theta>=0.5),nrow=signals,ncol=signals)
	  diag(theta_init) = 1
	  rownames(theta_init) = colnames(theta_init) = colnames(ratio)
  }else{
	  if(nrow(theta_init)!=signals || ncol(theta_init)!=signals){print("Error: the dimensions of theta_init are not consistent with the ratio matrix");return(list())}
	  if(length(which(theta_init!=0))!=length(which(theta_init==1))){print("Error: theta_init contains invalid values");return(list())}
  }
  if (is.null(prior.hidden)){
	  prior.hidden = matrix(1/(signals+1),nrow=signals+1,ncol=effects)
	  rownames(prior.hidden) = c(colnames(ratio),"notAttached")
  }else{
	  if(nrow(prior.hidden)!=(signals+1) || ncol(prior.hidden)!=effects){print("Error: the dimensions of prior.hidden are not consistent with the ratio matrix");return(list())}
	  if(min(prior.hidden)<0 || max(prior.hidden)>1){print("Error: prior.hidden contains invalid values");return(list())}
	  if(any(round(colSums(prior.hidden),digits=10)!=1)){print("Error: columns of prior.hidden do not sum up to one");return(list())}
  }
  
  # the EM-algorithm needs the logit prior
  logitprior.theta = log(prior.theta/(1-prior.theta))
  
  # define which edges in the signals graph are not fixed to be absent / present by the prior - and check for a reasonable value
  changeList = which(prior.theta!=0 & prior.theta!=1)
  if(length(changeList)==0){print("Error: all edges are fixed by the prior choice");return(list())}
  
  # start sampling process
  res = mcmc(ratio=ratio,nrRuns=nrRuns,theta_init=theta_init,changeList=changeList,logitprior.theta=logitprior.theta,prior.hidden=prior.hidden,maxsteps_eminem=maxsteps_eminem,sd_val=sd_val,signals=signals,effects=effects,changeHfreq=changeHfreq,probVal=probVal,ep=ep)
  
  # add additional information to the output
  res[["prior.theta"]] = prior.theta
  res[["changeList"]] = changeList
  
  return(res)

}





