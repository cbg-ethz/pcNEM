# ############################ #
# MCMC Sampling:  EM algorithm #
# ############################ #

# Method description: 	this method conducts the EM algorithm and offers a function to calculate the log-likelihood
#			for a detailed description of the algorithms, see the paper
# parameters: see runMCMC.m


# Dimensions of the matrices: theta: 	signals x signals (the same as in the paper)
#					ratio: effects x signals (different from the paper, transposed!)
#					logitprior.theta: signals x signals
#					prior.hidden: (signals+1) x effects



# apply EMiNEM to calculate the local maximum of theta - one iteration
onestep_EMiNEM = function(theta,ratio,logitprior.theta,prior.hidden){

    signals = nrow(theta)

    # calculate the closed update formula, as derived in the paper (for more informations on these formulas, see the paper)
    fkdn = exp(ratio %*% theta)
    Aj = apply(t(prior.hidden)*cbind(fkdn,rep(1,nrow(fkdn))),1,sum)
    thetaprob = t(ratio/Aj) %*% (fkdn*t(prior.hidden[1:signals,])) + logitprior.theta

    # the new, improved signals graphs has edges where thetaprob>0 and no edges, where thetaprob<0
    thetanew = as.numeric(thetaprob>0)
    dim(thetanew) = dim(theta)
    
    # make sure that edges with prior probabilities of 0 or 1 are kept at 0 or 1
    thetanew[which(logitprior.theta==(-Inf))]=0
    thetanew[which(logitprior.theta==(Inf))]=1

    return(thetanew)
    
}



# apply EMiNEM to calculate the local maximum of theta - framework
EMforMCMC = function(ratio,theta,logitprior.theta,prior.hidden,maxsteps){

    # initialization of some variables
    signals = ncol(ratio)
    effects = nrow(ratio)
    steps = 1
    protocol = list()
    
    # stepwise improvement of the signals graph - until maximum number of steps is reached or no improvement is made any more
    repeat{
    
	# add the current maximum to the results-list
	protocol[[steps]] = theta
	# stop, if maximum number of steps is reached
	if (steps==(maxsteps+1)) break()
	# do another improvement-step
	thetanew =  onestep_EMiNEM(theta,ratio,logitprior.theta,prior.hidden)
	# stop, if no improvement is made any more
	if (all(thetanew==theta)) break()
	# otherwise, update variables and continue
	rownames(thetanew) = rownames(theta)
	colnames(thetanew) = colnames(theta)		
	theta = thetanew
	steps = steps + 1
	
    }
    
    # return the final signals graphs and the number of steps needed
    return(list(res=protocol[[steps]],nrSteps=steps))
}



# calculate the log-likelihood for <theta>, based on the current prior for the effects graph - see the paper for a derivation of this formula
calcLikelihood = function(ratio,theta,prior.hidden){

    p1 = t(prior.hidden[1:nrow(prior.hidden)-1,])
    p2 = exp(ratio%*%theta)

    c = rowSums(p1*p2)

    return(sum(log(c)))
    
}

