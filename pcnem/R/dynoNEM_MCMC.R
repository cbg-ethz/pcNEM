# ###########################################################################
# Author:
# Date:
# Purpose: To run the MCMC based dynoNEM.
# 
# 
# Related files: wrapper.c, simulation.R, MCMCdocu.txt, functionsMCMC.h
#		 netlearn.h, LogfilesMCMCrun.txt, and ouputMCMCdynoNEM
# 		 directory
# ###########################################################################

#dyn.load("~/workingAt/trunk/dynoNEM/dynoNEMMCMC/wrapper.so")

dynoNEM_MCMC = function(D, SAMPLE=500000, BURNIN=1000000, initial=NULL, priorNet=NULL, priorE=NULL, delta=1, inv.nu=0.01, theta=1, type=c("CONTmLLBayes", "CONTmLL", "mLL"), nrep=4, alpha=0.1, beta=0.2, seed = 1234){	

	if(any(is.na(D) | is.infinite(D)))
		stop("data contains NA or Inf")
	T = dim(D)[1] # get the value of T
	nsgenes = dim(D)[3] # get the value of nsgenes
	negenes = dim(D)[2] # get the value of negenes
	type = match.arg(type, several.ok=FALSE) # accepting the agrguement for type
	if(type == "CONTmLLBayes")
		int.type = 0
	else if(type == "CONTmLL")
		int.type = 1
	else if(type == "mLL")
		int.type = 2
	
	#if(is.null(priorNet)) # The prior net (networkPrior)
	#	priorNet = matrix(0, ncol=nsgenes, nrow=nsgenes) # as an empty matrix
	#else{ # conditions for un-acceptability
	if(!is.null(priorNet)){
		if(NCOL(priorNet) != nsgenes | NROW(priorNet) != nsgenes)
			stop("priorNet has to be a # S-genes x # S-genes matrix")
		if(any(priorNet > T) | any(priorNet < 0))
			stop("All entries in priorNet have to be between 0 and T")
		if(inv.nu == 0)
			inv.nu = 0.01
	}
	if(is.null(priorE)){ # The E-gene prior (Egene_prior) uniform prior
		priorE = matrix(1/nsgenes, ncol=nsgenes, nrow=negenes)		
	}		
	else{ # conditions for un-acceptability
		if(NCOL(priorE) != nsgenes | NROW(priorE) != negenes)
			stop("priorE has to be a # E-genes x # S-genes matrix")
		if(any(priorE > 1) | any(priorE < 0))
			stop("All entries in priorE have to be probabilities")
	}
	if(ncol(priorE) == nsgenes){
		priorE = cbind(priorE, as.matrix(rep(delta/nsgenes, negenes))) # add a "null" S-gene
		priorE = t(apply(priorE, 1, function(x) x/sum(x))) # re-normalize E-gene prior
	}
	if (!is.null(initial)){ # initial nework matrix
		if(NCOL(initial) != nsgenes | NROW(initial) != nsgenes)
			stop("initial has to be a # S-genes x # S-genes matrix")
		if(any(initial < 0))
			stop("All entries in initial have to be positive integers")
	}
	if(any(inv.nu < 0)) # nu (priorScale) parameter
		stop("Parameter inv.nu (inverse prior scale) has to be non-negative")		
	net = matrix(0, ncol=nsgenes, nrow=nsgenes)
	dimnames(net) = list(dimnames(D)[[3]], dimnames(D)[[3]])
	if(is.null(initial)){
		if(type == "CONTmLL")
			logD = log(D/(1-D))
		else
			logD = D
		myD = matrix(0, nrow=dim(D)[2], ncol=dim(D)[3])					
		for(t in 1:T){
			myD = myD + logD[t,,]
		}									
		colnames(myD) = colnames(D[1,,])
		rownames(myD) = as.character(1:nrow(myD))					
		if(type != "mLL")
			control = set.default.parameters(unique(colnames(myD)), type="CONTmLLMAP", Pm=priorNet, Pe=priorE[,1:nsgenes], trans.close=FALSE, backward.elimination=TRUE, lambda=inv.nu, delta=delta)					
		else if(type == "mLL")
			control = set.default.parameters(unique(colnames(myD)), Pm=priorNet, Pe=priorE[,1:nsgenes], lambda=inv.nu, para=c(alpha,beta), delta=delta)
		initial = nem(myD, control=control, verbose=FALSE)		
		#sel = as.numeric(initial$selected) # ACHTUNG: automatische Feature Selection
		initial = as(initial$graph, "matrix")
		cat("initial network:\n")
		print(initial)		
	}
	sel = 1:negenes
# Call the MCMC function
	if(type == "CONTmLLBayes")
		DD = exp(D[,sel,]) # IMPORTANT: code works with original p-value densities
	else
		DD = D[,sel,]
	negenes = dim(DD)[2]		
	RetVec = .Call("MCMCrunWrapper",
		as.integer(SAMPLE), #1
		as.integer(BURNIN), #2	
		as.numeric(initial), # 5
		as.integer(nsgenes), #6
		as.integer(negenes), #7
		as.integer(T), #8
		as.numeric(DD), #9
		as.numeric(priorNet), #10
		as.numeric(priorE), #11
		as.numeric(inv.nu), # 12
		as.numeric(theta),
		as.integer(int.type), #13
		as.integer(nrep), #14
		as.double(alpha), #15
		as.double(beta), #16
		as.integer(seed), #17	
		PACKAGE="nem"
	)	
	all.likelihoods = RetVec$allLikelihoods
	network = matrix(RetVec$net.res, ncol=nsgenes)
	SDconf= matrix(RetVec$sdmat, ncol=nsgenes)
	dimnames(network) = dimnames(net)
	dimnames(SDconf) = dimnames(net)
	cat("Finished MCMC run....................\n")
	print(network)
	print(SDconf)
				
	return(list(network=network, all.likelihoods=all.likelihoods, SDconf=SDconf))	 
}

dynoNEM.perturb = function(Psi, T, k){
	if(NROW(Psi) != NCOL(Psi))
		stop("net has to be a quadratic, weighted adjacency matrix")
	if(k > NROW(Psi) | k < 1)
		stop("k has to be between 1 and #S-genes")
	if(T < 0)
		stop("T has to be > 0")
	nsgenes = nrow(Psi)
	perturb_prob = matrix(0, nrow=nsgenes, ncol=T)
	rownames(perturb_prob) = rownames(Psi)
	colnames(perturb_prob) = as.character(1:T)
	nsgenes = NCOL(Psi)
	for(t in 1:T){
		for(s in 1:nsgenes){
			perturb_prob[s, t] = 0; # default: s is unperturbed		
			perturb_prob[k, t] = 1; # the perturbed gene is always inactive		
			if(s != k){
				for(p in 1:nsgenes){
					if(Psi[p, s] != 0 && abs(Psi[p, s]) <= t){ # p is a parent
						if(t > 1){
							parent_perturb_prob = perturb_prob[p, t-1];						
						}
						else
							parent_perturb_prob = (p == k)*1;
						if(parent_perturb_prob == 1){
							perturb_prob[s, t] = 1;
							break;
						}
					}
				}
			}
		}
	}   	
	perturb_prob
}

dynoNEM.posteriorEGenePos = function(net, D, priorE=NULL, delta=1, type=c("CONTmLLBayes", "CONTmLL", "mLL"), nrep=4, alpha=0.1, beta=0.2){
	if(is.null(dimnames(D)))
		stop("Each dimension of your data array should have names (first dimension: time points, e.g. as.character(1:T); seond dimension: E-genes; third dimension: S-genes)")
	T = dim(D)[1]
	nsgenes = dim(D)[3]
	negenes = dim(D)[2]
	type = match.arg(type, several.ok=FALSE)	
	if(is.null(priorE))
		priorE = matrix(1/nsgenes, ncol=nsgenes, nrow=negenes)
	else{
		if(NCOL(priorE) != nsgenes | NROW(priorE) != negenes)
			stop("priorE has to be a # E-genes x # S-genes matrix")
		if(any(priorE > 1) | any(priorE < 0))
			stop("All entries in priorE have to be probabilities")
	}		
	priorE = cbind(priorE, as.matrix(rep(delta/nsgenes, negenes))) # add "null" S-gene
	priorE = t(apply(priorE, 1, function(x) x/sum(x))) # re-normalize S-gene prior
	if(type == "CONTmLLBayes")
		D = exp(D) # IMPORTANT: code below works with original p-value densities
	perturb_prob = array(0, dim=c(nsgenes,nsgenes, T))
	for(k in 1:nsgenes){
		perturb_prob[k,,] = dynoNEM.perturb(net, T, k)
	}
	epos = matrix(0, ncol=(nsgenes+1), nrow=negenes)
	loglik = 0	
	loglik0 = double(nsgenes + 1)
	for(i in 1:negenes){
		loglik_tmp = 0
		for(s in 1:(nsgenes + 1)){
			tmp = 0
			for(k in 1:nsgenes){
				for(t in 1:T){
					if(priorE[i, s] > 0 & s <= nsgenes){
						if(type == "CONTmLLBayes")
							tmp = tmp + log((D[t, i, k]*perturb_prob[k, s, t] + (1-perturb_prob[k, s, t])*1) * priorE[i, s] + 1e-100)
						else if(type == "CONTmLL")
							tmp = tmp + log((D[t, i, k]*perturb_prob[k, s, t] + (1-perturb_prob[k, s, t])*(1 - D[t, i, k])) * priorE[i, s] + 1e-100)
						else if(type == "mLL")
							tmp = tmp + log( ((1-beta)^(D[t, i, k]*perturb_prob[k, s, t]) * beta^((nrep-D[t, i, k])*perturb_prob[k, s, t]) +
												alpha^(D[t, i, k]*(1-perturb_prob[k, s, t])) * (1-alpha)^((nrep-D[t, i, k])*(1-perturb_prob[k, s, t]))) * priorE[i, s])
					}
					if(priorE[i, s] > 0 & s > nsgenes){ # log-likelihoods for "null" S-gene (unconnected to the rest of S-genes)
						if(type == "CONTmLLBayes")
							tmp = tmp + log(priorE[i, s] + 1e-100)
						else if(type == "CONTmLL")
							tmp = tmp + log((1 - D[t, i, k]) * priorE[i, s] + 1e-100)
						else if(type == "mLL")
							tmp = tmp + log(alpha^D[t, i, k] * (1-alpha)^(nrep-D[t, i, k]) * priorE[i, s])
					}
				}
			}		
			loglik0[s] = tmp					
			epos[i, s] = loglik0[s]
		}
		maxidx = which.max(loglik0)
		loglik = loglik + loglik0[maxidx] + log(sum(exp(loglik0 - loglik0[maxidx])))
	}
	LLperGene = rowSums(epos)
	epos = exp(epos)
	epos = epos / rowSums(epos)	
	dimnames(epos) = list(dimnames(D)[[2]], c(dimnames(D)[[3]], "null"))	
	map.pos = apply(epos, 1, function(x){ s=length(colnames(epos)[which(x == max(x))]); data.frame(attachment=colnames(epos)[which(x == max(x))], probability=rep(max(x), s))})		
	map.pos = as.data.frame(do.call("rbind", map.pos))
	map.pos = map.pos[map.pos$attachment != "null",]
	list(pos.probability=epos, map.pos=map.pos, LLperGene=LLperGene, mLL=loglik)
}

