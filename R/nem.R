set.default.parameters = function(Sgenes, ...){         
    map = as.list(Sgenes)
    names(map) = Sgenes 
    control = list(Sgenes=Sgenes,type="mLL",para=c(0.13,0.05),hyperpara=c(1,9,9,1),Pe=NULL,Pm=NULL,Pmlocal=NULL, Pm.frac_edges=0.2,
			local.prior.size=length(Sgenes),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE, selEGenes.method="regularization",
			trans.close=TRUE, backward.elimination=FALSE, mode="continous_Bayesian", lambda.intervention=4, lambda.no_intervention=4,  df.intervention=4.4, 
			df.no_intervention=4.4, nu.intervention=0.6, nu.no_intervention=0.95, scale.intervention=0.023, scale.no_intervention=0.023, map=map, outputdir="QualityControl", 
			debug=FALSE, mc.cores=8,
			mcmc.nsamples=1e6, mcmc.nburnin=1e6, mcmc.seed=1234, mcmc.hyperprior=1,
			eminem.maxsteps=1000, eminem.sdVal=1, eminem.changeHfreq=NULL,
			prob.cutoff=0.5,
            pcombi=FALSE, iterations = 2e4, stepsave = 1e2, revallowed = 1,
            AcceptRate = NULL, Temp = 50, AdaptRate = 0.3, noiseEst = TRUE,
            moveprobs = c(0.6,0.4), moveprobsNoise = c(0.5,0.5),
            sigma = diag(x=1,2)/10)
    args = list(...)
    na = names(args)
    for(a in na){
        control[[a]] = args[[match(a, na)]]
    }
    control
}

nem <- function(D, inference="nem.greedy",models=NULL,control=set.default.parameters(setdiff(unique(colnames(D)),"time")), verbose=FALSE){
#------------------------------
# Sanity checks              
if(inference == "nem.greedyMAP")
    control$type = "CONTmLLMAP"
  
if (!(inference %in% c("pairwise", "triples", "search","ModuleNetwork","ModuleNetwork.orig", "nem.greedy","nem.greedyMAP","BN.greedy","BN.exhaustive","mc.eminem","dynoNEM","AdaSimAnneal"))) 
    stop("\nnem> argument 'inference' is not valid\n")
if (!(control$type %in% c("mLL", "FULLmLL", "CONTmLL","CONTmLLBayes","CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio","depn"))) 
    stop("\nnem> argument 'type' is not valid")
if (is.null(control$para) & is.null(control$hyperpara) & !(control$type %in% c("CONTmLL", "CONTmLLBayes", "CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio","depn")))
    stop("\nnem> provide either 'para' or 'hyperpara'\n")
if (control$type == "mLL" & is.null(control$para)) 
    stop("\nnem> provide argument 'para'\n")
if (control$type == "FULLmLL" & is.null(control$hyperpara)) 
    stop("\nnem> provide argument 'hyperpara'\n")
if (!is.null(control$hyperpara)) {
    if (length(control$hyperpara) != 4) 
        stop("\nnem> 'hyperpara' is not a vector of length 4")
    if (!all(control$hyperpara > 0)) 
        stop("\nnem> 'hyperpara' must be >0")
}
if (!is.null(control$para)) {
    if (length(control$para) != 2) 
        stop("\nnem> 'para' is not a vector of length 2")
    if (any(control$para < 0) | any(control$para > 1)) 
        stop("\nnem> 'para' must be in [0,1]")
}
if(control$lambda < 0) control$lambda <- abs(control$lambda)

if(control$pcombi==TRUE){
    map<-control$map
    Sgenes<-colnames(map)
}else{
  if(is(D, "list"))
    Sgenes = setdiff(unique(colnames(D[[1]])), "time")
  else
    Sgenes = setdiff(unique(colnames(D)), "time")
}

if(control$selEGenes){  
    if(inference %in% c("BN.greedy","BN.exhaustive","depn", "search"))
        stop("No automatic tuning of reporter filtering implemented for this inference method so far!\n")
    return(nem.featureselection(D, inference, models, control, verbose))
}

#------------------------------
# PAIRWISE                     

if (inference == "pairwise"){   
    if (is.null(control$local.prior.size))
        Pmlocal <- NULL 
    else{ 
        if (control$local.prior.size <= 0 | control$local.prior.bias <= 0) 
            stop("\nnem> local prior parameters invalid")
            Pmlocal <- local.model.prior(control$local.prior.size,length(Sgenes),control$local.prior.bias)      
        }
        result <- pairwise.posterior(D,control,verbose)    
    }

#------------------------------
# MODULE NETWORK                       
else if(inference == "ModuleNetwork"){
    result <- moduleNetwork(D,control,verbose=verbose)  
}
else if(inference == "ModuleNetwork.orig"){
	result <- moduleNetwork.orig(D,control,verbose=verbose)  
}

#------------------------------
# TRIPLES                     

else if (inference == "triples"){

        result <- triples.posterior(D,control,verbose)
        #A.t <- transitive.lp(A)

    }

#------------------------------
# GREEDY                     
else if(inference == "nem.greedy"){
    result <- nem.greedy(D,initial=models[[1]],control, verbose=verbose)    
}

else if(inference == "nem.greedyMAP"){
    result <- nem.greedyMAP(D,control, verbose=verbose) 
}
# MC.EMINEM
else if(inference == "mc.eminem"){
	cat("mc.eminem algorithm\n\n")
	if(!control$type %in% c("CONTmLLBayes","CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio"))
		stop("Likelihood type has to be one of 'CONTmLLBayes', 'CONTmLLMAP', 'CONTmLLDens', 'CONTmLLRatio'")
	if(!is.null(control$Pe)){
		if(ncol(control$Pe) == length(Sgenes)){				
			control$Pe = cbind(control$Pe, double(nrow(D)))  			
			control$Pe[,ncol(control$Pe)] = control$delta/length(Sgenes)
			control$Pe = control$Pe/rowSums(control$Pe)		
		}
		Pe = t(control$Pe)
	}
	else
		Pe = NULL	
	result1 = runMCMC(D, nrRuns=control$mcmc.nsamples + control$mcmc.nburnin, theta_init=models[[1]], prior.theta=control$Pm, prior.hidden=Pe, maxsteps_eminem=control$eminem.maxsteps, sd_val=control$eminem.sdVal, probVal=control$Pm.frac_edges, ep=control$lambda, changeHfreq=control$eminem.changeHfreq)		
	result1$avg = getFinalTheta(result1$theta_list, burnin=control$mcmc.nburnin)
	dimnames(result1$avg) = list(control$Sgenes, control$Sgenes)
	Phi = (result1$avg > control$prob.cutoff)*1
	result = score(list(Phi), D, control=control, verbose=verbose)		
	result$mappos = result$mappos[[1]]
	result$local.maxima = result1$theta_list
	result$graphs.sampled = result1$graph2_list
	result$mLL = result1$pP_list
	result$EB = result1$hidden_list
	result$avg = result1$avg
	result$acc_list = result1$acc_list
	class(result) = "mc.eminem"
}
# dynoNEM with MCMC
else if(inference == "dynoNEM"){
	cat("DynoNEM algorithm for time series data\n\n")
	if(length(dim(D)) < 2)
		stop("Data has to be an array of dimension T x #E-genes x #S-genes for the continous and T x #E-genes x (#S-genes * #replicates) for the discrete case")
	if(control$type %in% c("mLL", "FULLmLL")){
		nreps = sapply(control$Sgenes, function(s) sum(dimnames(D)[[3]] == s))
		if(var(nreps) != 0)
			stop("dynoNEM currently only supports the same number of replicates per S-gene")
		nrep = nreps[1]
		if(nrep > 1){
			Dnew = array(dim=c(dim(D)[1], dim(D)[2], length(control$Sgenes)))
			dimnames(Dnew) = list(as.character(1:dim(D)[1]), dimnames(D)[[2]], control$Sgenes)
			for(s in control$Sgenes){
				for(t in 1:dim(D)[1]){
					Dnew[t,,s] = rowSums(D[t,, colnames(D[t,,]) == s])
				}			
			}
		}
		else
			Dnew = D
	}
	else{
		nrep=1
		Dnew = D
	}
	if(dim(Dnew)[3] != length(control$Sgenes))
		stop("Dimensions of 'D' disagree with provided list of S-genes. 'D' has to be an array of dimension T x #E-genes x #S-genes (after summarization of replicates)")
	result = dynoNEM_MCMC(Dnew, SAMPLE=control$mcmc.nsamples, BURNIN=control$mcmc.nburnin, initial=models[[1]], priorNet=control$Pm, priorE=control$Pe, delta=control$delta, inv.nu=control$lambda, theta=control$mcmc.hyperprior, type=control$type, nrep=nrep, alpha=control$para[1], beta=control$para[2], seed=control$mcmc.seed)
	names(result)[names(result) == "all.likelihoods"] = "mLL"
	take = which(result$mLL != Inf)
	if(length(take) > 0 & control$mcmc.nsamples*control$mcmc.nburnin > 0){
		try(plot(take, result$mLL[take], type="l", main=paste("posterior log-likelihood along MCMC sampling"), xlab="step", ylab="log likelihood"), silent=TRUE)
		try(abline(v=control$mcmc.nburnin, lty=3), silent=TRUE)		
	}	
	result$avg = result$network
	result$avg[result$avg - 2*result$SDconf < 0] = 0 # filter low confidence edges
	result$graph = as(result$avg, "graphNEL")
	result$ppost = exp(result$mLL)
	pos = dynoNEM.posteriorEGenePos(result$network, Dnew, priorE=control$Pe, delta=control$delta, type=control$type, nrep=nrep, alpha=control$para[1], beta=control$para[2])
	result$ep = pos$pos.probability
	Theta = apply(result$ep,1,function(e) e ==max(e))
	result$mappos = apply(Theta,1,which)
	result$control = control	
	result$selected = unique(unlist(result$mappos[control$Sgenes]))	
	result$para = control$para
	result$T = dim(Dnew)[1]
	class(result) = "dynoNEM"
}

#------------------------------
# ADAPTIVE SIMULATED ANNEALING

else if(inference == "AdaSimAnneal"){
    nSgenes <- length(Sgenes)
    if(is.null(models)){
        startDAG <- matrix(numeric(nSgenes*nSgenes),nrow=nSgenes) # starting DAG is empty say
        dimnames(startDAG) <- list(Sgenes,Sgenes)
    }else{
        startDAG <- models[[1]]
        diag(startDAG) <- 0
        dimnames(startDAG) <- list(Sgenes,Sgenes)
    }  
    tempRes <- nem.AdaSA(nSgenes,startDAG,D,control)
    control$para <- c(tempRes$typeI,tempRes$typeII)
    result <- c(score(list(tempRes$maxDAG),D,control=control, verbose),
    list(models = tempRes$graphs, LLscores = tempRes$LLscores,Temperatures = tempRes$Temp, AcceptRates = tempRes$AcceptRate,
    minIter = tempRes$minIter, transformAR=tempRes$transformAR, typeIEst = tempRes$typeI, typeIIEst = tempRes$typeII, typeIchain = tempRes$AlphaVals, typeIIchain = tempRes$BetaVals))
}

#------------------------------
# SEARCH                       

else if (inference == "search"){    				
    if (is.null(models)) models <- enumerate.models(length(Sgenes),Sgenes,trans.close=control$trans.close,verbose=verbose)
    if(control$pcombi==TRUE){
        models <- models[unlist(lapply(models,function(x){diag(x) <- 0
            is.DAG(x)}))]
    }
    result <- score(models,D,control,verbose)
}

# Bayesian network
else if(inference %in% c("BN.greedy","BN.exhaustive")){
    if(inference == "BN.greedy")
        result = nem.BN(D, "greedy", mode=control$mode, Pm=control$Pm, lambda=control$lambda, verbose=verbose)
    else
        result = nem.BN(D, "exhaustive", mode=control$mode, Pm=control$Pm, lambda=control$lambda, verbose=verbose)
}
else
    stop(paste("Unknown inference method", inference,"\n"))

#------------------------------
# OUTPUT                       
#class(result) <- "nem"
return(result)

}
