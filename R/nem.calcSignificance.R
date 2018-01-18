nem.calcSignificance <- function(D, x, N=1000, which.test=c("permutation", "rand.net", "rand.modify"), seed=1, mc.cores=8){	
	modify.rand = function(Phi0){	
		if(x$control$trans.close){
			Phired = transitive.reduction(Phi0)
			actions = c("insertion","deletion")				
		}
		else{
			Phired = Phi
			if(class(x) == "dynoNEM")
				actions = c("weight increase", "weight decrease", "weight swap")
			else
				actions = c("insertion","deletion","reversion")				
		}
		if(all(Phi0 - diag(ncol(Phi0)) == 0))
			actions = "insertion"	
		else if(!any(Phi0 == 0))
			actions = setdiff(actions,"insertion")
		if(length(actions) > 1)
			act = sample(actions,1)	
		else
			act = actions[1]
		if(act == "insertion"){		
			Phinew = Phi0
			idx = sample(which(Phi0 == 0),1)
			Phinew[idx] = 1
			if(x$control$trans.close)		
				Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE)					
		}
		else if(act == "deletion"){
			Phinew = Phired
			diag(Phired)=0
			idx = sample(which(Phired == 1),1)
			Phinew[idx] = 0
			if(x$control$trans.close)		
				Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE)			
		}	
		else if(act == "reversion"){
			Phinew = Phi0
			allidx = which((Phi0 + t(Phi0) == 1) & Phi0, arr.ind=TRUE)
			idx = sample(1:NROW(allidx),1)
			Phinew[allidx[idx,1],allidx[idx,2]] = 0
			Phinew[allidx[idx,2],allidx[idx,1]] = 1
		}
		else if(act == "weight increase"){
			Phinew = Phi0
			idx = sample(which(Phi0 < dim(D)[1]), 1)
			Phinew[idx] = Phinew[idx] + 1
		}
		else if(act == "weight decrease"){
			Phinew = Phi0
			idx = sample(which(Phi0 > 0), 1)
			Phinew[idx] = Phinew[idx] - 1
		}
		else if(act == "weight swap"){
			Phinew = Phi0
			allidx = which(Phi0 != t(Phi0), arr.ind=TRUE)
			idx = sample(1:NROW(allidx), 1)
			Phinew[allidx[idx,1], allidx[idx,2]] = Phi0[allidx[idx,2], allidx[idx,1]]	
			Phinew[allidx[idx,2], allidx[idx,1]] = Phi0[allidx[idx,1], allidx[idx,2]]
		}
		return(Phinew)
	}	

	conduct.test = function(which.test){
		if(which.test == "rand.net"){
			if(class(x) == "dynoNEM")
				stop("test w.r.t. random network is not yet available for dynoNEMs")
			net = sampleRndNetwork(Sgenes)		
			loglik <- nem(D,models=list(net),inference="search", x$control, verbose=FALSE)$mLL
		}
		else if(which.test == "rand.modify"){
			net = modify.rand(Phi)		
			if(class(x) != "dynoNEM")
				loglik <- nem(D,models=list(net),inference="search",x$control, verbose=FALSE)$mLL
			else{
				nreps = sapply(x$control$Sgenes, function(s) sum(dimnames(D)[[3]] == s))
				loglik = dynoNEM.posteriorEGenePos(net, D, priorE=x$control$Pe, delta=x$control$delta, type=x$control$type, nrep=nreps[1], alpha=x$control$para[1], beta=x$control$para[2])$mLL
			}
		}
		else if(which.test == "permutation"){
			rnd = sample(1:ncol(Phi), replace=FALSE)
			net = Phi[rnd,rnd]		
			dimnames(net) = list(Sgenes, Sgenes)
			if(class(x) != "dynoNEM")		
				loglik <- nem(D,models=list(net),inference="search",x$control, verbose=FALSE)$mLL
			else{
				nreps = sapply(x$control$Sgenes, function(s) sum(dimnames(D)[[3]] == s))
				loglik = dynoNEM.posteriorEGenePos(net, D, priorE=x$control$Pe, delta=x$control$delta, type=x$control$type, nrep=nreps[1], alpha=x$control$para[1], beta=x$control$para[2])$mLL
			}
		}
		return(likelihood - loglik)
	}
	
	which.test = match.arg(which.test, c("permutation", "rand.net", "rand.modify"), several.ok=FALSE)
	cat("Testing significance of NEM model\nTest: ", which.test, "\n")
	likelihood = x$mLL
	Phi = as(x$graph, "matrix")
	if(class(x) != "dynoNEM")
		Phi[Phi != 0] = 1	
	Sgenes = colnames(Phi)
	set.seed(seed)
	if(class(x) == "dynoNEM")
		x$control$trans.close = FALSE
				
	if ("doMC" %in% loadedNamespaces()){
		registerDoMC(mc.cores)
		test = foreach(i=1:N)%dopar% {		
			cat(".")
			BF = conduct.test(which.test)			
		}
		BF = unlist(test$BF)
	}
	else{
		BF = double(N)
		for(i in 1:N){
			cat(".")
			BF[i] = conduct.test(which.test)
		}		
	}
	p.value = permp(sum(BF <= 0), N, twosided=FALSE, total.nperm=N) 
	cat("\ndone.\n")
	
	p.value
}

