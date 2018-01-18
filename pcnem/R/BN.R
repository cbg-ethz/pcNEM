#param und score zusammenlegen

createBN = function(data, datatotal=NULL, poss.vertices = NULL,coregraph = NULL, nullnode  = TRUE, marginal = NULL){
	if (!is.character(colnames(data))) {
        colnames(data) = as.character(1:ncol(data))
    }
	if (!is.character(rownames(data))) {
        rownames(data) = as.character(1:nrow(data))
    }
	if(is.null(datatotal)){
		datatotal = matrix(1, nrow = nrow(data), ncol = ncol(data))
	}
	if(is.null(poss.vertices)){
		poss.vertices = 1:length(unique(colnames(data)))
	}
	dimnames(datatotal) = dimnames(data)
	cores = unique(colnames(data))
	if(nullnode){
		cores = c(cores, "null")
		}
	reporters = rownames(data)
	nrrep = length(reporters)
	nrcores = length(cores)
	marginals = rownames(data)
	if(is.null(coregraph)){
		coregraph = matrix(0,ncol = nrcores, nrow = nrcores)
		diag(coregraph) = 1
		}
	if(is.null(marginal)){
		marginal = matrix(0,ncol = nrrep, nrow = nrcores)
	}
	dimnames(coregraph) = list(cores, cores)
	dimnames(marginal) = list(cores, marginals)
	
	var.gr = diag(nrow(coregraph))
	dimnames(var.gr) = dimnames(coregraph)
	if(nullnode){var.gr["null",]=1
		var.gr[,"null"]=1
		coregraph["null","null"]=0
		}
	wh = which(var.gr == 0)
	parameters = matrix(1,nrow=5,ncol=nrrep)
	ve = 1
	k=1
	for(i in wh){
		row = i %% nrow(coregraph)
		if (row==0){row = nrow(coregraph)}
		col = ceiling(i/nrow(coregraph))
		common = (data[,row]%*% data[,col])
		if (common > 0){ve[k]=i
			k = k + 1}}
	BN = list(data = data, datatotal = datatotal, coregraph = coregraph, marginal = marginal, parameters = parameters, 
			poss.vertices = poss.vertices, variable.edges = ve, score = -Inf)
	return(BN)
}

fit.BN = function(BN, method = "greedy", verbose = FALSE, trace = FALSE, mode = "binary_ML", Pm=NULL, lambda=0,...){
	switch(method, exhaustive = {
            return(exhaustive_BN(BN, verbose = verbose, trace = trace, 
            mode = mode, Pm=Pm, lambda=lambda,...))
        }, greedy = {
            return(ingreed_BN(BN, verbose = verbose, trace = trace,
               mode = mode, Pm=Pm, lambda=lambda, ...))
        }, stop("This method does not exist"))
}

# Vor Aufruf der score Funktionen muss gesichert sein, dass der coregraph transitiv abgeschlossen ist und dass die parameters Parameter aktualisiert sind. 

score_discrete_ML = function(BN, loc = NULL, likelihood = TRUE, dummy = TRUE){
		if(!is.null(loc)){ind = loc}
		else{ind = 1:nrow(BN$data)}
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal[,ind]) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		if(is.null(dim(wholegraph))){wholegraph = matrix(wholegraph, ncol = 1)}
		repnr = 0
		for(reporter in ind){
			repnr = repnr + 1
			N1 = BN$datatotal[reporter,] %*% wholegraph[,repnr]
			predicted.pos = N1
			if (predicted.pos == 0){predicted.pos = 1}
			N11 = BN$data[reporter,]%*%wholegraph[,repnr]
			parameters[4,reporter] = max(0.5, N11/predicted.pos)
			parameters[3,reporter] = 1 - parameters[4,reporter]
			N0 = BN$datatotal[reporter, ] %*% (1 - wholegraph[,repnr])
			predicted.neg = N0
			if (predicted.neg == 0){predicted.neg = 1}	
			N00 = ((BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,repnr]))
			parameters[1, reporter] = max(0.5, N00 / predicted.neg)
			parameters[2, reporter] = 1 - parameters[1,reporter]
			score1 = log(parameters[4, reporter] ^ (N11) * (parameters[3, reporter]) ^ (N1 - N11))
			score0 = log(parameters[2, reporter] ^ (N0 - N00) * (parameters[1, reporter]) ^ N00)
			parameters[5,reporter] = score1 + score0
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("p00","p01","p10","p11", "score")
		return(list(parameters = parameters, score = sum(parameters[5,ind])))
}


score_discrete_Bayesian = function(BN, loc = NULL, alpha1 = 5, beta1 = 2, alpha0 = 2, beta0 = 5, likelihood = FALSE, get.parameters = FALSE){
	if(!is.null(loc)){ind = loc}
	else{ind = 1:nrow(BN$data)}
	parameters = BN$parameters
	wholegraph = ((BN$coregraph %*% BN$marginal[,ind]) > 0) * 1
	rownames(wholegraph) = rownames(BN$coregraph)
	wholegraph = wholegraph[colnames(BN$data),]
	if(is.null(dim(wholegraph))){wholegraph = matrix(wholegraph, ncol = 1)}
	repnr = 0
	for(reporter in ind){
		repnr = repnr + 1
		predicted.pos = BN$datatotal[reporter,] %*% wholegraph[,repnr]		
		parameters[4,reporter] = BN$data[reporter,] %*% wholegraph[,repnr]
		parameters[3,reporter] = predicted.pos - parameters[4,reporter]
		predicted.neg = BN$datatotal[reporter,] %*% (1 - wholegraph[,repnr])		
		parameters[1,reporter] = (BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,repnr])
		parameters[2,reporter] = predicted.neg - parameters[1,reporter]
		parameters[5,reporter] = lfactorial(parameters[4, reporter] + alpha1 - 1) + lfactorial(parameters[3, reporter] + beta1 - 1) + lfactorial(parameters[2, reporter] + alpha0 - 1) + lfactorial(parameters[1,reporter] + beta0 - 1) - (lfactorial(predicted.pos + alpha1 + beta1 - 1) + lfactorial(predicted.neg + alpha0 + beta0 - 1))
	}
	if(likelihood){
		parameters[5,] = parameters[5,] + lfactorial(alpha0 + beta0 - 1) + lfactorial(alpha1 + beta1 - 1) - lfactorial(alpha0 - 1)- lfactorial(alpha1 - 1) - lfactorial(beta0 - 1)- lfactorial(beta1 - 1)
		}
	if(get.parameters){
		for(reporter in 1:nrow(BN$data)){
			predicted.pos = BN$datatotal[reporter,] %*% wholegraph[,reporter]		
			parameters[4, reporter] = alpha1 + BN$data[reporter,] %*% wholegraph[,reporter]
			parameters[3, reporter] = beta1 + predicted.pos - parameters[4,reporter]
			predicted.neg = BN$datatotal[reporter,] %*% (1 - wholegraph[,reporter])		
			parameters[1, reporter] = beta0 +  (BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,reporter])
			parameters[2, reporter] = alpha0 + predicted.neg - parameters[1,reporter]
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("beta0", "alpha0", "beta1", "alpha1", "score")
	}
	return(list(score = sum(parameters[5, ind]), parameters = parameters))
}

score_continuous_ML = function(BN, loc = NULL, likelihood = TRUE, dummy = TRUE){
	if(!is.null(loc)){ind = loc}
	else{ind = 1:nrow(BN$data)}
	parameters = BN$parameters
	wholegraph = ((BN$coregraph %*% BN$marginal[,ind]) > 0) * 1
	rownames(wholegraph) = rownames(BN$coregraph)
	wholegraph = wholegraph[colnames(BN$data),]
	if(is.null(dim(wholegraph))){wholegraph = matrix(wholegraph, ncol = 1)}
	repnr = 0
	for(reporter in ind){
		repnr = repnr + 1
		predicted.pos = sum(wholegraph[,repnr])
		if (predicted.pos == 0){predicted.pos = 1}
		parameters[3,reporter] = (BN$data[reporter,] %*% wholegraph[,repnr])/predicted.pos
		parameters[4,reporter] = sum(((BN$data[reporter,] - parameters[3,reporter])^2) %*% wholegraph[,repnr]) / predicted.pos
		predicted.neg = sum(1 - wholegraph[,repnr])
		if(predicted.neg == 0){predicted.neg = 1}
		parameters[1,reporter] = (BN$data[reporter,]%*%(1-wholegraph[,repnr]))/predicted.neg
		parameters[2,reporter] = sum(((BN$data[reporter,] - parameters[1,reporter])^2) %*% (1 - wholegraph[,repnr])) / predicted.neg
		score1 = sum(dnorm(BN$data[reporter,which(wholegraph[,repnr] == 1)], parameters[3,reporter], sqrt(parameters[4,reporter]), log = TRUE))
		score0 = sum(dnorm(BN$data[reporter,which(wholegraph[,repnr] == 0)], parameters[1,reporter], sqrt(parameters[2,reporter]), log = TRUE))
		parameters[5, reporter] = score1 + score0
	}
	colnames(parameters) = rownames(BN$data)
	rownames(parameters) = c("mu0", "var0", "mu1", "var1", "score")
	return(list(parameters = parameters, score = sum(parameters[5,ind])))
}


score_continuous_Bayesian = function(BN, loc = NULL, mu0 = 10, mu1 = 10, alpha1 = 10, beta1 = 1, alpha0 = 10, beta0 =1, v0 = 0.001, v1 = 0.001, likelihood = FALSE, get.parameters = FALSE){
	if(!is.null(loc)){ind = loc}
	else{ind = 1:nrow(BN$data)}
	parameters = BN$parameters
	wholegraph = ((BN$coregraph %*% BN$marginal[,ind]) > 0) * 1
	rownames(wholegraph) = rownames(BN$coregraph)
	wholegraph = wholegraph[colnames(BN$data),]
	if(is.null(dim(wholegraph))){wholegraph = matrix(wholegraph, ncol = 1)}
	repnr = 0
	for(reporter in ind){
		repnr = repnr + 1
		N1 = sum(wholegraph[,repnr])
		div1 = N1
		if(div1 == 0){div1 = 1}	
		parameters[3,reporter] = (BN$data[reporter,] %*% wholegraph[,repnr])/div1
		parameters[4,reporter] = sum(((BN$data[reporter,] - parameters[3,reporter])^2) %*% wholegraph[,repnr])
		N0 = sum(1 - wholegraph[,repnr])
		div0 = N0
		if(div0 == 0){div0 = 1}
		parameters[1,reporter] = (BN$data[reporter,]%*%(1-wholegraph[,repnr]))/div0
		parameters[2,reporter] = sum(((BN$data[reporter,] - parameters[1,reporter])^2) %*% (1 - wholegraph[,repnr]))
		score1 = (v1 / (v1 + N1)) ^ (1 / 2)  * gamma((alpha1 + N1)/2) /  abs(beta1 + parameters[4,reporter] + v1 * N1 / (v1 + N1) * (parameters[3,reporter] - mu1) ^ 2) ^ ((alpha1 + mu1) / 2)
		score2 = (v0 / (v0 + N0)) ^ (1 / 2)  * gamma((alpha0 + N0)/2) /  abs(beta0 + parameters[2,reporter] + v0 * N0 / (v0 + N0) * (parameters[1,reporter] - mu0) ^ 2) ^ ((alpha0 + mu0) / 2)
		parameters[5,reporter] = score1 * score2
	}
	if(likelihood){
		parameters[5,reporter] = parameters[5,reporter] * 1/pi^(ncol(BN$data) / 2) * abs(beta0) ^ (alpha0 / 2) * abs(beta1) ^ (alpha1 / 2) / (gamma(alpha1/2) * gamma(alpha0/2))
		}
	parameters[5,reporter] = log(parameters[5,reporter])
	if(get.parameters){
		parameters = matrix(1, nrow = 9, ncol = nrow(BN$data))
		parameters[9, ] = parameters[5, ]
		for(reporter in 1:nrow(BN$data)){
			N1 = sum(wholegraph[,reporter])
			div1 = N1
			if(div1 == 0){div1 = 1}
			mean1 = (BN$data[reporter,] %*% wholegraph[,reporter])/div1
			s1 = ((BN$data[reporter,] - mean1)^2) %*% wholegraph[,reporter]
			parameters[5, reporter] = alpha1 + N1
			parameters[6, reporter] = beta1 + s1 + v1*N1/(v1+N1) * (mean1 - mu1)^2
			parameters[7, reporter] = (v1*mu1 + N1*mean1)/(v1+N1)
			parameters[8, reporter] = v1 + N1
			N0 = sum(1 - wholegraph[,reporter])
			div0 = N0
			if(div0 == 0){div0 = 1}
			mean0 = (BN$data[reporter,]%*%(1-wholegraph[,reporter]))/div0
			s0 = ((BN$data[reporter,] - mean0)^2) %*% (1 - wholegraph[,reporter])
			parameters[1, reporter] = alpha0 + N0
			parameters[2, reporter] = beta0 + s0 + v0 * N0 / (v0 + N0) * (mean0 - mu0)^2
			parameters[3, reporter] = (v0 * mu0 + N0 * mean0)/(v0 + N0)
			parameters[4, reporter] = v0 + N0
			
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("alpha0", "beta0", "mu0", "v0", "alpha1", "beta1", "mu1", "v1", "score")		
	}
	return(list(parameters = parameters, score = sum(parameters[5,ind])))
}

score_BN = function(BN, mode = "binary_ML", loc = NULL, Pm=NULL, lambda=0, ...){
	if(mode == "binary_ML"){sco = (score_discrete_ML(BN, loc))}
	else if(mode == "continuous_ML"){sco = (score_continuous_ML(BN, loc))}
	else if(mode == "binary_Bayesian"){sco = (score_discrete_Bayesian(BN, loc, ...))}
	else if(mode == "continuous_Bayesian"){sco = (score_continuous_Bayesian(BN, loc, ...))}
	else{stop("inexistant mode")}
	Phi = BN$coregraph
	if(!is.null(Pm)){
		if(lambda != 0)
  			sco <- sco - lambda*sum(abs(Phi - Pm))
		lpPhi <- sapply(Phi, PhiDistr, Pm, a=1, b=0.5)
		sco = sco + sco + lpPhi
	}
	sco
}


exhaustive_BN = function (BN, models = NULL, verbose = TRUE, trace = FALSE, mode="binary_ML", IC.mode = "Bayesian", nu = c(0.001,0.1, 1,Inf), Pm=NULL, lambda=0,...){	
    core = BN$coregraph
    BN$coregraph = transitive.closure(BN$coregraph, mat=TRUE, loops=TRUE)
    diag(BN$coregraph) = diag(core)
    if(is.null(models)){models = enumerate.models2(BN$coregraph, BN$variable.edges, acyclic = TRUE, name = colnames(BN$coregraph))}
    BN$coregraph = models[[1]]
    lastcoregraph = models[[1]]
    marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)
    BN$marginal = marginal$graph
    score = score_BN(BN, mode = mode, ...)
    bestpost = score$score
    BN$parameters = score$parameters
    nrofedges = ncol(BN$coregraph)
	if(trace){post = rep(-Inf, 2^nrofedges)
		post[1] = bestpost
	}
    if(verbose){k = 1}
	if(verbose){cat("score: ", bestpost," ")}
	bestBN = BN
	for (j in 2:length(models)) {
		if(verbose){k = k + 1}
		lastcoregraph = BN$coregraph
		BN$coregraph = models[[j]]
		marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, lastcoregraph = lastcoregraph, lastscore = marginal$lastscore, ...)
		BN$marginal = marginal$graph
		score = score_BN(BN, mode = mode, ...)
    	BN$score = score$score
    	BN$parameters = score$parameters
		if(verbose){
    		cat(" nr: ", k, "score: ", BN$score, " ")
        	if(k %% 5 == 0){cat("\n")}}
        if (BN$score > bestpost) {
            bestpost = BN$score
            bestBN   = BN
        }
        if(trace){post[j + 1] = BN$score}
		
    }
    score = score_BN(bestBN, mode = mode, get.parameters = TRUE, ...)
    bestBN$score = score$score
    bestBN$parameters = score$parameters	
    if (trace)
        return(list(bestBN,post))
    return(bestBN)
}

optimizemarginal = function(BN, mode,lastcoregraph = NULL, lastscore = NULL, nu = c(0,0.001,1,10,Inf), trace = FALSE, verbose = FALSE, IC.mode = "Bayesian", ...){
	list = as.list(1:length(nu))
	IC = numeric(length(nu))
	nr=1
	if(is.null(lastscore)){sc = as.list(nu)}
	else {sc = lastscore}
	pennr = 0
	for(penalty in nu){
		pennr = pennr + 1
		if(penalty == Inf){marginal = optimizemarginal.nuInf(BN, mode, lastcoregraph = lastcoregraph, lastscore = sc[[pennr]],...)
			BN$marginal = marginal$graph
			sc[[pennr]] = marginal$lastscore}
		else{
			if(is.null(lastcoregraph)){sc[[pennr]] = matrix(-Inf, nrow = nrow(BN$marginal), ncol = ncol(BN$marginal))
				change.sc = BN$poss.vertices}
			else{sc[[pennr]] = lastscore[[pennr]]
				change.sc = which(apply(BN$coregraph != lastcoregraph, 2, sum) != 0)}
			
			for (rep in 1:ncol(BN$marginal)){
				pr = numeric(4)
				BN$marginal[,rep]=0
				score.rep2 = -Inf
				continue = TRUE
				k = 1
				while(continue){
					if(k > 1){sc2 = rep(-Inf, nrow(BN$marginal))}
					score.rep1 = score.rep2
					if(k == 1){poss.vertices = intersect(change.sc,BN$poss.vertices[which(BN$marginal[BN$poss.vertices,rep] < 1)])}
					else{poss.vertices = BN$poss.vertices[which(BN$marginal[BN$poss.vertices,rep] < 1)]}
					for (sign in poss.vertices){
						BN$marginal[sign,rep]=1	
						score_para = score_BN(BN, mode = mode, loc = rep, ...)
    					score = score_para$score
				    	if(k == 1){sc[[pennr]][sign,rep] = score}
						else{sc2[sign] = score}
						BN$marginal[sign,rep] = 0
					}
					if(k == 1){
						score.rep2 = max(sc[[pennr]][,rep])
						new.edges = which(sc[[pennr]][,rep] == max(sc[[pennr]][,rep]))}
					else{score.rep2 = max(sc2)
						new.edges = which(sc2 == max(sc2))}
					BN$marginal[new.edges, rep] = 1
					if(length(new.edges > 1)){
						score = score_BN(BN, mode = mode, loc = rep)
    					score.rep2 = score$score
				    	BN$parameters = score$parameters}					
					if(k == 1){continue = (score.rep2 - score.rep1 > (length(new.edges) - 1) * penalty)}
					else{continue = (score.rep2 - score.rep1 > length(new.edges) * penalty)}
					if(!continue){BN$marginal[new.edges,rep] = 0
						if(k == 1){BN$marginal[which.max(sc[[pennr]][,rep]), rep] = 1}
					}
					k = k + 1
				}
			}
		}
		list[[nr]] = BN$marginal
		score = score_BN(BN, mode = mode, likelihood = TRUE, ...)
		likelihood = score$score
    	BN$parameters = score$parameters
		if(IC.mode == "Bayesian"){	
			IC[nr] =  2 * (likelihood) - log(sum(BN$datatotal)) * sum(BN$marginal)}
		else if(IC.mode == "Akaike"){	
			IC[nr] =  2 * (likelihood) - 2 * sum(BN$marginal)}
		else (stop("inexistand IC.mode"))
			nr = nr + 1
	}
	if (verbose)
		{cat("nu.opt = ",nu[which.max(IC)]," nr marg not one connection = ", sum(apply(list[[which.max(IC)]],2,sum)!=1), " ")}
	if(trace){return(list(graph = list[[which.max(IC)]], list = list,IC = IC,lastscore = sc))}
	return(list(graph = list[[which.max(IC)]],lastscore = sc))
}


optimizemarginal.nuInf = function(BN, mode, lastcoregraph = NULL, lastscore = NULL,...){
	Daten = BN$Daten
	if(is.null(lastcoregraph)){sc = matrix(-Inf, nrow = nrow(BN$marginal), ncol = ncol(BN$marginal))
		change.sc = BN$poss.vertices}
	else{sc = lastscore
		change.sc = which(apply(BN$coregraph != lastcoregraph, 2, sum) != 0)}
	
	for (rep in 1:ncol(BN$marginal)){
		pr = rep(1, times = 4)
		Daten = BN$Daten
		BN$marginal[,rep]=0
		for (sign in change.sc){
			BN$marginal[sign,rep]=1	
			score_para = score_BN(BN, mode = mode, loc = rep, ...)
			score = score_para$score
			sc[sign, rep] = score
			BN$marginal[sign,rep]=0
			}
		BN$marginal[which.max(sc[,rep]),rep] = 1
		}
	graph = BN$marginal
	return(list(graph = graph, lastscore = sc))
	}


ingreed_BN = function(BN, verbose = FALSE, mode = "binary_ML", trace = TRUE, max.iter = 100, IC.mode = "Bayesian", nu = c(0,0.001,1,0.1,Inf), Pm=NULL, lambda=0,...) 
{
    if (trace) {
        tracelist = list()
    }
    if (verbose) 
        cat("Number of iterations : ")
    counter = 0
    marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)
    BN$marginal = marginal$graph
    repeat {
        coregraphold = BN$coregraph
        BN$coregraph = optimizecoregraph(BN, mode, ...)
        if (identical(BN$coregraph, coregraphold))
            {break()}
        marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, lastcoregraph = coregraphold, lastscore = marginal$lastscore, ...)
		BN$marginal = marginal$graph
        counter = counter + 1
        if (verbose){
            cat("nr of iterations ", counter, " ")
            if ((counter%%5) == 0) 
                cat("\n")
        }
        if (trace)
            tracelist[[counter]] = BN
        if (counter == max.iter) 
            break()
    }
    score = score_BN(BN, mode = mode, ...)
    BN$score = score$score
    BN$parameters = score$parameters
	if (verbose) 
        cat("\n")
    if (!trace) 
        return(BN)
    return(list(BN = BN, tracelist = tracelist))
}

optimizecoregraph = function(BN, mode, best.marg=FALSE, verbose=FALSE, IC.mode = "Bayesian", nu = c(0,0.001,1,0.1,Inf), ...){
	bestcoregraph = BN$coregraph
	original = BN$coregraph
	original.red = transitive.reduction(original)
	lastcoregraph = NULL
	score = score_BN(BN, mode = mode, ...)
    bestpost = score$score
    BN$parameters = score$parameters
	for (i in BN$variable.edges){
	    core = original.red
	    core[i] = 1 - core[i]
		BN$coregraph = transitive.closure(core, mat=TRUE, loops=TRUE)
		diag(BN$coregraph) = diag(original)
		if(!is.dag(core) || identical(original,BN$coregraph)){next()}
		if(best.marg){BN$marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)}
		score = score_BN(BN, mode = mode, ...)
    	post = score$score
    	BN$parameters = score$parameters
		lastcoregraph = BN$coregraph
    	if (post > bestpost){
			bestcoregraph = BN$coregraph
			bestpost=post}
		}
    return(bestcoregraph)
}

is.dag=function(graph){
	diag(graph) = 0
	graph1 = graph
	for (i in 1:nrow(graph))
		{graph1 = graph1 %*% graph
		if(sum(graph1)==0){return(TRUE)}
	}
	return(FALSE)
}

# simulate_BN = function(vert = 4, edges = NULL, gamma = 1, nr_intven = 3, beta = 0.9, reporters = 10, nullnode = TRUE, maxOutDegree = vert, maxInDegree = vert){
# 	core = diag(vert)
# 	colnames(core) = as.character(1:vert)
# 	rownames(core) = as.character(1:vert)	
# 	if(is.null(edges)){
# 		for(i in 1:vert){ # connect network randomly
# 			maxOutDegree = min( i,maxOutDegree)	
# 			degprob = (0:maxOutDegree)^(-gamma)
# 			degprob[1] = 1
# 			degprob = degprob/sum(degprob)
# 			outdeg = sample(0:maxOutDegree, 1, prob = degprob)# power law for out-degree => scale-free network
# 			if(outdeg > 0){
# 				idx0 = which(core[i,] == 0 & 1:vert < i)
# 				if(length(idx0) > 1){
# 					idx = which(colSums(core[,idx0]) <= maxInDegree)
# 					if(length(idx) > 0){			
# 						idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
# 						core[i,idx] = 1					
# 					}			
# 				}
# 				else if(length(idx0) == 1){
# 					idx = which(sum(core[idx0]) <= maxInDegree)
# 					if(length(idx) > 0){			
# 						idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
# 						core[i,idx] = 1					
# 					}			
# 				}	
# 			}
# 		}							
# 	}
# 	else{
# 		ve = which(core ==0)
# 		veind = which(core ==0, arr.ind = TRUE)
# 		veinddiff = veind[,1]-veind[,2]
# 		wh = which(veinddiff > 0)
# 		wh = ve[wh]
# 		sample = sample(wh,size = edges, replace = FALSE)
# 		core[sample]=1
# 	}
# 	original = core
#     core = transitive.closure(core, mat=TRUE, loops=TRUE)
#     diag(core) = diag(original)
# 	intven = rep(1:vert, each = nr_intven)	
# 	Delta = t(core)[,intven]
# 	ind = sample(1:nrow(core), reporters, replace = TRUE)
# 	Delta = Delta[ind,]
# 	colnames(Delta) = paste("I", intven,sep=".")
# 	if(nullnode){
# 		l = nrow(core)
# 		coregraph = matrix(0,nrow = l+1, ncol = l+1)
# 		coregraph[1:l,1:l] = core
# 		rownames(coregraph) = c(unique(colnames(Delta)),"null")
# 		colnames(coregraph) = c(unique(colnames(Delta)),"null")
# 		core = coregraph}
# 	Theta = matrix(0, ncol=reporters, nrow =ncol(core))
# 	for(i in 1:ncol(Theta))
# 		{Theta[ind[i],i] = 1}
# 	if(nullnode){rownames(Theta) = c(unique(colnames(Delta)),"null")}
# 	else{rownames(Theta) = unique(colnames(Delta))} 
# 	colnames(Theta) = as.character(1:ncol(Theta))
# 	report = rbinom(length(Delta), size = 1, prob = beta)
# 	report[which(Delta == 0)] = 1 - report[which(Delta == 0)]
# 	dim(report) = dim(Delta)
# 	colnames(report)=colnames(Delta)
# 	rownames(report) = colnames(Theta)
# 	BN = createBN(data = report, coregraph = core, marginal = Theta, nullnode = nullnode)
# 	return(BN)
# }
# 
# plotcoregraph = function (coregraph, main = NULL, transitive.reduction = TRUE) 
# {
#     gr = coregraph
#     diag(gr) = 0
#     wh1 = which(apply(gr,1,sum) > 0)
#     wh2 = which(apply(gr,2,sum) > 0)
#     wh = union(wh1,wh2)
#     if (length(wh)<1)(stop("coregraph empty"))
#     gr = gr[wh,wh]
#     if(transitive.reduction){gr = transitive.reduction(gr)}
#     grph = as(gr, "graphNEL")
#     plot(grph, main = main)
# }


enumerate.models2 = function (model, variable.edges, name = NULL, verbose = TRUE, acyclic = FALSE) 
{
    bc <- bincombinations(length(variable.edges))
    fkt1 <- function(x, name) {
        M <- model
        M[variable.edges] <- x
        dimnames(M) <- list(name, name)
        M <- transitive.closure(M, mat = TRUE, loops = TRUE)
        diag(M) = diag(model)
        return(list(M))
    }
    models <- apply(bc, 1, fkt1, name)
    n = ncol (model)
    if(acyclic){
	    models = lapply(models, erase.cycles)
	    }
    models <- unique(matrix(unlist(models), ncol = n * n, byrow = TRUE))
    fkt2 <- function(x, n, name) {
        M <- matrix(x, n)
        dimnames(M) <- list(name, name)
        return(list(M))
    }
    models <- unlist(apply(models, 1, fkt2, n, name), recursive = FALSE)
    if (verbose) 
        cat("Generated", length(models), "unique models ( out of", 
            2^(length(variable.edges)), ")\n")
    return(models)
}

erase.cycles = function(model){
	w = which(model[[1]]==1, arr.ind = TRUE)
	w1 = which(w[,1]-w[,2]!=0)
	w2 =w[w1,]
	if(is.null(dim(w2))){return(model)}
	all = rbind(cbind(w2[,2],w2[,1]), w2)
	x = length(all)
	y = length(unique(all))
	if(x>y){model[[1]][,] = diag(nrow(model[[1]]))}
	return(model)
}
