nem.BN = function(D, inference="greedy", mode="binary_Bayesian", Pm=NULL, lambda=0, verbose=TRUE){
	if(!(inference %in% c("greedy", "exhaustive"))) stop("Inference method has to be either 'greedy' or 'exhaustive'")
	if(!(mode %in% c("binary_ML", "continuous_ML","binary_Bayesian","continuous_Bayesian"))) stop(paste("Unknown mode '", mode, "'\n"))
	Sgenes = unique(colnames(D))
	BN = createBN(D, nullnode=FALSE)	
	best = fit.BN(BN, method=inference, mode=mode, Pm=Pm, lambda=lambda, trace=TRUE)				
	diag(best$BN$coregraph) = 0
	graph = as(best$BN$coregraph, "graphNEL")		
	mLL = best$BN$score
	mappos = apply(best$BN$marginal>0,1,which)
	if(!is.null(rownames(D)))
		mappos = sapply(mappos, names) 		
	selected = unique(unlist(mappos[Sgenes]))
	res = list(graph=graph, mLL=mLL, mappos=mappos, selected=selected, lam=lambda, type=mode)
	class(res) <- "nem.BN"
	if(verbose)
		cat("log-likelihood of model = ",res$mLL,"\n")
	return(res)
}
