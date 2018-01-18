nem.featureselection <- function(D,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))),verbose=TRUE, tol=1e-4){
		
	control$selEGenes=FALSE	
	if(control$selEGenes.method == "regularization"){
		Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
		nrS = length(Sgenes)
		 if (is.null(control$Pe)){ 	
			control$Pe <- matrix(1/nrS,nrow=nrow(D),ncol=nrS)
			colnames(control$Pe) <- Sgenes  		
  		}    		
		deltaseq = 0:10
		if ("doMC" %in% loadedNamespaces()){
			registerDoMC(control$mc.cores)
			results = foreach(d = deltaseq) %dopar% {					
				net = nem(D, inference, models,control,verbose=verbose)	
				if(verbose){
					cat(length(net$selected), " selected E-genes (delta = ", d, ", BIC = ", -2*net$mLL + log(nrow(D))*length(net$selected), "):", sort(net$selected)[1:min(20, length(net$selected))]," ...\n")				
				}			
				net
			}
		}
		else{
			results = sapply(deltaseq, function(d){					
						net = nem(D, inference, models,control,verbose=verbose)	
						if(verbose){
							cat(length(net$selected), " selected E-genes (delta = ", d, ", BIC = ", -2*net$mLL + log(nrow(D))*length(net$selected), "):", sort(net$selected)[1:min(20, length(net$selected))]," ...\n")				
						}			
						net
					})	
		}		    
 		s = -2*unlist(results["mLL",]) + log(nrow(D))*sapply(results["selected",], length) # BIC model selection
# 		s = -unlist(results["mLL",])/sapply(results["selected",], length) # Achim's original approach		
 		winner = which.min(s)				
		net = results[,winner]	
		class(net) = inference
	}
	else{
		if(inference %in% c("dynoNEM", "mc.eminem"))
			stop("This feature selection method is not applicable to 'dynoNEM' and 'mc.eminem'")
		converged = FALSE
		while(!converged){
			net = nem(D, inference, models, control,verbose=verbose)
			control$delta = 1e-10
			sel = getRelevantEGenes(as(net$graph,"matrix"), D, control)			
			if(verbose)
				cat("Old Likelihood = ", net$mLL, "new likelihood = ", sel$mLL, "\n", length(sel$selected), " selected E-genes:", sort(sel$selected)[1:min(20, length(sel$selected))]," ...\n")			
			converged = (abs(max(net$mLL) - max(sel$mLL))/abs(max(net$mLL)) < tol)
			if(converged)
				break
			D = D[sel$selected,]
		}
		net$selected = sel$selected
		net$mLL = sel$mLL		
	}
	if(verbose)
		cat("Converged! ", length(net$selected),"selected E-genes\n")
	net
}
