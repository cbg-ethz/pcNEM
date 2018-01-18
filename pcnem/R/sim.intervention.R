sim.intervention = function(x, int, T=1){	
	M = as(x$graph, "matrix")
	if(length(x$mLL) > 1 & class(x) != "dynoNEM"){
		winner <- which.max(x$mLL)
		pos <- x$mappos[[winner]]
	}
	else{
		pos <- x$mappos
	}		  
	diag(M) = 1	 
  if(class(x) != "dynoNEM"){
  	Sgenes.effected = lapply(int, function(i) colnames(M)[M[i,] == 1])		
  }
  else{
  	if(is.na(T))
  		stop("Number of time steps has to be provided for dynoNEM simulation")
  	Sgenes.effected = lapply(int, function(i) apply(dynoNEM.perturb(M, T, i), 2, function(x) names(x)[x==1]))		
  }	
  Egenes.effected = lapply(Sgenes.effected, function(i) unique(unlist(pos[i])))
  Egenes.effected.prob = lapply(Sgenes.effected, function(S){EP=rowSums(x$pos[, S]); names(EP) = rownames(x$pos); EP})  
  names(Egenes.effected) = names(Egenes.effected.prob) = names(Sgenes.effected) = paste("intervention", 1:length(int))
	list(Sgenes.effected=Sgenes.effected, Egenes.effected=Egenes.effected, Egenes.effected.prob)
}
