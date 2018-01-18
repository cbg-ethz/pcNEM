closest.transitive.greedy = function(Phi, verbose=TRUE){
	prPhi = prune.graph(Phi, verbose=verbose)	
	Phi0 = transitive.closure(prPhi$graph)	
	Phi0 = as(Phi0, "matrix")
	diag(Phi0) = 1
	diag(Phi) = 1
	sco0 = mean(abs(Phi0 - Phi))
	if(verbose)
		cat("initial difference = ", sco0, "\n")
	finished = FALSE
	while(!finished){
		toset = which((Phi-Phi0) == 1)		
		scotmp = double(length(toset))		
		for(t in 1:length(toset)){
			Phi1 = Phi0
			Phi1[toset[t]] = 1
			Phi1 = transitive.closure(Phi1, mat=TRUE, loops=TRUE)			
			scotmp[t] = mean(abs(Phi1 - Phi))			
		}					
		sconew = min(scotmp)
		if(sconew < sco0){
			t = toset[which(scotmp == sconew)]
			Phi0[t] = 1
			Phi0 = transitive.closure(Phi0, mat=TRUE, loops=TRUE)
			sco0 = sconew
		}
		else
			finished = TRUE
		if(verbose)
			cat("difference = ", sco0, "\n")
	}
	Phi0
}
