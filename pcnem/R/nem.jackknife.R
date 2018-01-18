nem.jackknife <- function(D, thresh=0.5, inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE){
	if(inference == "dynoNEM")
		stop("nem.jackknife is not applicable for dynoNEMs")
	
	inferNetwork <- function(sgenesk){	
		if(verbose)
			cat("S-genes used: ", sgenesk,"\n")
		Dtmp = D[,sapply(sgenesk,function(x) which(colnames(D) %in% x))]
		controltmp = control
		controltmp$Pe = control$Pe[,sgenesk]	
		controltmp$Pm = control$Pm[sgenesk,sgenesk]
		if(!is.null(control$Pm) & length(control$lambda) > 1)
			res = as(nemModelSelection(control$lambda,Dtmp,inference,models,controltmp,verbose)$graph,"matrix")
		else
			res = as(nem(Dtmp,inference,models,controltmp,verbose)$graph,"matrix")
		overlapJack[sgenesk,sgenesk] = overlapJack[sgenesk,sgenesk] + res
		overlapJack
	}
	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
	n = length(Sgenes)		
	overlapJack = matrix(0,ncol=n,nrow=n)	
	colnames(overlapJack) = Sgenes
	rownames(overlapJack) = Sgenes	
	for(i in 1:n)
		overlapJack = inferNetwork(Sgenes[-i])
	overlapJack = overlapJack / (n-2)	
	overlapJack = round(overlapJack,digits=2)			
	res = nem(D,models=list((overlapJack>thresh)*1),inference="search", control, verbose=verbose)	
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	res$para = res$para[[1]]
	res$control= control	
	g = res$graph
	edgeDataDefaults(g, "label") = 1	
	edgeDataDefaults(g, "weight") = 1
	for(s1 in Sgenes){
		for(s2 in Sgenes){
			if(s2 %in% unlist(adj(g, s1))){
				edgeData(g, from = s1, to = s2, attr = "weight") = overlapJack[s1,s2]			
				edgeData(g, from = s1, to = s2, attr = "label") = overlapJack[s1,s2]
			}
		}
	}
	res$graph = g
	class(res) <- "nem.jackknife"
	res
}
