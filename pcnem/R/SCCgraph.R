SCCgraph <- function(x,name=TRUE,nlength=20){

	if (!is(x, "graphNEL") & !is(x,"matrix")) stop("Wrong class of argument 'x': must be 'graphNEL' or 'matrix'")  
	if (is(x,"matrix"))  x <- as(x,"graphNEL")
	edgemode(x) = "directed"
	scc   <- strongComp(x)
	N     <- length(scc)
	
	# concatenate node names in same scc
	V <- as.character(1:N)
	if (name==TRUE){
		for (i in 1:N){        
			v <- paste(scc[[i]],collapse=":")
			if (nchar(v)>nlength) v <- paste(substr(v,1,nlength-3),"...",sep="")
			V[i] <- v     
		}
		names(scc) <- V
	}
	
	# which node is in which scc?
	which.scc <- numeric(length(nodes(x)))
	names(which.scc)<-nodes(x)
	for (i in names(scc)) which.scc[scc[[i]]] <- i
	
	
	# build scc graph
	Phi = matrix(0,length(scc), length(scc))
	dimnames(Phi) = list(V, V)
	for (i in names(scc)){		
		Phi[i, which.scc[unlist(adj(x, scc[[i]]))]] = 1
	}
	diag(Phi) = 0
	gR <- as(Phi, "graphNEL")	

	if(is(x,"graphNEL") && (length(edgeDataDefaults(x)) != 0) && (numEdges(gR) > 0)){ # there exist weights in the original graph: transfer them!
		edgeDataDefaults(gR, "label") <- 1
		edgeDataDefaults(gR, "weight") <- 1
		M = as(x, "matrix")
		for(e in 1:length(edges(gR))){
			f = names(edges(gR))[e]
			if(length(edges(gR)[[e]]) > 0){
				for(t in edges(gR)[[e]]){
					MM = M[scc[[f]], scc[[t]],drop=FALSE]
					idx = which.max(abs(MM))
					tmp = MM[idx]
					edgeData(gR, from=f, to=t, attr="weight") = tmp
					edgeData(gR, from=f, to=t, attr="label") = ifelse(abs(tmp) > 1, abs(tmp)-1, abs(tmp))
				}
			}
		}		
	}
		
	return(list(graph=gR,scc=scc,which.scc=which.scc))
}
