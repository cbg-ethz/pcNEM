prune.graph <- function(g,cutIN=NULL,cutOUT=NULL,quant=.95,verbose=TRUE){

if (class(g)=="matrix") g <- as(g,"graphNEL")
# compute degree of nodes
edgemode(g) = "directed"
dg <- degree(g)
gc = transitive.closure(as(g,"matrix"), loops=FALSE)
edgemode(gc) = "directed"
dc = degree(gc)
nr = 0
# compute missing edges in and out
miss.in  <- dc$inDegree  - dg$inDegree
miss.out <- dc$outDegree - dg$outDegree
if(numEdges(gc) > 0){
	if (is.null(cutIN )) cutIN  <- ceiling(quantile(miss.in ,quant))
	if (is.null(cutOUT)) cutOUT <- ceiling(quantile(miss.out,quant))	
		
	if (verbose) cat("cutIN:",cutIN," -  cutOUT:",cutOUT,"\n")
	if (cutIN==0 | cutOUT==0) warning("ALL edges removed since at least one cutoff is 0")
	
	# remove in-edges
	#----------------------------
	iE <- inEdges(g)
	nrIN <- 0
	removeIN <- nodes(g)[which(miss.in >= cutIN)]
	# make sure there are edges
	removeIN <- removeIN[!unlist(lapply(iE[removeIN],function(x) length(x)==0))]
	for(i in removeIN){ 
		g <- removeEdge(iE[[i]],i,g)
		nrIN <- nrIN + length(iE[[i]])
	}
	
	# remove out-edges
	#----------------------------
	E <- edges(g)
	nrOUT <- 0
	removeOUT <- nodes(g)[which(miss.out >= cutOUT)]
	# make sure that edges were not removed in last step
	removeOUT <- removeOUT[!unlist(lapply(E[removeOUT],function(x) length(x)==0))]
	for(j in removeOUT){ 
		g <- removeEdge(j,E[[j]],g)
		nrOUT <- nrOUT + length(E[[i]])
	}
	
	# output
	nr <- nrIN + nrOUT
	if (verbose) if (nr==1) cat("Removed 1 edge\n") else cat("Removed",nr,"edges\n")
}
return(list(graph=g,removed=nr,missing.in=miss.in,missing.out=miss.out))
}
