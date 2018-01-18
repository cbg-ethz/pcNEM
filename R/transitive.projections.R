# library(Rgraphviz)
# library(e1071)
# library(caTools)
# library(combinat)

## The following function computes the transitive projections (approximations) of the input graph. Input to the function is either a graphNEL object or the adjacency matrix

transitive.projections <- function(adjmat)
{
   if (!(class(adjmat)%in%c("graphNEL","matrix")))
   stop("Input must be either graphNEL object or adjacency matrix")

## if input is an adjacency matrix
   if (class(adjmat)=="matrix")
    {
    n <- ncol(adjmat)
    }

## If input is a graph
## convert the graph into its adjacency matrix
    graph2adj <- function(gR) {
#        adj.matrix <- matrix(0,
#                        length(nodes(gR)),
#                        length(nodes(gR))
#                        )
#        rownames(adj.matrix) <- nodes(gR)
#        colnames(adj.matrix) <- nodes(gR)
#        for (i in 1:length(nodes(gR))) {
#        adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
#         }
#        return(adj.matrix)
	as(gR, "matrix")
     }

   if (class(adjmat)=="graphNEL")
    {
      adjmat <- graph2adj(adjmat)
      n <- ncol(adjmat)
    }

## If input graph is transitive, then return the same graph with minimum distance 0
   if (is.transitive(adjmat,nrow(adjmat)))
    {
    return(list(Optimum=adjmat, MinDist=0))
    }

## To speed up the process one could start with some transitive subgraph of the initial graph and proceed. The following are the initialization steps. In case no transitive subgraphs are found, then we start with the empty graph.
        imat <- transSubGr(adjmat)
	s0 <- list(imat) # list of transitive approximations
	s1 <- OneNeighborhood(imat) # suboptimal transapprox of dist 1
	s2 <- TwoNeighborhood(imat) # suboptimal transapprox of dist 2
	s3 <- ThreeNeighborhood(imat) # suboptimal transapproxs of dist 3
	s4 <- FourNeighborhood(imat)  # suboptimal transapproxs of dist 4
	edgesofg <- as.vector(adjmat) # edges of the givengraph
	currentgraph <- as.vector(matrix(imat,nrow=n,ncol=n))
	mindist = 0
	NrRemainingEdges <- setdiff(which(adjmat==1),which(imat==1))

##
CheckEdge <- function(listgraph, edge)
  { 	l <- list()
  	for (i in 1:length(listgraph))
      	l[[i]] <- listgraph[[i]][edge]==1
      	return(any(unlist(l)==TRUE))
  }
##

	  for (j in 1:length(NrRemainingEdges))
	   { cat("j=",j,"\n")
		currentgraph[NrRemainingEdges[j]] <- edgesofg[NrRemainingEdges[j]]
     	if (CheckEdge(s0,NrRemainingEdges[j]))
      	  	{ mindist = mindist - 1
           	s0c <- s0
	   	s1c <- s1
	   	s2c <- s2
	   	s0 <- distdecrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$optimal
	   	s1 <- distdecrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal1
	   	s2 <- distdecrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal2
	   	s30 <- lapply(s0, function(x) ThreeNeighborhood(x))
	   	if (class(s30)=="list")
	   	s30 <- unique(unlist(s30, recursive=FALSE))
	   	s40 <- lapply(s0, function(x) FourNeighborhood(x))
	   	if (class(s40)=="list")
	   	s40 <- unique(unlist(s40, recursive=FALSE))
	   	s31 <- lapply(s1, function(x) TwoNeighborhood(x))
	   	if (class(s31)=="list")
	   	s31 <- unique(unlist(s31, recursive=FALSE))
	   	s41 <- lapply(s1, function(x) ThreeNeighborhood(x))
           	if (class(s41)=="list")
	   	s41 <- unique(unlist(s41, recursive=FALSE))
	   	s32 <- lapply(s2, function(x) OneNeighborhood(x))
	   	if (class(s32)=="list")
	   	s32 <- unique(unlist(s32, recursive=FALSE))
	   	s42 <- lapply(s2, function(x) TwoNeighborhood(x))
	   	if (class(s42)=="list")
		s42 <- unique(unlist(s42, recursive=FALSE))
	   	s3 <- unique(c(s30,s31,s32))
           	hd <- sapply(s3, function(x) hamming.distance(currentgraph,as.vector(x)))
	   	hd <- which(hd==mindist + 3)
	   	s3 <- s3[hd]
	   	s43 <- lapply(s3, function(x) OneNeighborhood(x))
	   	if (class(s43)=="list")
	   	s43 <- unique(unlist(s43, recursive=FALSE))
	   	s4 <- unique(c(s40,s41,s42,s43))
	   	hd <- sapply(s4, function(x) hamming.distance(currentgraph,as.vector(x)))
           	hd <- which(hd==mindist + 4)
	   	s4 <- s4[hd]
	  	}
 	 else 
         if (CheckEdge(s1,NrRemainingEdges[j]))
            	{ mindist = mindist 
	    	s0c <- s0
 	    	s1c <- s1
	    	s2c <- s2
	    	s0 <- distsame(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$optimal
            	s1 <- distsame(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal1
	    	s2 <- distsame(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal2
	    	s30 <- lapply(s0, function(x) ThreeNeighborhood(x))
	   	if (class(s30)=="list")
	   	s30 <- unique(unlist(s30, recursive=FALSE))
	   	s40 <- lapply(s0, function(x) FourNeighborhood(x))
	   	if (class(s40)=="list")
	   	s40 <- unique(unlist(s40, recursive=FALSE))
	   	s31 <- lapply(s1, function(x) TwoNeighborhood(x))
	   	if (class(s31)=="list")
	   	s31 <- unique(unlist(s31, recursive=FALSE))
	   	s41 <- lapply(s1, function(x) ThreeNeighborhood(x))
           	if (class(s41)=="list")
	   	s41 <- unique(unlist(s41, recursive=FALSE))
	   	s32 <- lapply(s2, function(x) OneNeighborhood(x))
	   	if (class(s32)=="list")
	   	s32 <- unique(unlist(s32, recursive=FALSE))
	   	s42 <- lapply(s2, function(x) TwoNeighborhood(x))
	   	if (class(s42)=="list")
	   	s42 <- unique(unlist(s42, recursive=FALSE))
	   	s3 <- unique(c(s30,s31,s32))
           	hd <- sapply(s3, function(x) hamming.distance(currentgraph,as.vector(x)))
	   	hd <- which(hd==mindist + 3)
	   	s3 <- s3[hd]
	   	s43 <- lapply(s3, function(x) OneNeighborhood(x))
	   	if (class(s43)=="list")
	   	s43 <- unique(unlist(s43, recursive=FALSE))
	   	s4 <- unique(c(s40,s41,s42,s43))
	   	hd <- sapply(s4, function(x) hamming.distance(currentgraph,as.vector(x)))
           	hd <- which(hd==mindist + 4)
	   	s4 <- s4[hd]
	   	}
	 else 
	 if (CheckEdge(s2,NrRemainingEdges[j]))
           	{ mindist = mindist + 1
	   	s0c <- s0
 	    	s1c <- s1
	    	s2c <- s2
	    	s0 <- distincrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$optimal
	    	s1 <- distincrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal1
	    	s2 <- distincrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal2
	    	s30 <- lapply(s0, function(x) ThreeNeighborhood(x))
	   	if (class(s30)=="list")
	   	s30 <- unique(unlist(s30, recursive=FALSE))
	   	s40 <- lapply(s0, function(x) FourNeighborhood(x))
	   	if (class(s40)=="list")
	   	s40 <- unique(unlist(s40, recursive=FALSE))
	   	s31 <- lapply(s1, function(x) TwoNeighborhood(x))
	   	if (class(s31)=="list")
	   	s31 <- unique(unlist(s31, recursive=FALSE))
	   	s41 <- lapply(s1, function(x) ThreeNeighborhood(x))
           	if (class(s41)=="list")
	   	s41 <- unique(unlist(s41, recursive=FALSE))
	   	s32 <- lapply(s2, function(x) OneNeighborhood(x))
	   	if (class(s32)=="list")
	   	s32 <- unique(unlist(s32, recursive=FALSE))
	   	s42 <- lapply(s2, function(x) TwoNeighborhood(x))
	   	if (class(s42)=="list")
	   	s42 <- unique(unlist(s42, recursive=FALSE))
	   	s3 <- unique(c(s30,s31,s32))
           	hd <- sapply(s3, function(x) hamming.distance(currentgraph,as.vector(x)))
	   	hd <- which(hd==mindist + 3)
	   	s3 <- s3[hd]
	   	s43 <- lapply(s3, function(x) OneNeighborhood(x))
	   	if (class(s43)=="list")
	   	s43 <- unique(unlist(s43, recursive=FALSE))
	   	s4 <- unique(c(s40,s41,s42,s43))
	   	hd <- sapply(s4, function(x) hamming.distance(currentgraph,as.vector(x)))
           	hd <- which(hd==mindist + 4)
	   	s4 <- s4[hd]
	   	}
	  else
	 	{ mindist = mindist + 1
	    	s0c <- s0
 	    	s1c <- s1
	    	s2c <- s2
	    	s0 <- distincrease1(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$optimal
	    	s1 <- distincrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal1
	    	s2 <- distincrease(s0c,s1c,s2c,s3,s4, NrRemainingEdges[j])$suboptimal2
	    	s30 <- lapply(s0, function(x) ThreeNeighborhood(x))
	   	if (class(s30)=="list")
	   	s30 <- unique(unlist(s30, recursive=FALSE))
	   	s40 <- lapply(s0, function(x) FourNeighborhood(x))
	   	if (class(s40)=="list")
	   	s40 <- unique(unlist(s40, recursive=FALSE))
	   	s31 <- lapply(s1, function(x) TwoNeighborhood(x))
	   	if (class(s31)=="list")
	   	s31 <- unique(unlist(s31, recursive=FALSE))
	   	s41 <- lapply(s1, function(x) ThreeNeighborhood(x))
           	if (class(s41)=="list")
	   	s41 <- unique(unlist(s41, recursive=FALSE))
	   	s32 <- lapply(s2, function(x) OneNeighborhood(x))
	   	if (class(s32)=="list")
	   	s32 <- unique(unlist(s32, recursive=FALSE))
	   	s42 <- lapply(s2, function(x) TwoNeighborhood(x))
	   	if (class(s42)=="list")
	   	s42 <- unique(unlist(s42, recursive=FALSE))
	   	s3 <- unique(c(s30,s31,s32))
           	hd <- sapply(s3, function(x) hamming.distance(currentgraph,as.vector(x)))
	   	hd <- which(hd==mindist + 3)
	   	s3 <- s3[hd]
	   	s43 <- lapply(s3, function(x) OneNeighborhood(x))
	   	if (class(s43)=="list")
	   	s43 <- unique(unlist(s43, recursive=FALSE))
	   	s4 <- unique(c(s40,s41,s42,s43))
	   	hd <- sapply(s4, function(x) hamming.distance(currentgraph,as.vector(x)))
           	hd <- which(hd==mindist + 4)
	   	s4 <- s4[hd]
	  	}
	    
	 }
   return(list(Optimum=s0, MinDist=mindist))
}

## end of transitive.projection function
######################################################################################
## Other functions used in transitive.projection
remTwoEdges <- function(mat)
{
   # remove 2 edge
      n <- length(which(mat==1))
      rnbhd <- list()
      if (n > 2)
      {
      rep <- t(combn(which(mat==1),2))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 0 # remove 2 edge
       if (is.transitive(matswap, nrow(mat)))
      {
      rnbhd[[length(rnbhd) + 1]] <- matswap
      }
      }
      }
return(rnbhd)
}
##--------------------------------------------------------------------------------
transSubGr <- function(adjmat)
{       
        n <- ncol(adjmat)
	## transGr <- TwoNeighborhood(adjmat)
        transGr <- remTwoEdges(adjmat)
	edg <- lapply(transGr,function(x) which(x==1))
        x1 <- which(adjmat==1)
        if (any(sapply(edg,function(x) length(intersect(x,x1))==length(x)))==TRUE) return(transGr[[sample(which(sapply(edg, function(x) length(intersect(x,x1))==length(x))==TRUE),1)]])
        else 
        return(matrix(0,ncol=n,nrow=n))
}
##--------------------------------------------------------------------------------
EdgeEk <- function(adjmat, ek)
{
g <- sapply(adjmat, function(x) x[ek]==1)
resTRUE <- adjmat[which(g==TRUE)]
resFALSE <- adjmat[which(g==FALSE)]
return(list(resT=resTRUE,resF=resFALSE))
}
##--------------------------------------------------------------------------------
# distance decreases then s0 s1 s2
distdecrease <- function(S0, S1, S2, S3, S4, ek)
{
 s0 <- EdgeEk(S0,ek)$resT
 s1 <- EdgeEk(S1,ek)$resT
 s2 <- c(EdgeEk(S0,ek)$resF, EdgeEk(S2,ek)$resT)
 return (list(optimal=s0,suboptimal1= s1,suboptimal2= s2))
}
##--------------------------------------------------------------------------------
# distance remains the same 
distsame <- function(S0, S1, S2, S3, S4, ek)
{
 s0 <- EdgeEk(S1,ek)$resT
 s1 <- c(EdgeEk(S0,ek)$resF, EdgeEk(S2,ek)$resT)
 s2 <- c(EdgeEk(S1,ek)$resF, EdgeEk(S3,ek)$resT)
return (list(optimal=s0,suboptimal1= s1,suboptimal2= s2))
}
##--------------------------------------------------------------------------------
# distance increases by 1
distincrease <- function(S0, S1, S2, S3, S4, ek)
{
 s0 <- c(EdgeEk(S0,ek)$resF, EdgeEk(S2,ek)$resT)
 s1 <- c(EdgeEk(S1,ek)$resF, EdgeEk(S3,ek)$resT)
 s2 <- c(EdgeEk(S2,ek)$resF, EdgeEk(S4,ek)$resT)
 return (list(optimal=s0,suboptimal1= s1,suboptimal2= s2))
}
##--------------------------------------------------------------------------------
# distance increase by 1
distincrease1 <- function(S0, S1, S2, S3, S4, ek)
{
 s0 <- c(EdgeEk(S0,ek)$resF)
 s1 <- c(EdgeEk(S1,ek)$resF, EdgeEk(S3,ek)$resT)
 s2 <- c(EdgeEk(S2,ek)$resF, EdgeEk(S4,ek)$resT)
 return (list(optimal=s0,suboptimal1= s1,suboptimal2= s2))
}
##--------------------------------------------------------------------------------
# Compute the one neighborhood of the input graph
OneNeighborhood <- function(mat)
{
      # add 1 edges
      matd1 <- mat
      diag(matd1) <- 1
      # n <- length(which(matd1==0))
      anbhd <-list()
      if (length(which(matd1==0)) > 0)
      {
      rep <- t(combn(which(matd1==0),1))
      NrComb <- nrow(rep)
      #fmat <- sapply(1:NrComb, function(x) mat[x,] <-1)
      #ffmat <- sapply(fmat, function(x) is.transitive(x, nrow(fmat(x))))
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 1 # add 1 edge
      if (is.transitive(matswap, nrow(mat)))
      {
      anbhd[[length(anbhd) + 1]] <- matswap
      } 
      }
      }
      # remove 1 edge
      # n <- length(which(mat==1))
      rnbhd <- list()
      if (length(which(mat==1)) > 0)
      {
      rep <- t(combn(which(mat==1),1))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 0 # remove 1 edge
       if (is.transitive(matswap, nrow(mat)))
      {
      rnbhd[[length(rnbhd) + 1]] <- matswap
      }
      }
      }
 # return(list(anbhd=anbhd,rnbhd=rnbhd))
   return(c(anbhd,rnbhd))

}
##--------------------------------------------------------------------------------
# Compute two neighborhoods of the input graph
TwoNeighborhood <- function(mat)
{
 # add 2 edges
      matd1 <- mat
      diag(matd1) <- 1
      n <- length(which(matd1==0))
      anbhd <-list()
      if (n > 1)
      {
      rep <- t(combn(which(matd1==0),2))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <-mat
      matswap[rep[i,]] <- 1 # add 2 edge
      if (is.transitive(matswap, nrow(mat)))
      {
      anbhd[[length(anbhd) + 1]] <- matswap
      } 
      }
      }
      # remove 2 edge
      n <- length(which(mat==1))
      rnbhd <- list()
      if (n > 2)
      {
      rep <- t(combn(which(mat==1),2))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 0 # remove 2 edge
       if (is.transitive(matswap, nrow(mat)))
      {
      rnbhd[[length(rnbhd) + 1]] <- matswap
      }
      }
      }
# add 1 edge remove 1 edges
      matd1 <- mat
      diag(matd1) <- 1
      r <- length(which(mat==1))
      a <- length(which(matd1==0))
      ranbhd <-list()
      if (r > 0 & a > 0)
      {
      rep <- t(combn(which(mat==1),1))
      rep1 <- which(matd1==0)
      NrComb <- nrow(rep)
      NrComb1 <- length(rep1)
      for (i in 1:NrComb)
      { 
      matswap <- mat
      matswap[rep[i,]] <- 0 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j]] <- 1
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      ranbhd[[length(ranbhd) + 1]] <- matswapz
      }
      }
      }
      }
return(c(anbhd,rnbhd, ranbhd))
}
##--------------------------------------------------------------------------------
# Compute three neighborhoods of the input graph
ThreeNeighborhood <- function(mat)
  {
  SwapIndexOld <- function(vec,ind)
      {
        res=vec
        if (res[ind]==1) res[ind]=0 else
        res[ind]=1
        return(res)
      }
      # add 3 edges
      matd1 <- mat
      diag(matd1) <- 1
      n <- length(which(matd1==0))
      anbhd <-list()
      if (n > 2)
      {
      rep <- t(combn(which(matd1==0),3))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <-mat
      matswap[rep[i,]] <- 1 # add 3 edges
      if (is.transitive(matswap, nrow(mat)))
      {
      anbhd[[length(anbhd) + 1]] <- matswap
      } 
      }
      }
      # remove 3 edges
      n <- length(which(mat==1))
      rnbhd <- list()
      if (n > 2)
      {
      rep <- t(combn(which(mat==1),3))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 0 # remove 3 edges
       if (is.transitive(matswap, nrow(mat)))
      {
      rnbhd[[length(rnbhd) + 1]] <- matswap
      }
      }
      }
      # add 2 edges remove 1 edge
      matd1 <- mat
      diag(matd1) <- 1
      a <- length(which(matd1==0))
      r <- length(which(mat==1))
      arnbhd <-list()
      if (a > 1 & r > 0)
      {
      rep <- t(combn(which(matd1==0),2))
      rep1 <- which(mat==1)
      NrComb <- nrow(rep)
      NrComb1 <- length(rep1)
      for (i in 1:NrComb)
      { 
       matswap <- mat
      matswap[rep[i,]] <- 1 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j]] <- 0
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      arnbhd[[length(arnbhd) + 1]] <- matswapz
      }
      }
      }
      }
      # add 1 edge remove 2 edges
      matd1 <- mat
      diag(matd1) <- 1
      r <- length(which(mat==1))
      a <- length(which(matd1==0))
      ranbhd <-list()
      if (r > 1 & a > 0)
      {
      rep <- t(combn(which(mat==1),2))
      rep1 <- which(matd1==0)
      NrComb <- nrow(rep)
      NrComb1 <- length(rep1)
      for (i in 1:NrComb)
      { 
      matswap <- mat
      matswap[rep[i,]] <- 0 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j]] <- 1
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      ranbhd[[length(ranbhd) + 1]] <- matswapz
      }
      }
      }
      }

           return(c(anbhd,rnbhd, arnbhd, ranbhd))
	   #s <- list(anbhd,rnbhd, arnbhd, ranbhd)
	   #return(list(anbhd,rnbhd, arnbhd, ranbhd))
	   
  }
##--------------------------------------------------------------------------------
# Compute the four neighborhood of the input graph
FourNeighborhood <- function(mat)
  {
  SwapIndexOld <- function(vec,ind)
      {
        res=vec
        if (res[ind]==1) res[ind]=0 else
        res[ind]=1
        return(res)
      }
      # add 4 edges
      matd1 <- mat
      diag(matd1) <- 1
      n <- length(which(matd1==0))
      anbhd <-list()
      if (n > 3)
      {
      rep <- t(combn(which(matd1==0),4))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <-mat
      matswap[rep[i,]] <- 1 # add 4 edges
      if (is.transitive(matswap, nrow(mat)))
      {
      anbhd[[length(anbhd) + 1]] <- matswap
      } 
      }
      }
      # remove 4 edges
      n <- length(which(mat==1))
      rnbhd <- list()
      if (n > 3)
      {
      rep <- t(combn(which(mat==1),4))
      NrComb <- nrow(rep)
      for (i in 1:NrComb)
      {
      matswap <- mat
      matswap[rep[i,]] <- 0 # remove 4 edges
       if (is.transitive(matswap, nrow(mat)))
      {
      rnbhd[[length(rnbhd) + 1]] <- matswap
      }
      }
      }
      # add 3 edges remove 1 edge
      matd1 <- mat
      diag(matd1) <- 1
      a <- length(which(matd1==0))
      r <- length(which(mat==1))
      arnbhd <-list()
      if (a > 2 & r > 0)
      {
      rep <- t(combn(which(matd1==0),3))
      rep1 <- which(mat==1)
      NrComb <- nrow(rep)
      NrComb1 <- length(rep1)
      for (i in 1:NrComb)
      { 
       matswap <- mat
      matswap[rep[i,]] <- 1 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j]] <- 0
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      arnbhd[[length(arnbhd) + 1]] <- matswapz
      }
      }
      }
      }
      # add 1 edge remove 3 edges
      matd1 <- mat
      diag(matd1) <- 1
      r <- length(which(mat==1))
      a <- length(which(matd1==0))
      ranbhd <-list()
      if (r > 2 & a > 0)
      {
      rep <- t(combn(which(mat==1),3))
      rep1 <- which(matd1==0)
      NrComb <- nrow(rep)
      NrComb1 <- length(rep1)
      for (i in 1:NrComb)
      { 
      matswap <- mat
      matswap[rep[i,]] <- 0 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j]] <- 1
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      ranbhd[[length(ranbhd) + 1]] <- matswapz
      }
      }
      }
      }
 # add 2 edges remove 2 edge
      matd1 <- mat
      diag(matd1) <- 1
      a <- length(which(matd1==0))
      r <- length(which(mat==1))
      aarrnbhd <-list()
      if (a > 1 & r > 1)
      {
      rep <- t(combn(which(matd1==0),2))
      rep1 <- t(combn(which(mat==1),2))
      NrComb <- nrow(rep)
      NrComb1 <- nrow(rep1)
      for (i in 1:NrComb)
      { 
       matswap <- mat
      matswap[rep[i,]] <- 1 
      for (j in 1:NrComb1)
      {
      matswapz <- matswap
      matswapz[rep1[j,]] <- 0
      
       if (is.transitive(matswapz, nrow(mat)))
      {
      aarrnbhd[[length(aarrnbhd) + 1]] <- matswapz
      }
      }
      }
      } 
      return(c(anbhd,rnbhd, arnbhd, ranbhd, aarrnbhd))
  }
##--------------------------------------------------------------------------------
## to check for transitivity
VecToMat <- function(vec,n)
{
  mat=matrix(vec,nrow=n)
  return(mat)
}
  
is.transitive <- function(graph,n, verb=FALSE)
  {
    # change to matrix representation
    if(is.vector(graph)) mat=VecToMat(graph,n)
    if(is.matrix(graph)) mat=graph
    if(!is.matrix(graph) & !is.vector(graph)) stop("please input a  vector or an adjacency matrix \n")
    
    
    # check transitivity from this node
    node.trans<-function(row)
      {
        if(sum(row)==0) return(TRUE)
        tmat<- mat[as.logical(row),,drop=FALSE]
        for(i in 1:dim(tmat)[1])
        {
          if( any(row < tmat[i,]) ) return(FALSE)
        }
        return(TRUE)
      }
    
     # do the check for all nodes
     for(i in 1:n)
       {
          if( !node.trans(mat[i,]) ) return(FALSE)
          if(verb) print(paste("checking node:", i))
       }
     return(TRUE)
  }
##--------------------------------------------------------------------------------









