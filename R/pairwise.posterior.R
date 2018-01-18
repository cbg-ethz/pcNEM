##########################################
# Estimate the hierarchy edge by edge
# in each step only a pair of nodes is involved

pairwise.posterior = function (D, control, verbose = TRUE)
{
  # Sgenes
  	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
    nrS <- length(Sgenes)
    nrTest <- nrS*(nrS-1)/2
    if(verbose) cat(nrS,"perturbed genes ->", nrTest, "pairwise tests (lambda = ", control$lambda,")\n")
    
    # priors
    if (is.null(control$Pm)) control$Pm <- matrix(0, ncol=nrS, nrow=nrS, dimnames=list(Sgenes,Sgenes))		
    if (is.null(control$Pmlocal)) control$Pmlocal <- rep(0.25,4)
    if(is.null(control$Pe)){ 
	      control$Pe <- matrix(1/nrS,nrow=nrow(D),ncol=nrS, dimnames=list(rownames(D),Sgenes))  	
    }
    
     # init output
    graph <- diag(nrS)
    dimnames(graph) <- list(Sgenes,Sgenes)
    scores <- matrix(nrow=nrTest,ncol=5)
    dimnames(scores) <- list(as.character(1:nrTest),c("..","->","<-","<->","support")) 
    ix <- 1
        
    for (i in 1:(nrS - 1)) {
        for (j in (i + 1):nrS) {
            # get data
            x <- Sgenes[i]
            y <- Sgenes[j]
            sel <- which(colnames(D)==x | colnames(D)==y)            
            D.xy <- D[, sel]
            sel2 <- which(rowSums(D.xy) != 0)
            D.xy <- D.xy[sel2, , drop = FALSE]            
            
            # get local priors
            sel <- c(i,j)
	          controltmp = control
            controltmp$Pe <- control$Pe[sel2, sel, drop=FALSE] 
            controltmp$Pe[rowSums(controltmp$Pe) == 0,] = 1e-10 
            controltmp$Pm <- control$Pm[sel, sel, drop=FALSE]                        
            support <- nrow(D.xy)            
            
            # four models per edge: x..y  x->y  x<-y  x<->y
      	    models <- enumerate.models(2,name=c(x,y),trans.close=control$trans.close, verbose=FALSE)
            if(support > 0){            	
            	# score            	
          		ss <- score(models, D.xy, controltmp, verbose = FALSE)
          		post <- exp(ss$mLL) *control$Pmlocal
          		post <- post/sum(post)		
          		post[is.na(post)] = 0           
	           }
      	    else
      	    	post = matrix(0,nrow=1,ncol=4)
      	    	
	     # winner
            winner <- models[[which.max(post)]]
            graph[i, j] <- winner[1, 2]
            graph[j, i] <- winner[2, 1]
            
            # scores
            scores[ix, ] <- c(post, support)
            rownames(scores)[ix] <- paste(c(x, y), collapse = "~")
            ix <- ix + 1
            
            # counter
            if (verbose) cat(".")
        }
    }
    if (verbose) cat("\n")        
    
     # estimate effect positions
    if (verbose) cat("estimating effect positions\n")
    if(control$trans.close)
	    graph = transitive.closure(graph,mat=TRUE)
    ep <- score(list(graph), D,     	       
               control,
               verbose=FALSE)  
     
    # output
#     graph <- graph - diag(nrS)          
#     graph <- as(graph,"graphNEL")
#     res <- list(graph=graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],scores=scores,type=type,para=para,hyperpara=hyperpara,lam=lambda,selected=ep$selected,delta=delta)
	res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control, selected=ep$selected, LLperGene=ep$LLperGene[[1]], scores=scores, para=ep$para[[1]])	# output: data likelihood under given model!	
    class(res) <- "pairwise"
    return(res)
}

