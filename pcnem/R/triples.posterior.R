triples.posterior <- function(D, control,verbose=TRUE){

  # Sgenes
  Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
  nrS    <- length(Sgenes)
  nrTest <- choose(nrS,3) 
  cat(nrS,"perturbed genes ->", nrTest, "triples to check (lambda = ",control$lambda,")\n")

  # local model prior
  if (is.null(control$Pmlocal)) control$Pmlocal <- rep(1/29,29)
 
 
  ## 
  ## 1. learn tripel models
  ##
 
  ##=== create all (ordered) triples
  triples = subsets(nrS,3)

  ##=== 29 candidate models for each triple
  cnd.models   = enumerate.models(3,trans.close=control$trans.close,verbose=FALSE)
  
  fkt2 <- function(x,name){	  
	  dimnames(x) <- list(name,name)
	  x
  }
  ##=== store maximum model in a list
  mll.models   = list()
  if (verbose) cat(".") 
  for(i in 1:nrow(triples)){

	sel <- which(colnames(D)%in%Sgenes[triples[i,]])
    D.tmp <- D[,sel]	
	cnd.models <- lapply(cnd.models,fkt2, Sgenes[triples[i,]])
#        D.tmp <- D.tmp[rowSums(D.tmp)!=0,]	
#        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,1]]] = "a"
#        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,2]]] = "b"
#        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,3]]] = "c"	
	controltmp = control
    controltmp$Pe <- control$Pe[,sel,drop=FALSE]            
	if(!is.null(control$Pm))    
        controltmp$Pm <- control$Pm[sel,sel,drop=FALSE]				
   	tmp.mdl = score(cnd.models,D.tmp, controltmp,verbose=FALSE)	
 
        ##== do prior stuff
    post = tmp.mdl$mLL + log(control$Pmlocal)
    winner <- cnd.models[[which.max(post)]]
    mll.models[[i]] = list()
    mll.models[[i]]$graph = winner
    mll.models[[i]]$posterior = post
    dimnames(mll.models[[i]]$graph) = list(Sgenes[triples[i,]],Sgenes[triples[i,]]) #=== rename node to inds
    if (verbose) cat(".")
  }

  ##
  ## 2. break down to pairwise-data by model averaging
  ##

  A = matrix(NA, ncol=nrS, nrow=nrS)
  dimnames(A) = list(Sgenes,Sgenes)

  diag(A) = 1
  for(i in 1:(nrS-1)){  
  for(j in (i+1):nrS){
        ##=== find contributing triples
        inds = apply(triples,1,function(x) sum( x %in% c(i,j)) == 2)
        inds = which(inds)
        ##=== count edges and set A_ij accordingly
        tmp = list()
        for(k in 1:length(inds)){ ##=== for each contributing triple
                tmp[[k]] = mll.models[[inds[k]]]$graph #=== get adjacency matrix
        }
        ##=== mean number of edges from i to j (including doubles...)
        A[i,j] = mean(unlist(lapply(tmp,function(x) x[Sgenes[i],Sgenes[j]])))
        ##=== mean number of edges from j to i (including doubles...) 
        A[j,i] = mean(unlist(lapply(tmp,function(x) x[Sgenes[j],Sgenes[i]])))
  }}
  B = (A>=control$triples.thrsh)*1-diag(nrow(A)) 
  if(verbose) cat("\n")
  ##
  ## 3. estimate effect positions
  ##
  if (verbose) cat("Estimating effect positions in combined graph\n")    
    if(control$trans.close)
	B = transitive.closure(B,mat=TRUE)
  diag(B) = 0
  graph <- as(B,"graphNEL")
  ep <- score(list(B), D,     	       
               control,
               verbose=FALSE)  
  ##
  ## 4. output
  ##
  res <- list(graph=graph,avg=A,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])
  class(res) <- "triples"
  if(verbose)
	cat("log-likelihood of model = ",res$mLL,"\n")
  return(res)
}
