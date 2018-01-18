score.aux <- function(models, D, control, verbose=TRUE, graphClass="graphNEL") {

  #if single model as input
  if (class(models)=="matrix") models <- list(models)    

  # Which Sgenes were silenced?
  if(control$pcombi== FALSE){
    Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
    nrS <- length(Sgenes)
  }else{
    map<-control$map
    nexp<-rownames(map)
    Sgenes<-colnames(map)
    nrS <- length(Sgenes)
  }
  
  # check that all models have S-genes as names
  fkt <- function(x,s){
     ss <- sort(s)
     c1 <- all(sort(setdiff(colnames(x), "unknown"))==ss)
     c2 <- all(sort(setdiff(rownames(x), "unknown"))==ss)
     return(c1 & c2)
  }  
  if (!all(sapply(models,fkt,s=Sgenes))) stop("\nnem:score> models must have same names as data")
  if (control$pcombi== TRUE){
    if(length(nexp) != ncol(D)){
      stop("\nnem:score> the number of experiments in KO map do not match with the data")
      }
    }
  #nrS <- length(Sgenes)
  # make probability/density matrices D0 and D1  
  # nrow=#E-genes and ncol=#S-genes        
  if(control$type %in% c("mLL", "FULLmLL")){
      if(control$pcombi==TRUE){
        D1 = sapply(nexp, function(s) rowSums(D[,colnames(D) == s,drop=FALSE]))  #Sgenes
        D0 = sapply(nexp, function(s) sum(colnames(D) == s)) - D1
      }else{
        D1 = sapply(Sgenes, function(s) rowSums(D[,colnames(D) == s,drop=FALSE]))
        D0 = sapply(Sgenes, function(s) sum(colnames(D) == s)) - D1
      }
  }else{
	D1 = D
	D0 = NULL
  }
    # if no prior is supplied:
  # assume uniform prior over E-gene positions      
  if (is.null(control$Pe)){ 	
	  control$Pe <- matrix(1/nrS,nrow=NROW(D1),ncol=nrS)
	  colnames(control$Pe) <- Sgenes  		
  }  
  if("unknown" %in% colnames(models[[1]])){
    control$Pe = cbind(control$Pe, 1/nrS)
    control$Pe = control$Pe/rowSums(control$Pe)  	
    colnames(control$Pe)[ncol(control$Pe)] = "unknown"
  }    
  if(control$selEGenes.method == "regularization" && ncol(control$Pe) == nrS){
  	control$Pe = cbind(control$Pe, double(nrow(D1)))
  	control$Pe[,ncol(control$Pe)] = control$delta/nrS
    control$Pe = control$Pe/rowSums(control$Pe)		
    colnames(control$Pe)[ncol(control$Pe)] = "null"
  }    
  if(is.null(control$Pm) & (control$lambda != 0)){
	cat(">>> Regularization parameter non-zero: Generating sparsity prior automatically! <<<\n")
	control$Pm = diag(length(Sgenes))
  }
     
  if (control$type=="FULLmLL"){ # FULL log marginal likelihood of all models
    if (verbose==TRUE) cat("Computing FULL (marginal) likelihood for",length(models),"models\n")
   	if(control$lambda != 0)
    	results <- sapply(models,FULLmLL,D1,D0,control, verbose)
	else{
		if ("doMC" %in% loadedNamespaces()){
			registerDoMC(control$mc.cores)
			results = foreach(m = models) %dopar%
				FULLmLL(m, D1,D0,control, verbose)
		}
		else{
			results <- sapply(models,FULLmLL,D1,D0,control, verbose) 
		}
	}
  }
  else{   # log marginal likelihood of all models		
	if (verbose==TRUE) cat("Computing (marginal) likelihood for",length(models),"models\n")
	if(control$lambda != 0)
		results <- sapply(models,mLL,D1,D0,control, verbose)     	
	else{
		if ("doMC" %in% loadedNamespaces())
			results = foreach(m = models) %dopar%
				mLL(m, D1,D0,control, verbose)
		else
			results <- sapply(models,mLL,D1,D0,control, verbose)
	}
  }
  if(control$lambda != 0 | !("doMC" %in% loadedNamespaces())){	  
	  s       <- unlist(results["mLL",])
	  ep      <- results["pos",]
	  map     <- results["mappos",] 	  
	  LLperGene = results["LLperGene",]
	  para = results["para",]
	  selected = results["mappos",which.max(s)][[1]]
  }
  else{
	  s = sapply(results, function(r) r$mLL)
	  ep   <- lapply(results, function(r) r$pos)
	  map = lapply(results, function(r) r$mappos)	
	  LLperGene = lapply(results, function(r) r$LLperGene)
	  para = lapply(results, function(r) r$para)	
	  selected = map[[which.max(s)]] 
  }  
  selected = unique(unlist(selected[Sgenes]))
#   if(!is.null(Pm)){  	
#   	log_pD_cond_Phi <- s  	  	  		
#   	if(is.null(control$lambda) || (control$lambda == 0)){  			
# 		if(verbose) cat("--> Using Bayesian model averaging to incorporate prior knowledge\n")
#   		lpPhi <- sapply(models, PhiDistr, Pm, a=1, b=0.5)		  		
#   		s <- log_pD_cond_Phi + lpPhi 		
#   	}
#   	else
# 		s = s - control$lambda*sapply(models, function(M) sum(abs(M - control$Pm))) + nrS^2*log(control$lambda*0.5)  
#   }     
  if(verbose){
	if(length(s) > 1){		
		mLL.sorted = sort(s, decreasing=TRUE)
		cat("((Marginal) posterior likelihood difference of best vs. second best model for ", Sgenes, ":", mLL.sorted[1] - mLL.sorted[2],")\n")
	}
  }  
   # winning model       
  winner <- models[[which.max(s)]]  
  diag(winner) <- 0  
  if(graphClass == "graphNEL"){
  	gR <- new("graphAM",adjMat=winner[Sgenes,Sgenes],edgemode="directed")  
  	gR <- as(gR,"graphNEL")    
  }
  else
	gR <- winner  
  res <- list(graph=gR, mLL=s, pos=ep, mappos=map, control=control, selected=selected, LLperGene=LLperGene, para=para)
  class(res) <- "score"   
  return(res)  
}

score = function(models, D, control, verbose=TRUE, graphClass="graphNEL") {	
	if(is(D, "matrix")){
		return(score.aux(models, D, control, verbose, graphClass))
	}
	cat("Cauton: NEM inference with several datasets is experimental so far!\n")
	if(is(D, "list")){
		if(!is.null(control$Pe) & !(is(control$Pe, "list") & length(control$Pe) == length(D)))
			stop("There has to be one E-gene prior for each data set")
		mLLs = 0
		pos = list()
		mappos = list()
		selected = list()
		LLperGene = list()
		scs = list()
		para = list()
		for(i in 1:length(D)){
			control.tmp = control
			if(i > 1){ # S-gene prior is added just once
				control.tmp$lambda = 0
				control.tmp$Pm = NULL
			}			
			control.tmp$Pe = control$Pe[[i]]
			scs[[i]] = score.aux(models, D[[i]], control.tmp, verbose, graphClass)
			mLLs = mLLs + scs[[i]]$mLL
			pos[[i]] = scs[[i]]$pos
			mappos[[i]] = scs[[i]]$mappos
			LLperGene[[i]] = scs[[i]]$LLperGene
			selected[[i]] = scs[[i]]$selected
			para[[i]] = scs[[i]]$para
		}		
		winner <- models[[which.max(mLLs)]]  
		diag(winner) <- 0  
		if(graphClass == "graphNEL"){
			gR <- new("graphAM",adjMat=winner,edgemode="directed")  
			gR <- as(gR,"graphNEL")    
		}
		else
			gR <- winner  		
		sc.consensus = list(graph=gR, mLL=mLLs, pos=pos, mappos=mappos, control=control, selected=selected, LLperGene=LLperGene, para=para)
		class(sc.consensus) = "score.list"
		return(sc.consensus)
	}
	else
		stop("data has to be either a list of matrices or one matrix")
}
