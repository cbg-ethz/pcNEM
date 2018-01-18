# log marginal likelihood of model
mLL <- function(Phi,D1,D0=NULL,control, verbose=FALSE) {
  if(control$pcombi==FALSE){
    if (!all(diag(Phi)==1))
       diag(Phi) <- 1
    Phi <- Phi[colnames(D1),colnames(D1)]
    if(control$selEGenes.method == "regularization"){
       Phi2 = cbind(Phi, double(NROW(Phi)))
       colnames(Phi2)[ncol(Phi2)] = "null"
    }
    else
       Phi2 = Phi
  }else{
    KOmap <- control$map
    Sgenes <- colnames(KOmap)
    cmap <- 1-KOmap
    Phi2 <- getperturb.matrices(cmap,Phi,control)$PropMat
  }
  
  para=NULL
  if(control$type == "mLL"){
    if(control$pcombi==TRUE){
       L <- matrix(NA,nrow = nrow(D1),ncol=ncol(Phi2))
       colnames(L) <- colnames(Phi2)
       if(!is.null(rownames(D1))){rownames(L) <- rownames(D1)}
       Pred1 <- (control$para[1]*(1-Phi2) + (1-control$para[2])*Phi2)
       Pred0 <- (control$para[2]*Phi2 + (1-control$para[1])*(1-Phi2))
       for(l in 1:nrow(D1)){
          L[l,] <- apply(((D1[l,]*Pred1) + (D0[l,]* Pred0)),2,prod)
       }
    }else
        L  <- control$para[1]^(D1 %*% (1-Phi2)) * (1-control$para[1])^(D0 %*% (1-Phi2)) * (1-control$para[2])^(D1 %*% Phi2) * control$para[2]^(D0 %*% Phi2)
  }
  else if(control$type %in% c("CONTmLLBayes", "CONTmLLDens")){		 
	  L <- exp(D1%*%Phi2)	
	#Q=apply(Phi,2, function(x) rowSums((matrix(rep(x, nrow(D1)), byrow=TRUE, ncol=4)*control$Pe[,1:NCOL(D1)])*D1))			
  }  
  else if(control$type == "CONTmLL")
	  L <- exp(log(D1)%*%Phi2 + log((1-D1))%*%(1-Phi2)) 
  else if(control$type %in% c("CONTmLLMAP", "CONTmLLRatio")){			
  	ep = D1%*%Phi2 + log(control$Pe)	
  	Theta = apply(ep,1,function(e) e ==max(e))
  	L = t((Phi2%*%(Theta*1)>0)*1)*D1
  	LLperGene=rowSums(L)		
  	s = sum(LLperGene)			
  	map = apply(Theta,1,which)	
  }
  else if(control$type == "depn"){	
  	res = score.network(D1, Phi, control)
  	s = res$loglik
  	LLperGene = res$LLperSample
  	ep = effect.likelihood(D1, res$net)		
  	Theta = apply(ep,1,function(e) e ==max(e))
  	map = apply(Theta,1,which)
  	para = res$net$parameters
  }
  if(!(control$type %in% c("CONTmLLMAP","CONTmLLRatio", "depn"))){	    
  	if(!is.null(control$Pe))
  		LP <- L*control$Pe[,1:NCOL(L)]
  	else
  		LP <- L	
  	LLperGene = log(rowSums(LP))
  	#LLperGene = rowSums(Q)
  	ep <- LP/(rowSums(LP))  	  			
  	Theta = apply(ep,1,function(e) e ==max(e))
  	s  <- sum(LLperGene)	
    	map = apply(Theta,1,which)			
    }   
    if(!is.null(rownames(D1))){
      map = sapply(map, names)}
    if(!is.null(control$Pm)){
    	if(control$lambda != 0){
  		if(verbose) cat("--> Using regularization to incorporate prior knowledge\n")			
    		s <- s - control$lambda*sum(abs(Phi - control$Pm)) + ncol(Phi)^2*log(control$lambda*0.5)    
  	}
  	else{
  		if(verbose) cat("--> Using Bayesian model averaging to incorporate prior knowledge\n")
  		s = s + PhiDistr(Phi, control$Pm)
  	}
  }
  list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene, para=para)
}

