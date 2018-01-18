nem.greedyMAP <- function(D, control, verbose=TRUE){		
	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
	n <- length(Sgenes)		
	cat("Alternating optimization for",n,"S-genes (lambda =", control$lambda,")...\n\n")
	Theta = apply(D,1, function(e) (e == max(e)) & (e>0))*1	  
	best = NULL
	converged = FALSE
	control$type="CONTmLLMAP"
	i = 1		
	while(!converged){
# 		if(trans.close & (n < 9)){
# 			cvec = as.vector(t(Theta%*%D))									
# 			T = diag(n)			
# 			combis = combn(which(T==0),3)
# 			bvec = rep(1,ncol(combis)+n)
# 			Amat = matrix(0,nrow=ncol(combis),ncol=length(cvec))
# 			for(c in 1:ncol(combis)){
# 				Amat[c,combis[,c]] = c(1,1,-1)
# 			}					
# 			ones = which(T==1)			
# 			for(c in 1:length(ones)){				
# 				tmp = double(n^2)
# 				tmp[ones[c]] = 1
# 				Amat = rbind(Amat, tmp)
# 			} 
# 			
# 			Phi0 = lp("max", cvec, Amat, c(rep("<=",ncol(combis)),rep("=",n)), bvec, all.bin=TRUE)		
# 			Phi = matrix(Phi0$solution, byrow=FALSE, ncol=n)			
# 			dimnames(Phi) = list(Sgenes, Sgenes)	
# 			print(Phi)
# 			readline()					
# 		}
# 		else
			Phi = apply(Theta%*%D, 1, function(e) (e > 0))*1		
		sco <- nem(D,models=list(Phi), inference="search", control=control,verbose=verbose)
		if(verbose){
			cat("iteration ", i, ": likelihood = ", sco$mLL, "\n")
		}
		if(i > 1)		
			converged = all(Phi==Phiold)
		if(converged)
			break;	
		Theta = apply(sco$pos[[1]],1,function(e) (e == max(e)))*1				
		Theta = Theta[1:(nrow(Theta)-1),]				
		Phiold = Phi		
		i = i + 1
	}			
	Phi = as(sco$graph, "matrix")	
	if(control$trans.close)	
		Phi = closest.transitive.greedy(Phi, verbose)
	ep <- nem(D,models=list(Phi), inference="search", control=control,verbose=FALSE)
    	res <-
	list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])	# output: data likelihood under given model!	)	# output: data likelihood under given model!	
	class(res) <- "nem.greedyMAP"	
	return(res)
}
