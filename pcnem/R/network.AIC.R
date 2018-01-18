# 	get AIC for network
network.AIC = function(network, Pm=NULL, k=length(nodes(network$graph)), verbose=TRUE){	
	M <- as(network$graph,"matrix")
	diag(M) = 0
	if(is.null(Pm))
		Pm = matrix(0,ncol=ncol(M),nrow=nrow(M))
	diag(Pm) = 0
	npar <- sum(abs(M - Pm)>0)
	AIC <- -2*network$mLL + k*npar
	if(verbose) 
		cat(paste("==> AIC ( lambda = ",network$lam,") = ",AIC,"( #param =", npar,")===============\n"))
	return(AIC)
}
