PhiDistr <- function(Phi, Pm, a=1, b=0.1){		# essentially an exponential distribution over lambda, which is integrated out
	d <- abs(Phi - Pm)
	pPhi = a/(2*b)*(1 + d/b)^(-a-1)	
	sum(log(pPhi))
}