# density function
dbum <- function(x,a,lambda){		
	lambda[1] +  lambda[2]*dbeta(x,a[1],1) + lambda[3]*dbeta(x,1,a[2])
}

# gradient of BUM w.r.t. parameters: WRONG -> we need the deriv of the LOG-likelihood!
#bum.gradient = function(x, a, lambda){
#	lambda[2]*(sum(log(x)) - length(x)*(-digamma(a[1] + a[2]) + digamma(a[1]))) + lambda[3]*(sum(log(1-x)) - length(x)*(-digamma(a[1] + a[2]) + digamma(a[2])))
#}

# cdf
pbum <- function(x,a,lambda){					
	lambda[1]*x + lambda[2]*pbeta(x,a[1],1) + lambda[3]*pbeta(x,1,a[2])
}

getComponent <- function(x,i,a){	
	if(i == 1)
		return(double(length(x))+1)	
	else if(i == 2)
		return(dbeta(x,a[1],1))
	else if(i == 3)
		return(dbeta(x,1,a[2]))
}

# quantile function
qbum <- function(p,a,lambda,nbisect=20){			
	n <- length(p)
	mid <- rep(0,n)
	top <- rep(1,n)
	bot <- rep(0,n)
	gohigher <- rep(FALSE,n)
	for (j in 1:nbisect){
		mid <- (top+bot)/2
		gohigher <- (pbum(mid,a,lambda)<p)		
		bot[gohigher] <- mid[gohigher]
		top[!gohigher] <- mid[!gohigher]
	}
	return(mid)
}

rbum <- function(n,a,lambda){
	u<-runif(n)
	return(qbum(u,a,lambda))
}

logit <- function(x){
	return(log(x)-log(1-x))
}

inv.logit <- function(x){
	return(exp(x)/(1+exp(x)))
}

# negative log-likelihood for given data set X with given hyperparameters
bum.negLogLik <- function(hyper,X,lambda){		
	a = inv.logit(hyper[1]) + 1e-10
	b = exp(hyper[2])+2 + 1e-10		
	nLL = -sum(log(dbum(X,c(a,b),lambda))) - dexp(b, 0.1, log=TRUE) - dbeta(a,1,2, log=TRUE)
	#cat("a = ", a, "b = ", b, "\n")	  	
	nLL
}

# MLE estimates of parameters 
bum.mle <- function(X,a,lambda){			
	res <- optim(c(logit(a[1]),log(a[2]-2+1e-10)), bum.negLogLik,X=X,lambda=lambda, method="L-BFGS-B", upper=10, lower=-10)
	list(a=c(inv.logit(res$par[1]) + 1e-10,exp(res$par[2])+2 + 1e-10),lambda=lambda,logLik=-res$value+1e-50)
}

bum.EM <- function(X,starta=c(0.3,10),startlam=c(0.6,0.1,0.3), tol=1e-4){			
	converged <- FALSE		
	a <- starta	
	lambda <- startlam
	ncomp <- length(startlam)
	cat("start values: a = ", a, "lambda = ",lambda,"\n")
	logLik <- 0
	Z <- matrix(0,nrow=length(X),ncol=ncomp)	
	while(!converged){	
# 	E-step		
		for(i in 1:ncomp)
			Z[,i] <- lambda[i]*getComponent(X,i,a)/dbum(X,a,lambda)		
		lambdanew <- apply(Z,2,mean)						
# 	M-step		
		paras <- bum.mle(X,a,lambdanew)	
							
		converged <- (abs(logLik/paras$logLik - 1) < tol)
		a <- paras$a
		lambda <- lambdanew
		logLik <- paras$logLik						
		cat("logLik = ", logLik, "a = ",a, "lambda = ",lambda,"\n")
	}
	print("converged!")	
	return(list(a=a,lambda=lambda,Z=Z,logLik=logLik))	
}

qqbum<-function(pvals,a,lambda,main="BUM QQ Plot",xlab="BUM Expected p-value",ylab="Observed p-value"){	
	n <- length(pvals)
	pvals <- sort(pvals)		
	plot(c(0,1),c(0,1),main=main,xlab=xlab,ylab=ylab,type="n")
	lines(qbum((rank(pvals)-.5)/n,a,lambda),pvals,lty=2)		
	lines(c(0,1),c(0,1))					
}

# histogram plot with imposed distribution plotted atop
bum.histogram <- function(pvalues,a,lambda,main="Histogram",xlab="p-value",ylab="Density",breaks="Sturges"){	
	hist(pvalues,probability=TRUE,main=main,xlab=xlab,ylab=ylab,breaks)
	x <- 1:100/100
	lines(x,dbum(x,a,lambda),lwd=3)
	lines(x,bum.dalt(x,a,lambda),lwd=2,col="red")
}


# density function of alternative distribution
bum.dalt <- function(x,a,lambda,normalize=TRUE){		
	pihat = dbum(1,a,lambda) # uniform component	
	dalt = dbum(x,a,lambda) - pihat
	if(normalize)
		dalt = dalt/(1-pihat)
	dalt
}

# cdf of alternative distribution
bum.palt <- function(x,a,lambda){
	pihat = dbum(1,a,lambda) # uniform component
	(pbum(x,a,lambda) - pihat*x)/(1-pihat)
}

# quantile function
bum.qalt <- function(p,a,lambda,nbisect=20){				
	n <- length(p)
	mid <- rep(0,n)
	top <- rep(1,n)
	bot <- rep(0,n)
	gohigher <- rep(FALSE,n)
	for (j in 1:nbisect){
		mid <- (top+bot)/2
		gohigher <- (bum.palt(mid,a,lambda)<p)		
		bot[gohigher] <- mid[gohigher]
		top[!gohigher] <- mid[!gohigher]
	}
	return(mid)
}


# draw n random samples from alternative distribution
bum.ralt <- function(n,a,lambda){
	u <- runif(n)
	return(bum.qalt(u,a,lambda))
}