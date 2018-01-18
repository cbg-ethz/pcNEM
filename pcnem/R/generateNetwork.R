sampleRndNetwork = function(Sgenes, scaleFree=TRUE, gamma=2.5, maxOutDegree=length(Sgenes), maxInDegree=length(Sgenes), trans.close=TRUE, DAG=FALSE){		
	n = length(Sgenes)
	S = diag(n) # network of S-genes	
	maxOutDegree = min(n,maxOutDegree)	
	degprob = (0:maxOutDegree)^(-gamma)
	degprob[1] = 1
	degprob = degprob/sum(degprob)	
	for(i in 1:n){ # connect network randomly
		if(scaleFree)
			outdeg = sample(0:maxOutDegree,1,prob=degprob)# power law for out-degree => scale-free network
		else
			outdeg = sample(0:maxOutDegree,1)		
		if(outdeg > 0){
			if(!DAG)
				idx0 = which(S[i,] == 0)
			else
				idx0 = which(S[i,] == 0 & 1:n < i) # sample on lower triangle matrix
			if(length(idx0) > 0){
				idx = which(colSums(S[,idx0,drop=FALSE]) <= maxInDegree)
				if(length(idx) > 0){			
					idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
					S[i,idx] = 1					
				}			
			}	
		}
	}						
	if(trans.close)
		S = transitive.closure(S, mat=TRUE,loops=FALSE)					
	diag(S) = 0			
	colnames(S) = Sgenes
	rownames(S) = Sgenes
	S
}

sampleData = function(Phi, m, prob=NULL, uninformative=0, type="binary", replicates=4, typeI.err=0.05, typeII.err=0.2, alpha=sample(seq(0.1,0.9,by=0.1),ncol(Phi),replace=TRUE), beta=sample(5:50,ncol(Phi),replace=TRUE), lambda=matrix(sample(seq(0.01,0.49,by=0.01),ncol(Phi)*2,replace=TRUE),ncol=2), meansH1=rep(0.5, ncol(Phi)), meansH0=rep(-0.5, ncol(Phi)), sdsH1=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE), sdsH0=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE)){
	Sgenes = colnames(Phi)
	n = length(Sgenes)
	epos = sample(1:n,m,replace=TRUE,prob=prob)
	Theta = matrix(0, nrow=length(Sgenes), ncol=m)
	for(i in 1:ncol(Theta))
		Theta[epos[i],i] = 1
	Phi = transitive.closure(Phi, mat=TRUE, loops=TRUE)			
	M = t((Phi%*%Theta > 0)*1)	
	if(type == "binary"){
		D = matrix(0, ncol=n*replicates, nrow=m)
		k2 = 1
		for(i in 1:n){			
			D[M[,i] == 1, k2:(k2+replicates-1)] = matrix(sample(c(0,1),replicates*sum(M[,i]),replace=TRUE,prob=c(typeII.err,1-typeII.err)),ncol=replicates)# effected genes => aus H1 ziehen
			D[M[,i] == 0, k2:(k2+replicates-1)] =  matrix(sample(c(0,1),replicates*sum(M[,i] == 0),replace=TRUE,prob=c(1-typeI.err,typeI.err)), ncol=replicates) # ... and not effected ones => aus H0 ziehen
			k2 = k2 + replicates
		}
		if(uninformative > 0)
			D = rbind(D, matrix(sample(c(0,1),n*replicates*uninformative,replace=TRUE,prob=c(1-typeI.err,typeI.err)), nrow=uninformative, ncol=n*replicates)) # ... and not effected ones => aus H0 ziehen)
	}
	else if(type %in% c("density")){
		lambda = cbind(1-rowSums(lambda), lambda)				
		palt = sapply(1:n, function(i) bum.ralt(m, c(alpha[i], beta[i]), lambda[i,]))		
		p0 = matrix(runif((m+uninformative)*n), ncol=n)
		P = M*palt + (1-M)*p0[1:m,]
		if(uninformative > 0)
			P = rbind(P, p0[(m+1):nrow(p0),])
		D = sapply(1:n, function(i) bum.dalt(P[,i], c(alpha[i], beta[i]), lambda[i,]))
		D = log(D)
	}
	else if(type %in% c("lodds")){
		palt = sapply(1:n, function(i) rnorm(m, mean=meansH1[i], sd=sdsH1[i]))
		p0 = sapply(1:n, function(i) rnorm(m+uninformative, mean=meansH0[i], sd=sdsH0[i]))		
		D = M*palt + (1-M)*p0[1:m,]
		if(uninformative > 0)
			D = rbind(D, p0[(m+1):nrow(p0),])		
	}	
	else
		stop(paste("unknown type", type, "\n"))	
	if(type == "binary")
		colnames(D) = rep(Sgenes, each=replicates)
	else
		colnames(D) = Sgenes
	list(D=D, epos=Sgenes[epos])
}

sampleData.pcnem = function(Phi, M, prob=NULL, map,uninformative=0,typeI.err=0.05, typeII.err=0.2){
    Sgenes <- colnames(Phi)
    N <- length(Sgenes)
    epos = sample(1:N, M,replace=TRUE,prob=prob)
    Theta = matrix(0, nrow=length(Sgenes), ncol=M)
    for(i in 1:ncol(Theta))
    Theta[epos[i],i] = 1
    cmap <- 1-map # complementary probabilities
    Phi2 <- matrix(NA,nrow=nrow(map),ncol=ncol(map))
    if(any(diag(Phi) == 1)){ diag(Phi) <- 0} # no loops
    PathMat <- matrix(NA,nrow=N,ncol=N) # Number of paths from a node(rows) to target node(columns)
    dimnames(PathMat) <- dimnames(Phi)
    for(i in 1:N){
        for(j in 1:N){
            PathMat[i,j] <- PathNumber(Sgenes[i],Sgenes[j],Phi)
        }
    }
    for(k in 1:nrow(map)){
        for(n in 1:ncol(map)){
            g <- which(PathMat[,n] != 0)
            Phi2[k,n] <- 1-prod((cmap[k,g]^PathMat[g,n]))
        }
    }
    
    colnames(Phi2) <- colnames(Phi)
    rownames(Phi2) <- rownames(map)
    
    F.star <- t(Phi2%*%Theta)
    
    D = matrix(0, ncol=nrow(Phi2), nrow=ncol(Theta))
    
    for(x in 1:nrow(D)){
        for(y in 1:ncol(D)){
            D[x,y] <- sample(c(0,1),1,
            prob=c((F.star[x,y]*typeII.err+ (1-F.star[x,y])*(1-typeI.err)),
            ((F.star[x,y]*(1-typeII.err))+ ((1-F.star[x,y])*typeI.err))))
        }
    }
    
    if(uninformative > 0){
        D = rbind(D, matrix(sample(c(0,1),nrow(Phi2)*uninformative,replace=TRUE,
        prob=c(1-typeI.err,typeI.err)), nrow=uninformative, ncol=nrow(Phi2)))
    }
    colnames(D) <- rownames(Phi2)
    list(D = D, epos =Sgenes[epos])
}


# sampleData.gnem = function(net, m, int, map, t){	
# 	nnodes = length(net$measure.nodes)
# 	t = as.character(t)
# 	D = matrix(0, ncol=nnodes, nrow=m)
# 	if(length(int) == 0){
# 		for(i in 1:nnodes){		
# 			cond = net$parameters[[i]]									
# 			D[,i] = rnorm(m, mean=cond[[t]]$no_intervention["mu"], sd=cond[[t]]$no_intervention["sd"])	
# 		}
# 		return(D)
# 	}
# 	else{		
# 		effected = simulate.interventions(net, int, map)		
# 		for(i in 1:nnodes){		
# 			cond = net$parameters[[i]]			
# 			if(net$measure.nodes[i] %in% effected)
# 				D[,i] = rnorm(m, mean=cond[[t]]$intervention["mu"], sd=cond[[t]]$intervention["sd"])
# 			else
# 				D[,i] = rnorm(m, mean=cond[[t]]$no_intervention["mu"], sd=cond[[t]]$no_intervention["sd"])			
# 		}	
# 	}	
# 	D
# }


# sampleData.BN = function(core, reporters=40, nr_intven=3, beta1=0.9, nullnode=FALSE){		
# 	vert = ncol(core)	
# 	original = core
# 	core = transitive.closure(core, mat=T, loops=T)
#     	diag(core) = diag(original)
# 	intven = rep(1:vert, each = nr_intven)	
# 	Delta = t(core)[,intven]
# 	ind = sample(1:nrow(core), reporters, replace = TRUE)
# 	Delta = Delta[ind,]
# 	colnames(Delta) = paste("I", intven,sep=".")
# 	if(nullnode){
# 		l = nrow(core)
# 		coregraph = matrix(0,nrow = l+1, ncol = l+1)
# 		coregraph[1:l,1:l] = core
# 		rownames(coregraph) = c(unique(colnames(Delta)),"null")
# 		colnames(coregraph) = c(unique(colnames(Delta)),"null")
# 		core = coregraph}
# 	Theta = matrix(0, ncol=reporters, nrow =ncol(core))
# 	for(i in 1:ncol(Theta))
# 		{Theta[ind[i],i] = 1}
# 	if(nullnode){rownames(Theta) = c(unique(colnames(Delta)),"null")}
# 	else{rownames(Theta) = unique(colnames(Delta))} 
# 	colnames(Theta) = as.character(1:ncol(Theta))
# 	report = rbinom(length(Delta), size = 1, prob = beta1)
# 	report[which(Delta == 0)] = 1 - report[which(Delta == 0)]
# 	dim(report) = dim(Delta)
# 	colnames(report)=colnames(Delta)
# 	rownames(report) = colnames(Theta)
# 	BN = createBN(data = report, coregraph = core, marginal = Theta, nullnode = nullnode)
# 
# 	D = BN$data	
# 	colnames(D) = rep(colnames(core), each=nr_intven)
# 	epos = apply(Theta>0,2,which)	
# 	list(D=D, epos=epos)
# }
