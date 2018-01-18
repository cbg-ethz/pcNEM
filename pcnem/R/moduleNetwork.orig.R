moduleNetwork.orig <- function(D,control,verbose=TRUE){              
	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")  
    n <- length(Sgenes) 
    P <- sapply(Sgenes, function(x) rowSums(D[, colnames(D) == x, drop=FALSE])) 
    cat("Estimating module network of",n,"S-genes (lambda =", control$lambda,")...\n\n")    
    C = cor(P)          
    if(is.null(control$Pe)){
        control$Pe <- matrix(1/n,nrow=nrow(D),ncol=n)   
        colnames(control$Pe) <- Sgenes                  
    }       
    modeltotal <- matrix(0,ncol=n,nrow=n)
    colnames(modeltotal) <- Sgenes
    rownames(modeltotal) <- Sgenes  
    
    Dis0 = as.dist(1-C)             
    hc <- hclust(Dis0,method="average")         
    clusterData <- function(h){
        cl <- cutree(hc,h=h)
        T = table(cl)
        list(clust = cl, sizes = T, height=h)
    }   
    heights = sort(hc$height,decreasing=TRUE)   
    clust <- lapply(heights, clusterData)       
    maxsizes = sapply(clust, function(x) max(x$sizes))
    minsizes = sapply(clust, function(x) min(x$sizes))
    hmin = min(which(maxsizes < 5))
    hmax = min(which(minsizes < 5)) 
    idx = list()    
    for(h in hmax:hmin){        
        if(verbose){                
            plot(hc)    
            if(h > 1)
                rect.hclust(hc,h=clust[[h-1]]$height)                           
        }                                   
        small = which(clust[[h]]$sizes < 5)     
        for(c in small){
            variables = which(clust[[h]]$clust == c)            
            variables = setdiff(variables,unlist(idx)) # remove those, which have already been part of a module
            if(length(variables) > 0){
                if(verbose){
                    if(h > 1)
                        rect.hclust(hc,h=clust[[h-1]]$height,which=c,border="green")    
                }
                idx = c(idx, list(variables))       
                vars <- colnames(P)[variables]                          
                varidx <- sapply(vars,function(x) which(colnames(D) %in% x))    
                controltmp = control
                controltmp$Pe = control$Pe[,variables]
                controltmp$Pm = control$Pm[variables,variables]
                modeltotal <- moduleNetwork.aux(D[,varidx,drop=FALSE],modeltotal,variables,controltmp,verbose=verbose)
            }   
        }       
    }
    if(length(idx) > 1)
        res <- connectModules.orig(D, modeltotal, idx, control,verbose=verbose)
    else{
        if(control$trans.close)
            modeltotal <- transitive.closure(modeltotal, mat=TRUE,loops=TRUE)    
        diag(modeltotal) <- 0   
        ep <- score(list(modeltotal),D,control,verbose=FALSE)       
    # output        
        res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])   # output: data likelihood under given model!    
        class(res) <- "ModuleNetwork"
        if(verbose){
            cat("log-likelihood of model = ",res$mLL,"\n")
        }
    }
    return(res)
}

connectModules.aux = function(D, Phi, mod1, mod2, control, verbose){
	models = list()
	k = 1
	variables = colnames(Phi)[c(mod1, mod2)]
	controltmp = control
	controltmp$Pe = control$Pe[,variables]
	controltmp$Pm = control$Pm[variables,variables]  
	varidx <- sapply(variables,function(x) which(colnames(D) %in% x)) 
	sco0 <- score(list(Phi[variables,variables]),D[,varidx,drop=FALSE],controltmp,verbose=verbose,graphClass="matrix")$mLL # onyl focus on the subnetwork between modules mod1 and mod2
	bi = bincombinations(length(mod1)*length(mod2))
	
	aux2 = function(a, b){
		fkt1 <- function(x, a, b) {
			M = Phi[variables,variables]
			M[a, b] <- x		
			if(control$trans.close)    
				M <- transitive.closure(M,mat=TRUE,loops=TRUE)    
			return(list(M))
		}
    a = match(a, c(mod1, mod2))
    b = match(b, c(mod1, mod2))
		models <- apply(bi,1,fkt1, a, b) 
		models <- unique(matrix(unlist(models),ncol=length(variables)^2,byrow=TRUE))
		
		fkt2 <- function(x){
			M <- matrix(x,length(variables))
			dimnames(M) <- list(variables,variables)
			return(list(M))
		}
		models <- unlist(apply(models,1,fkt2),recursive=FALSE)
		
		sconew <- score(models,D[,varidx,drop=FALSE], controltmp, verbose=verbose,graphClass="matrix")         
		maidx = which.max(sconew$mLL)
		if(sconew$mLL[maidx] > sco0){		
			sco0 <- sconew$mLL[maidx]
			modellocal <- sconew$graph			
		}	
		else
			modellocal = matrix(0, ncol=NCOL(sconew$graph), nrow=NROW(sconew$graph))
		modellocal
	}
	
	modellocal1 = aux2(mod1, mod2)
	modellocal2 = aux2(mod2, mod1)
  allmod = c(mod1, mod2)
	Phi[mod1, mod2] = modellocal1[match(mod1, allmod), match(mod2, allmod)]
	Phi[mod2, mod1] = modellocal2[match(mod2, allmod), match(mod1, allmod)]
	Phi
}

connectModules.orig <- function(D, Phi, modules, control, verbose=TRUE){
	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
    n <- length(Sgenes) 
    if(verbose){    
        dev.off()
        cat("Connecting module pairs ...\n\n")  
    }              
    for(i in 1:length(modules)){ # nur Kanten in eine Richtung!
        mod1 = modules[[i]]  
		if(i < length(modules)){
        	for(j in (i+1):length(modules)){
          		mod2 = modules[[j]]             
           		Phi = connectModules.aux(D, Phi, mod1, mod2, control, verbose)
			}
        }        
    }       
    ep <- score(list(Phi),D, control, verbose=FALSE,graphClass="graphNEL")
    res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])  # output: data likelihood under given model!    
    class(res) <- "ModuleNetwork"
    if(verbose)
        cat("log (posterior) (marginal) likelihood of model = ",res$mLL,"\n")
    return(res)
}
