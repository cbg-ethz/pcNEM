moduleNetwork <- function(D,control,verbose=TRUE){              
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
                if(typeof(varidx)=='list'){
                  modeltotal <- moduleNetwork.aux(D[,unlist(varidx),drop=FALSE],modeltotal,variables,controltmp,verbose=verbose)
                }else
                  modeltotal <- moduleNetwork.aux(D[,varidx,drop=FALSE],modeltotal,variables,controltmp,verbose=verbose)
                
            }   
        }       
    }
    if(length(idx) > 1)
        res <- connectModules2(D, modeltotal, idx, control,verbose=verbose)
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

moduleNetwork.aux <- function(D,modeltotal, variables, control, verbose=TRUE){
    if(verbose){
        cat("estimating network of genes:\n")
        cat(variables,"\n")
    }
    n <- length(variables)
    if(n > 1){          
        models <- enumerate.models(n,name=unique(colnames(D)), trans.close=control$trans.close, verbose=verbose)
        sco <- score(models,D,control,verbose=verbose,graphClass="matrix")
        modellocal <- sco$graph
    }
    else
        modellocal <- 0
    modeltotal[variables,variables] <- modellocal
    modeltotal
}

connectModules2 <- function(D, Phi, modules, control, verbose=TRUE){
	Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
    n <- length(Sgenes) 
    if(verbose){    
        dev.off()
        cat("Connecting module pairs via constraint hill climbing ...\n\n")  
    }
    sco0 <- score(list(Phi),D,control,verbose=verbose,graphClass="matrix")$mLL          
    BetweenModules = matrix(0,n,n)
    dimnames(BetweenModules) = dimnames(Phi)
    for(i in 1:length(modules)){
        mod1 = modules[[i]] 
        if(i < length(modules)){
            for(j in (i+1):length(modules)){
                mod2 = modules[[j]]             
                BetweenModules[mod1,mod2] = 1
                BetweenModules[mod2,mod1] = 1
            }
        }
    }   
    BetweenModules0 = BetweenModules
    finished <- FALSE
    while(!finished){       
#       propose new edges between modules                                   
        idx = which(BetweenModules == 1)
        if(length(idx) > 0){    
            models <- list()
            for(i in 1:length(idx)){ # test all possible new edges
                Phinew = Phi
                Phinew[idx[i]] = 1
                if(control$trans.close)
                    Phinew = transitive.closure(Phinew, mat=TRUE, loops=TRUE)   
                models[[i]] <- Phinew               
            }               
            models <- unique(models)
            sconew <- score(models,D, control, verbose=verbose,graphClass="matrix")         
            maidx = which.max(sconew$mLL)
            if(sconew$mLL[maidx] > sco0){
                newedges = which(sconew$graph - Phi == 1)
                BetweenModules[newedges] = 0
                sco0 <- sconew$mLL[maidx]
                Phi <- sconew$graph
            }
            else # otherwise no improving edge could be inserted
                finished <- TRUE
        }
        else
            finished <- TRUE    
    }
    if(!control$backward.elimination & !control$trans.close){
        if(verbose)
            cat("Backward elimination step:\n\n")
        finished <- FALSE
        while(!finished){       
    #       delete new edges between modules                                
            idx = which(BetweenModules0*Phi == 1)
            if(length(idx) > 0){    
                models <- list()
                for(i in 1:length(idx)){ # test all possible edge deletions
                    Phinew = Phi
                    Phinew[idx[i]] = 0                  
                    models[[i]] <- Phinew               
                }               
                models <- unique(models)
                sconew <- score(models,D, control, verbose=verbose,graphClass="matrix")         
                maidx = which.max(sconew$mLL)
                if(sconew$mLL[maidx] > sco0){
                    newedges = which(sconew$graph - Phi == 1)
                    BetweenModules[newedges] = 0
                    sco0 <- sconew$mLL[maidx]
                    Phi <- sconew$graph
                }
                else # otherwise no improving edge could be deleted
                    finished <- TRUE
            }
            else
                finished <- TRUE    
        }
    }
    ep <- score(list(Phi),D, control, verbose=FALSE,graphClass="graphNEL")
    res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])  # output: data likelihood under given model!    
    class(res) <- "ModuleNetwork"
    if(verbose)
        cat("log (posterior) (marginal) likelihood of model = ",res$mLL,"\n")
    return(res)
}
