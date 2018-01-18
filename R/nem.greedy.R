get.insertions = function(Phi, trans.close=TRUE){
    idx = which(Phi == 0)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible new edges
            Phinew = Phi
            Phinew[idx[i]] = 1
            if(trans.close)
                Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE) 
            models[[i]] <- Phinew
        }
    } 
    models       
}

get.deletions = function(Phi){
    Phi = Phi - diag(ncol(Phi))
    idx = which(Phi == 1)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible edge deletions
            Phinew = Phi
            Phinew[idx[i]] = 0
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    } 
    models       
}

get.reversions = function(Phi){
    idx = which(Phi + t(Phi) == 1, arr.ind=TRUE)
    models = list()
    if(NROW(idx) > 0){
        for(i in 1:NROW(idx)){ # test all possible edge reversions
            Phinew = Phi
            Phinew[idx[i,1],idx[i,2]] = 0
            Phinew[idx[i,2],idx[i,1]] = 1
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    } 
    models       
}

nem.greedy <- function(D, initial=NULL,control, verbose=TRUE){
    if(control$pcombi==TRUE){
        map<-control$map
        Sgenes<-colnames(map)
    }else{
      if(is(D, "list"))
        Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D[[1]]))]),"time")
      else
        Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
    }
        
	
    n <- length(Sgenes)     
    cat("Greedy hillclimber for",n,"S-genes (lambda =", control$lambda,")...\n\n")
    if(is.null(initial))
        Phi <- matrix(0,nrow=n,ncol=n)      
    else
        Phi = initial               
    if(control$pcombi==TRUE){
        diag(Phi) <- 0
    }else{
        diag(Phi) <- 1
    }
    dimnames(Phi) <- list(Sgenes,Sgenes)            
    sco0 <-
    nem(D, models=list(Phi), inference="search", control=control,verbose=verbose)$mLL
    finished <- FALSE
    while(!finished){
        models <- list()
#       propose new edges
        models = get.insertions(Phi, control$trans.close)
        if(control$type %in% c("CONTmLLMAP", "CONTmLLRatio") & !control$trans.close){ # these graphs are NOT transitively closed necessarily
               models = c(models, get.deletions(Phi), get.reversions(Phi))               
        }
        if(control$pcombi==TRUE){
            models <- models[unlist(lapply(models,function(x){is.DAG(x)}))]
        }
        models <- unique(models)
        if(verbose)
            cat(length(models), " local models to test ...\n")
        if(length(models) > 0){
            sconew <-
            nem(D, models=models, inference="search", control=control, verbose=verbose)          
            if(max(sconew$mLL) > sco0){
                if(verbose)
                    cat("--> Edge added, removed or reversed\n")
                sco0 <- max(sconew$mLL)
                Phi <- as(sconew$graph,"matrix")            
            }
            else # otherwise no improving edge could be inserted
                finished <- TRUE
        }else
            finished <- TRUE    
    }
    if(control$backward.elimination & !control$trans.close){               
        if(verbose)
            cat("Backward elimination step:\n\n")
        finished <- FALSE
        while(!finished){
    #       delete edges        
            idx = which(Phi - diag(n) == 1)
            if(length(idx) > 0){
                models <- list()
                for(i in 1:length(idx)){ # test all possible deletions
                    Phinew = Phi
                    Phinew[idx[i]] = 0
                    models[[i]] <- Phinew
                }
                if(control$pcombi==TRUE){
                    models <- models[unlist(lapply(models,function(x){is.DAG(x)}))]
                }
                models <- unique(models)
                sconew <- nem(D, models=models, inference="search", control=control, verbose=verbose)   
                if(max(sconew$mLL) > sco0){
                    if(verbose)
                        cat("--> Edge deleted\n")
                    sco0 <- max(sconew$mLL)
                    Phi <- as(sconew$graph,"matrix")            
                }
                else # otherwise no improving edge could be deleted
                    finished <- TRUE
            }else
                finished <- TRUE    
        }
    }
    ep <- nem(D, models=list(Phi), inference="search", control=control,verbose=FALSE)
        res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],control=control,selected=ep$selected, LLperGene=ep$LLperGene[[1]], para=ep$para[[1]])  # output: data likelihood under given model!    
    class(res) <- "nem.greedy"
    if(verbose)
        cat("log-likelihood of model = ",res$mLL,"\n")
    return(res)
}
