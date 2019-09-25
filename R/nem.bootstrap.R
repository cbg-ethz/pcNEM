nem.bootstrap <- function(D, thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE){
	if(inference == "dynoNEM")
		stop("nem.bootstrap is not applicable for dynoNEMs")
	if(is(D, "list")){
		if(length(D) > 1)
			return(nem.bootstrap.list(D, thresh, nboot, inference, models, control, verbose))
		else
			D = D[[1]]
	}
	#inferNetwork <- function(idx.orig=1:nrow(D), boot){				
	#	controltmp = control				
	#	controltmp$Pe = control$Pe[boot,]		
	#	Dtmp = D[boot,]
	#	if(!is.null(control$Pm) & length(control$lambda) > 1)
	#		res = as.vector(as(nemModelSelection(control$lambda,Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
	#	else
	#		res = as.vector(as(nem(Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
	#	res
	#}
	#results = bootstrap(1:nrow(D),nboot,theta=inferNetwork)$thetastar	
  inferNetwork <- function(idx.orig=1:nrow(D), boot){				
    controltmp = control				
    controltmp$Pe = control$Pe[boot,]		
    Dtmp = D[boot,]
    if(control$pcombi==TRUE){
      restmp = nem(Dtmp,inference,models,controltmp,verbose)
      res <- as.vector(as(restmp$graph,"matrix"))
      resnoise <- c(restmp$typeIEst, restmp$typeIIEst)
      c(res, resnoise)
    }else{
      if(!is.null(control$Pm) & length(control$lambda) > 1)
        res = as.vector(as(nemModelSelection(control$lambda,Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
      else
        res = as.vector(as(nem(Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
      res
    }
  }
	if("multicore" %in% loadedNamespaces())
		use.parallel = "multicore"
	else if("snow" %in% loadedNamespaces()){
		use.parallel = "snow"
		cl <- makeCluster(control$mc.cores, type = "SOCK")
	}
	else
		use.parallel = "no"
	cat("--> parallel bootstrap: ", use.parallel, "\n")	
	if(!is.null(rownames(D)) & "time"  %in% colnames(D)){		
		group = as.factor(paste(rownames(D), D[,"time"],sep=""))		
		res.boot = boot::boot(1:nrow(D), inferNetwork, nboot, strata=group, parallel=use.parallel, ncpus=control$mc.cores, cl=cl)
	}
	else
		res.boot = boot::boot(1:nrow(D), inferNetwork, nboot, parallel=use.parallel, ncpus=control$mc.cores, cl=cl)
	results =  res.boot$t
	if(control$pcombi==TRUE){
	  map<-control$map
	  Sgenes<-colnames(map)
	}else{
	  Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
	}
	n = length(Sgenes)
	overlapBoot = colMeans(results)
	if(control$pcombi==TRUE){
	  avgNoise <- tail(overlapBoot,2)
	  overlapBoot = matrix(round(overlapBoot[1:n^2],digits=2),ncol=n,nrow=n)
	  control$para <- avgNoise
	}else{
	  overlapBoot = matrix(round(overlapBoot,digits=2),ncol=n,nrow=n)
	}
	#overlapBoot = matrix(round(overlapBoot,digits=2),ncol=n,nrow=n)
	colnames(overlapBoot) = Sgenes
	rownames(overlapBoot) = Sgenes
	print(overlapBoot)
	control$lambda = 0
	control$Pm = NULL
	control$Pe = NULL
	if(control$pcombi==TRUE){
	  if(is.DAG((overlapBoot>thresh)*1)){# Since the bootstrap graph need not necessarily be a DAG and pcnem can't handle cycles
	    res = nem(D,models=list((overlapBoot>thresh)*1),inference="search",control, verbose=verbose) 
	  }else{
	    print("The consensus graph is not a DAG. Searching for the nearest DAG")
	    startDAG <- as(nem(D,inference="AdaSimAnneal",control=control,verbose=FALSE)$graph,'matrix')
	    startDAG <- startDAG[Sgenes,Sgenes]
	    resNearestDAG <- SearchNearestDAG(n,startDAG,D,control,(overlapBoot>thresh)*1)
	    print(paste0("The distance between consensus and nearest DAG from the search is ",resNearestDAG$bestDistance))
	    res = nem(D,models=list(resNearestDAG$maxDAG),inference="search",control, verbose=verbose)
	    res$bestDistance <- resNearestDAG$bestDistance
	  }
	}else
	  res = nem(D,models=list((overlapBoot>thresh)*1),inference="search",control, verbose=verbose)
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	res$para = res$para[[1]]
	res$control= control
	res$overlapBoot = overlapBoot
	res$allRuns = results
	g = res$graph
	edgeDataDefaults(g, "label") = 1	
	edgeDataDefaults(g, "weight") = 1
	for(s1 in Sgenes){
		for(s2 in Sgenes){
			if(s2 %in% unlist(adj(g, s1))){
				edgeData(g, from = s1, to = s2, attr = "weight") = overlapBoot[s1,s2]			
				edgeData(g, from = s1, to = s2, attr = "label") = overlapBoot[s1,s2]
			}
		}
	}
	res$graph = g
	class(res) <- "nem.bootstrap"
	res
}

nem.bootstrap.list = function(D, thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE){
  if(control$pcombi==TRUE){
    map<-control$map
    Sgenes<-colnames(map)
  }else{
    Sgenes <- setdiff(unlist(control$map[intersect(names(control$map), colnames(D))]),"time")
  }
  
	n = length(Sgenes)
	res.all = matrix(0, ncol=nboot, nrow=n^2)
	for(b in 1:nboot){
		Dsam = list()
		Pe.sam = list()
		for(i in 1:length(D)){
			sam = sample(1:NROW(D[[i]]), NROW(D[[i]]), replace=TRUE)
			Dsam[[i]] = D[[i]][sam, ]			
			Pe.sam[[i]] = control$Pe[[i]][sam, ]			
		}			
		control.tmp = control
		if(!is.null(control$Pe))
			control.tmp$Pe = Pe.sam		
		if(!is.null(control$Pm) & length(control$lambda) > 1)
			res.all[,b] = as.vector(as(nemModelSelection(control$lambda,Dsam,inference,models,control.tmp,verbose)$graph,"matrix"))
		else
			res.all[,b] = as.vector(as(nem(Dsam,inference,models,control.tmp,verbose)$graph,"matrix"))		
	}
	overlapBoot = rowMeans(res.all)
	overlapBoot = matrix(round(overlapBoot,digits=2),ncol=n,nrow=n)
	colnames(overlapBoot) = Sgenes
	rownames(overlapBoot) = Sgenes
	print(overlapBoot)
	control$lambda = 0
	control$Pm = NULL
	res = nem(D,models=list((overlapBoot>thresh)*1),inference="search",control, verbose=verbose)
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	res$para = res$para[[1]]
	res$control= control
	res$overlapBoot = overlapBoot
	g = res$graph
	edgeDataDefaults(g, "label") = 1	
	edgeDataDefaults(g, "weight") = 1
	for(s1 in Sgenes){
		for(s2 in Sgenes){
			if(s2 %in% unlist(adj(g, s1))){
				edgeData(g, from = s1, to = s2, attr = "weight") = overlapBoot[s1,s2]			
				edgeData(g, from = s1, to = s2, attr = "label") = overlapBoot[s1,s2]
			}
		}
	}
	res$graph = g
	class(res) <- "nem.bootstrap"
	res
}
