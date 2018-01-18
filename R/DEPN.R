# Deterministic Effects Propagation Networks for RPPA data with missing values (Fröhlich et al., BMC Bioinformatics, 2009)

encode.interventions = function(datall, I, map){    
    antibodies = setdiff(colnames(datall), "time")
    interventions = paste("I", antibodies, sep="_")
    dat = cbind(datall, matrix(0, ncol=length(antibodies), nrow=nrow(datall)))
    colnames(dat) = c(colnames(datall), interventions)
    rownames(dat) = rownames(datall)
    exps = unique(rownames(dat))
    for(c in exps){     
        effected = unlist(sim.interventions(I, c, map))         
        effected = intersect(effected, antibodies)
        if(length(effected) > 0)
            dat[rownames(dat) %in% c,paste("I",effected,sep="_")] = 1   
    }   
    dat
}

sim.interventions = function(M, int, map){
    diag(M) = 1 
    int = unlist(map[int])
    int = intersect(int, rownames(M))   
    effected = lapply(int, function(i) c(i, colnames(M)[M[i,] == 1]))
}

# dngamma = function(mu, tau, nu, lambda, alpha, beta){
#   beta^alpha/(gamma(alpha)*sqrt(2*pi*lambda)) * tau^(alpha - 0.5) * exp(-beta*tau) * exp(-tau*(mu - nu)^2/(2*lambda))
# }

learn.conditionals = function(node, datall, control){       
    if("time" %in% colnames(datall))
        time = unique(as.character(datall[,"time"]))
    else
        time = "0"  
    cond = list(method=control$mode)    
    if(all(is.na(datall[, node])))
        return(cond)
    if(control$debug)       
        pdf(file=paste(control$outputdir,"/",node,sep=""))      
    for(t in time){
        if("time" %in% colnames(datall))
            denst = datall[datall[,"time"] == t, ]
        else
            denst = datall      
        cond[[t]] = list(intervention=c(), no_intervention=c())                     
        dens1 = denst[denst[,paste("I",node,sep="_")]==1, node]
        dens0 = denst[denst[,paste("I",node,sep="_")]==0, node]         
        if(control$mode == "continous_ML"){             
            if(length(dens1) > 0)
                cond[[t]]$intervention = c(mu=mean(dens1, na.rm=TRUE), sd=sd(dens1, na.rm=TRUE))
            if(length(dens0) > 0)
                cond[[t]]$no_intervention = c(mu=mean(dens0, na.rm=TRUE), sd=sd(dens0, na.rm=TRUE))
        }
        else{
            if(length(dens1) > 0){          
                # the sample parameters
                m = mean(dens1, na.rm=TRUE)
                s2 = var(dens1, na.rm=TRUE)
                n = length(dens1)                       
                #  parameters of the posterior normal scaled inverse gamma distribution for (mu, sigma^2): Gelman, S. 79
                mu_n = control$lambda.intervention/(control$lambda.intervention + n) * control$nu.intervention + n/(control$lambda.intervention + n)*m
                kappa_n = control$lambda.intervention + n
                nu_n = control$df.intervention + n
                nu_n_sigma2 = (control$df.intervention*control$scale.intervention + (n-1)*s2 + (control$lambda.intervention*n)/(control$lambda.intervention + n)*(m - control$nu.intervention)^2)               
                cond[[t]]$intervention=c(mu=mu_n, sd=sqrt((nu_n_sigma2)/(nu_n-2)), mu.sd=sqrt((nu_n_sigma2)/(kappa_n*(nu_n+2))), var.sd=sqrt(2*nu_n_sigma2^2/((nu_n-2)^2*(nu_n-4)))) # for mu|x and sigma^2|x plug in the expectations and the variances        
            }
            if(length(dens0) > 0){
                # the sample parameters
                m = mean(dens0, na.rm=TRUE)
                s2 = var(dens0, na.rm=TRUE)
                n = length(dens0)                           
                # parameters of the posterior normal scaled inverse gamma distribution for (mu, sigma^2): Gelman, S. 79
                mu_n = control$lambda.no_intervention/(control$lambda.no_intervention + n) * control$nu.no_intervention + n/(control$lambda.no_intervention + n)*m
                kappa_n = control$lambda.no_intervention + n
                nu_n = control$df.no_intervention + n
                nu_n_sigma2 = (control$df.no_intervention*control$scale.no_intervention + (n-1)*s2 + (control$lambda.no_intervention*n)/(control$lambda.no_intervention + n)*(m - control$nu.no_intervention)^2)
                cond[[t]]$no_intervention=c(mu=mu_n, sd=sqrt((nu_n_sigma2)/(nu_n-2)), mu.sd=sqrt((nu_n_sigma2)/(kappa_n*(nu_n+2))), var.sd=sqrt(2*nu_n_sigma2^2/((nu_n-2)^2*(nu_n-4)))) # for mu|x and sigma^2|x plug in the expectations and the variances     
            }           
        }
        if(control$debug){                              
            if(length(dens0) > 0){
                hist(dens0, probability=TRUE, main=paste(node, paste("(time=",t,", no intervention)",sep="")))
                x = seq(min(dens0, na.rm=TRUE), max(dens0, na.rm=TRUE), by=0.01)
                lines(x, dnorm(x, mean=cond[[t]]$no_intervention["mu"], sd=cond[[t]]$no_intervention["sd"]), lwd=2,col="red")
            }
            
    #       x11()       
            if(length(dens1) > 0){
                hist(dens1, probability=TRUE, main=paste(node, paste("(time=",t,", intervention)",sep="")))
                x = seq(min(dens1, na.rm=TRUE), max(dens1, na.rm=TRUE), by=0.01)
                lines(x, dnorm(x, mean=cond[[t]]$intervention["mu"], sd=cond[[t]]$intervention["sd"]), lwd=2,col="red")
            }   
    #       readline()
    #       graphics.off()          
        }       
    }       
    if(control$debug)       
        dev.off()
    cond        
}

sample.effect.likelihood = function(d, net){    
    nnodes = length(net$measure.nodes)
    if("time" %in% names(d))        
        t = as.character(d["time"])
    else
        t = "0"
    l = double(nnodes)
    names(l) = net$measure.nodes        
    for(i in net$measure.nodes){        
        cond = net$parameters[[i]]
        if(!is.null(cond[[t]]$intervention))
            l[i] =  dnorm(d[i], mean=cond[[t]]$intervention["mu"], sd=cond[[t]]$intervention["sd"]) 
    }       
    l
}

effect.likelihood = function(datall, net){
    nnodes = length(net$measure.nodes)
    L = matrix(0, ncol=nnodes, nrow=nrow(datall))
    for(i in 1:nrow(datall))
        L[i,] = sample.effect.likelihood(datall[i,], net)
    dimnames(L) = list(rownames(datall), net$measure.nodes)
    L
}

sample.likelihood = function(d, net){       
    nnodes = length(net$measure.nodes)
    if("time" %in% names(d))        
        t = as.character(d["time"])
    else
        t = "0"
    l = 0   
    for(i in net$measure.nodes){        
        cond = net$parameters[[i]]      
        if(!is.null(cond[[t]]$intervention) & !is.null(cond[[t]]$no_intervention)){
            int_i = paste("I",i,sep="_")    
            if(is.na(d[i])){            
                if(d[int_i] == 1)
                    d[i] = cond[[t]]$intervention["mu"]         
                else if(d[int_i] == 0)
                    d[i] = cond[[t]]$no_intervention["mu"]          
            }
            if(d[int_i] == 1)
                l = l + dnorm(d[i], mean=cond[[t]]$intervention["mu"], sd=cond[[t]]$intervention["sd"], log=TRUE)       
            else if(d[int_i] == 0)
                l = l + dnorm(d[i], mean=cond[[t]]$no_intervention["mu"], sd=cond[[t]]$no_intervention["sd"], log=TRUE)
        }       
    }       
    list(likelihood=l, x=d)
}

learn = function(datall, net, control){
    antibodies = net$measure.nodes
    net$parameters = lapply(antibodies, learn.conditionals, datall, control=control)    
    names(net$parameters) = antibodies          
    net 
}

data.likelihood = function(datall, net){    
    L = apply(datall, 1, sample.likelihood, net)
    LLperSample = sapply(L, function(l) l$likelihood)
    likelihood = sum(LLperSample)
    dall = t(sapply(L, function(l) l$x))
    list(likelihood=likelihood, data=dall, LLperSample=LLperSample)
}

score.network = function(datall, I, control){
    antibodies = setdiff(colnames(datall), "time")
    datall = encode.interventions(datall, I, control$map)       
    net = list(measure.nodes=antibodies, network=I)
    latent = apply(datall, 2, function(x) all(is.na(x)))
    missing = apply(datall, 2, function(x) any(is.na(x)))
    missing = setdiff(which(missing), which(latent))
#    cat("Scoring network with ", ncol(I), "nodes (", sum(latent), "latent nodes; ", length(missing), " nodes with missing values)\n")
    if(any(missing)){
        converged = FALSE
        loglik = -Inf
#        cat("Missing value imputation via EM algorithm:\n")
        dat.new = datall
		iter = 0
        while(!converged){
            net = learn(dat.new, net, control)  #M-step     
            res = data.likelihood(datall, net) # E-step            
			if(is.na(res$likelihood))
				res$likelihood = -Inf
			logliktmp = res$likelihood
            dat.new = res$data          
            if((logliktmp == loglik) || abs(logliktmp/loglik - 1) < 1e-4 || iter == 100)
                converged = TRUE
            loglik = logliktmp
 #           cat("log likelihood = ", loglik, "\n")
			iter = iter + 1
        }
  #      cat("converged!\n")             
    }
    else{
        net = learn(datall, net, control)   
        res = data.likelihood(datall, net)
        dat.new = datall
    }       
    ret = list(loglik=res$likelihood, LLperSample=res$LLperSample, net=net)
}

# learn.network = function(datall, init, map, method="ML", shape=0.1, lambda=0.1, nu.intervention=0.7, nu.no_intervention=0.95, rate.intervention=1, rate.no_intervention=0.1,...){
#   net = init  
#   score = score.network(datall, net, map, method=method, shape=shape, lambda=lambda, nu.intervention=nu.intervention, nu.no_intervention=nu.no_intervention, rate.intervention=rate.intervention, rate.no_intervention=rate.no_intervention, ...)
#   cat("initial score = ", score,"\n")
# #     add edges
#   converged = FALSE
#   while(!converged){
#       toset = which(net == 0)
#       scoretmp = double(length(toset))
#       if(length(toset) > 0){
#           for(i in 1:length(toset)){
#               nettmp = net
#               nettmp[toset[i]] = 1
#               scoretmp[i] = score.network(datall, nettmp, map, method=method, shape=shape, lambda=lambda,nu.intervention=nu.intervention, nu.no_intervention=nu.no_intervention, rate.intervention=rate.intervention, rate.no_intervention=rate.no_intervention, ...)
#           }
#           best = which.max(scoretmp)
#           if(scoretmp[best] > score){
#               net[toset[best]] = 1
#               score = scoretmp[best]
#               cat("edge added: score = ", score,"\n")
#           }
#           else
#               converged = TRUE
#       }
#       else
#           converged = TRUE
#   }
# #     delete edges
#   converged = FALSE
#   while(!converged){
#       toset = which(net - diag(ncol(net)) == 1)
#       scoretmp = double(length(toset))
#       if(length(toset) > 0){
#           for(i in 1:length(toset)){
#               nettmp = net
#               nettmp[toset[i]] = 0
#               scoretmp[i] = score.network(datall, nettmp, map, method=method, shape=shape, lambda=lambda, nu.intervention=nu.intervention, nu.no_intervention=nu.no_intervention, rate.intervention=rate.intervention, rate.no_intervention=rate.no_intervention, ...)
#           }
#           best = which.max(scoretmp)
#           if(scoretmp[best] > score){
#               net[toset[best]] = 0
#               score = scoretmp[best]
#               cat("edge deleted: score = ", score,"\n")
#           }
#           else
#               converged = TRUE
#       }
#       else
#           converged = TRUE
#   }
# #     net = as(transitive.closure(as(net,"graphNEL")), "matrix")  
#   net
# }

# bootstrap.network = function(datall, init, nboot=1000, ...){
#   B = matrix(0, ncol=ncol(init), nrow=nrow(init))
#   dimnames(B) = dimnames(init)    
#   for(b in 1:nboot){
#       sub = sample(1:nrow(datall), nrow(datall), replace=TRUE)
#       B = B + learn.network(datall[sub,], init, ...)
#   }
#   B = B/nboot
#   B
# }
