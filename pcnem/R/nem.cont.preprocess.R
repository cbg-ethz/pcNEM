# D - expression data
nem.cont.preprocess <- function (D, neg.control = NULL, pos.control = NULL, nfold = 2,  
    influencefactor=NULL, empPval=.05, verbose = TRUE) {

    # check input
    if (is.null(neg.control) & is.null(pos.control)) 
        stop("\nnem> provide at least one control")        
    if (class(neg.control) == "matrix") 
        if (nrow(neg.control) != nrow(D)) 
            stop("\nnem> control and data must have the same number of rows")
    if (class(pos.control) == "matrix") 
        if (nrow(pos.control) != nrow(D)) 
            stop("\nnem> control and data must have the same number of rows")
    if (class(neg.control) %in% c("integer", "numeric") & class(pos.control) %in% 
        c("integer", "numeric") & !all(c(neg.control, pos.control)) %in% 
        1:ncol(D)) 
        stop("\nnem>controls not in data 'D'")
    if (class(neg.control) %in% c("integer", "numeric") & class(pos.control) %in% 
        c("integer", "numeric") & any(neg.control %in% pos.control)) 
        stop("\nnem>negative and positive controls overlap")
    if (!is.null(neg.control) & !is.null(pos.control)) {
        setting <- "twocontrols"
    }
    else {
        setting <- "onecontrol"
    }
    
    
    
    
    # pos and neg control available
    if (setting == "twocontrols") {
        if (verbose) 
            cat("preprocessing with respect to POS and NEG controls\n")
        if (class(neg.control) == "matrix") 
            neg <- neg.control
        if (class(neg.control) %in% c("integer", "numeric")) 
            neg <- D[, neg.control]
        if (class(pos.control) == "matrix") 
            pos <- pos.control
        if (class(pos.control) %in% c("integer", "numeric")) 
            pos <- D[, pos.control]
        if (class(neg.control) %in% c("integer", "numeric")) {
            d <- neg.control
        }
        else {
            d <- NULL
        }
        if (class(pos.control) %in% c("integer", "numeric")) 
            d <- c(d, pos.control)                            # controll indices
        if (!is.null(d)) 
            dat <- D[, -d]
            
        # select effector genes, here only high regulated genes
        sel <- which(exp(rowMeans(pos) - rowMeans(neg)) > nfold)
        
        dat.sel <- dat[sel, ]
        pos.sel <- pos[sel, ]
        neg.sel <- neg[sel, ]
        
        nrEgenes <- nrow(dat.sel)
        nrArrays <- ncol(dat.sel)
        
        prob.influenced   <- matrix(nrow=nrEgenes, ncol=nrArrays)        
        prob.pos.control <- matrix(nrow=nrEgenes, ncol=length(pos.control))
        
        
        # estimate influencefactor
        if (is.null(influencefactor)) {
            prob.neg.control <- matrix(nrow=nrEgenes, ncol=length(neg.control))
            influencefactor = 1
            flag = TRUE
            while (flag) {
                # calculate negative control probabilities
                for (i in seq(1, nrEgenes)) {
                    mean.neg = mean(neg.sel[i,])
                    mean.pos = mean(pos.sel[i,])
                    std.neg = sd(neg.sel[i,])
                    std.pos = sd(pos.sel[i,])
                    if ((std.neg+std.pos) < mean.pos-mean.neg) {
                        x <- ((mean.pos-mean.neg)-(std.neg+std.pos))/2
                        std.neg <- std.neg+x
                        std.pos <- std.pos+x
                    }

                    # neg control
                    for (j in seq(1:length(neg.control))) {
                        prob.neg.control[i,j] <- (dnorm(neg.sel[i,j], mean.neg, std.neg) / 
                            (dnorm(neg.sel[i,j], mean.neg, std.neg) + 
                            dnorm(neg.sel[i,j], mean.pos, std.pos))) * influencefactor
                        if (prob.neg.control[i,j] > 0.9999) 
                            prob.neg.control[i,j] = 0.9999
                        if (prob.neg.control[i,j] < 0.0001)
                            prob.neg.control[i,j] = 0.0001
                    }
                }
                if (sum(prob.neg.control < 0.5)>0)
                    influencefactor = influencefactor + 0.1
                else
                    flag = FALSE
            }
        }
        
        
        # for all E-genes calculate probabilities
        for (i in seq(1, nrEgenes)) {
            mean.neg = mean(neg.sel[i,])
            mean.pos = mean(pos.sel[i,])
            
            std.neg  = sd(neg.sel[i,])            
            std.pos  = sd(pos.sel[i,])
            std.neg = sd(neg.sel[i,])
            std.pos = sd(pos.sel[i,])
            if ((std.neg+std.pos) < mean.pos-mean.neg) {
                x <- ((mean.pos-mean.neg)-(std.neg+std.pos))/2
                std.neg <- std.neg+x
                std.pos <- std.pos+x
            }
            
            # for all arrays
            for (j in seq(1, nrArrays)) {
                
                if (dat.sel[i,j] <= mean.neg)
                    prob.influenced[i,j] = 1
                else if (dat.sel[i,j] >= mean.pos)
                    prob.influenced[i,j] = 0
                else {                
                    prob.influenced[i,j] <- (dnorm(dat.sel[i,j], mean.neg, std.neg) / 
                        (dnorm(dat.sel[i,j], mean.neg, std.neg) + 
                        dnorm(dat.sel[i,j], mean.pos, std.pos))) * influencefactor
                }
                
                if (prob.influenced[i,j] > 0.9999) 
                    prob.influenced[i,j] = 0.9999
                if (prob.influenced[i,j] < 0.0001)
                    prob.influenced[i,j] = 0.0001
                
            }
            
        } 
        
        dimnames(prob.influenced)    <- dimnames(dat.sel)
        
        
        result <- list(dat=dat.sel, pos=pos.sel, neg=neg.sel, sel=sel, prob.influenced=prob.influenced, influencefactor=influencefactor)
    }
    
    
    
    if (setting == "onecontrol") {
        if (verbose) 
            cat("discretizing with respect to one control\n")
        if (!is.null(pos.control)) {
            if (class(pos.control) == "matrix") {
                W <- pos.control
                M <- D
            }
            if (class(pos.control) %in% c("integer", "numeric")) {
                W <- D[, pos.control]
                M <- D[, -pos.control]
            }
        }
        if (!is.null(neg.control)) {
            if (class(neg.control) == "matrix") {
                W <- neg.control
                M <- D
            }
            if (class(neg.control) %in% c("integer", "numeric")) {
                W <- D[, neg.control]
                M <- D[, -neg.control]
            }
        }
        Wecdf <- apply(W, 1, ecdf)
        Mp <- matrix(0, ncol = ncol(M), nrow = nrow(M))
        for (i in 1:nrow(W)) {
            Pi <- Wecdf[[i]](M[i, ])
            Mp[i, ] <- ifelse(Pi <= 0.5, Pi, 1 - Pi)
        }
        Mt <- (Mp <= empPval) * 1
        dimnames(Mt) <- dimnames(M)
        result <- list(dat = Mt)
    }
    
    
    return(result)
    
}
