nem.discretize <- function(D,neg.control=NULL,pos.control=NULL,nfold=2,cutoff=0:10/10, pCounts=20, empPval=.05, verbose=TRUE){

#-------------------------
# sanity checks           

if (is.null(neg.control) & is.null(pos.control))                stop("\nnem> provide at least one control")
if (class(neg.control)=="matrix") if (nrow(neg.control)!=nrow(D))  stop("\nnem> control and data must have the same number of rows")
if (class(pos.control)=="matrix") if (nrow(pos.control)!=nrow(D))  stop("\nnem> control and data must have the same number of rows")
if (class(neg.control)%in%c("integer","numeric") & class(pos.control)%in%c("integer","numeric") & !all(c(neg.control,pos.control))%in%1:ncol(D))    stop("\nnem>controls not in data 'D'")
if (class(neg.control)%in%c("integer","numeric") & class(pos.control)%in%c("integer","numeric") & any(neg.control %in% pos.control))                stop("\nnem>negative and positive controls overlap")

if (!is.null(neg.control) & !is.null(pos.control)) { setting <- "twocontrols" } else { setting <- "onecontrol"}


#-------------------------
# two controls scenario   
# (Markowetz et al, 2005) 

if (setting=="twocontrols"){
if (verbose) cat("discretizing with respect to POS and NEG controls\n")
if (class(neg.control)=="matrix")                   neg <- neg.control 
if (class(neg.control)%in%c("integer","numeric"))   neg <- D[,neg.control] 
if (class(pos.control)=="matrix")                   pos <- pos.control
if (class(pos.control)%in%c("integer","numeric"))   pos <- D[,pos.control]

if (class(neg.control)%in%c("integer","numeric")){d <- neg.control}else{d<-NULL}
if (class(pos.control)%in%c("integer","numeric")) d <- c(d,pos.control)
if (!is.null(d)) dat <- D[,-d]

# select diff - maybe as extra input function??
# should also do downregulation ...
sel <- which(exp(rowMeans(pos) - rowMeans(neg)) > nfold)
dat.sel <- dat[sel,]
pos.sel <- pos[sel,]
neg.sel <- neg[sel,]

# count false decisions for different cutoff levels
count.false.decisions <- function(x){
thrsh    <- x*rowMeans(pos.sel) + (1-x)*rowMeans(neg.sel)
pos.disc <- (pos.sel <= thrsh)*1
neg.disc <- (neg.sel <= thrsh)*1        
a        <- round(sum(pos.disc)/length(pos.disc),2)
b        <- 1-round(sum(neg.disc)/length(neg.disc),2)                   
return(c(a,b))
}
false <- sapply(cutoff,count.false.decisions)
dimnames(false) <- list(c("a","b"),as.character(cutoff))

mycutoff <- cutoff[which.min(false[2,])]

# apply chosen cutoff
thrsh  <- mycutoff*rowMeans(pos.sel) + (1-mycutoff)*rowMeans(neg.sel)
dat.disc <- (dat.sel <= thrsh)*1
pos.disc <- (pos.sel <= thrsh)*1
neg.disc <- (neg.sel <= thrsh)*1        
a <-   round((sum(pos.disc)+pCounts)/(length(pos.disc)+pCounts),2)
b <- 1-round( sum(neg.disc)         /(length(neg.disc)+pCounts),2)                   
para <- c(a,b)
names(para) <- c("a","b")

# output
disc <- list(dat=dat.disc,pos=pos.disc,neg=neg.disc,sel=sel,cutoff=false,para=para)

}


#-------------------------
# one control scenario    

if (setting=="onecontrol"){
if (verbose) cat("discretizing with respect to one control\n")

if (!is.null(pos.control)){
    if (class(pos.control)=="matrix"){                   
        W <- pos.control
        M <- D
        }
    if (class(pos.control)%in%c("integer","numeric")){
        W <- D[,pos.control]
        M <- D[,-pos.control]
        }
    }
if (!is.null(neg.control)){
    if (class(neg.control)=="matrix"){
        W <- neg.control
        M <- D
        }
    if (class(neg.control)%in%c("integer","numeric")){
        W <- D[,neg.control]
        M <- D[,-neg.control]
        }
    }


# empirical distr. function
Wecdf <- apply(W,1,ecdf) 
Mp <- matrix(0,ncol=ncol(M),nrow=nrow(M))   
for (i in 1:nrow(W)){
    Pi <- Wecdf[[i]](M[i,])    
    Mp[i,] <- ifelse(Pi<=.5,Pi,1-Pi)   
} 
Mt <- (Mp <= empPval)*1   
dimnames(Mt) <- dimnames(M)

# output
disc <- list(dat=Mt) 

}


#-------------------------
# output                  
return(disc)

}
