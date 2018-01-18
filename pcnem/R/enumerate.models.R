enumerate.models <- function(x,name=NULL, trans.close=TRUE, verbose=TRUE) {

if (length(x) == 1) {
            n <- as.numeric(x)
	    if(is.null(name))
            	name <- letters[1:n]
        } else {
           n <- length(x)
           name <- x
        }

#------------------
# Sanity checks    

 if (n==1) stop("nem> choose n>1!")
 if (n>5)  stop("nem> exhaustive enumeration not feasible with more than 5 perturbed genes")            
 if (n==5) cat ("nem> this will take a while ... \n") 

#------------------

  bc <- bincombinations(n*(n-1))
  fkt1 <- function(x,n,name) {
    M <- diag(n)
    M[which(M==0)]<-x
    dimnames(M) <- list(name,name)	
    if(trans.close)    
    	M <- transitive.closure(M,mat=TRUE,loops=TRUE)    
    return(list(M))
  }
  
  models <- apply(bc,1,fkt1,n,name) 
  models <- unique(matrix(unlist(models),ncol=n*n,byrow=TRUE))
  
  fkt2 <- function(x,n,name){
     M <- matrix(x,n)
     dimnames(M) <- list(name,name)
     return(list(M))
  }
  models <- unlist(apply(models,1,fkt2,n,name),recursive=FALSE)
    
  if (verbose) cat("Generated",length(models),"unique models ( out of", 2^(n*(n-1)), ")\n")

  return(models)
}
