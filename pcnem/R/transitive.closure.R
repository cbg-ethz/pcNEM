transitive.closure <- function(g,mat=FALSE,loops=TRUE){

    if (!(class(g)%in%c("graphNEL","matrix"))) stop("Input must be either graphNEL object or adjacency matrix")
    g <- as(g, "matrix")
    
    #-- adjacency matrix
#     if (class(g)=="matrix"){
		n <- ncol(g)
		matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
		while(i < pow) {
			z <- z %*% x
			y <- y+z
			i <- i+1
		}
		return(y)
		}
	
		h <- matExpIterativ(g,n)
		h <- (h>0)*1   
		dimnames(h) <- dimnames(g)
		if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
		if (!mat) h <- as(h,"graphNEL")	
#     }

# #     -- graphNEL object
#     if (class(g)=="graphNEL"){
#         tc <- RBGL::transitive.closure(g)    
#         if (loops) tc$edges <- unique(cbind(tc$edges,rbind(tc$nodes,tc$nodes)),MARGIN=2)
# 
#         h <- ftM2graphNEL(ft=t(tc$edges),V=tc$nodes)
#         if (mat) h <- as(h, "matrix")
#     } 
       
    return(h)
}
