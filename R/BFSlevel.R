BFSlevel <- function(g,verbose=TRUE){

    ## 1. I should remove loops in g! that could cause problems ....
    ## 2. what if the graph is not one connected component?
    
    level <- rep(NA,length(nodes(g)))
    names(level) <- nodes(g)
    notvisited   <- nodes(g)

    ## level 1
    i <- 1
    thislevel <- nodes(g)[sapply(edgeL(g),function(x) length(x$edges))==0]    
    level[thislevel] = i
    notvisited <- setdiff(notvisited,thislevel)

    ## iteratively build next levels
    while(length(notvisited)>0){
    cat(i,"\n")
        i <- i + 1
        parents <- unique(unlist(sapply(thislevel, function(x) inEdges(x,g))))
        thislevel <- intersect(notvisited,parents)
        level[thislevel] = i
        notvisited <- setdiff(notvisited,thislevel)    
    }
    if(verbose) cat("Generalized hierarchy of",length(nodes(g)),"nodes in", i , "levels\n")
    return(level)
}
