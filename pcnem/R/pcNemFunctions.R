# Supplementary functions for inferring networks using adaptive simulated annealing
#
# Author: Sumana Srivatsa
###############################################################################

##############################
# SUPPLEMENTARY FUNCTIONS #
##############################
#' @title Supplementary functions for adaptive simulated annealing
#'
#' @description The functions are from the paper titled 'Partition MCMC for Inference 
#' on Acyclic Digraphs' by Jack Kuipers & Giusi Moffa which is further based on the code 
#' from the Dortmund course programmed by Miriam Lohr

#### Add citation here ####

### calculation of the first ancestor matrix:
ancestor <- function(incidence){
  incidence1 <- incidence
  incidence2 <- incidence
  k <- 1
  while (k < nrow(incidence)){
    incidence1 <- incidence1%*%incidence
    incidence2 <- incidence2 + incidence1
    k <-k+1
  }
  incidence2[which(incidence2[,]>0)] <- 1
  return(t(incidence2))
}

top_order <- function(incidence){
  n <- nrow(incidence)
  #print(n)
  #browser()
  Order <- numeric(n)
  fan_in <- numeric(n)
  no_fan_in <- numeric(0)
  m <- 1
  for (p in 1:n){                                       # number of parent nodes at the beginning
    fan_in[p] <- sum(incidence[,p])
  }
  no_fan_in <- which(fan_in==0)
  while (length(which(Order==0))>0){                    # as long as there is a node without an order
    fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
    no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
    Order[m] <- no_fan_in[1]
    no_fan_in <- no_fan_in[-1]
    m <- m+1
  }
  return(Order)
}

### assign the topological order of the descendants of the child
des_top_order <- function(incidence, ancest1,child){
  n <- nrow(incidence)
  #print(n)
  #browser()
  top <- top_order(incidence)
  position_child <- which(top==child)
  top_all_after <- top[position_child:n]                # top. order without the "first" nodes
  desc <- which(ancest1[,child]==1)                     # descendants of the child
  inter_step <- c(child,desc,top_all_after)
  des_top <- inter_step[which(duplicated(inter_step))]
  return(des_top)
}

# This function samples an element from a vector properly

propersample <- function(x){if(length(x)==1) x else sample(x,1)} 

################################################################################
# Computation of propagation matrix and scores

############# Functions for propagation matrix #############

#' @title Path count
#' 
#' @description a function to calculate the number of paths between two vertices
#' 
#' @param u source node
#' @param t target node
#' @param phi the S-gene model
#' 
#' @return number of paths (direct + indirect) between two nodes
#' 
#' @author Sumana Srivatsa

PathNumber <- function(u,t,phi){
  if(u == t){
    return(1)
  }else{
    u.children = names(which(phi[u,]==1))
    if(length(u.children)==0){
      return(0)
    }else{
      npaths = sum(unlist(lapply(u.children, function(c){PathNumber(c,t,phi)})))
      return(npaths)
    }
  }
}


#' @title Path count matrix and perturbation matrix
#' 
#' @description a function to calculate the number of paths between all the nodes in phi and
#' compute the propagated perturbation probabilities in each experiment
#' 
#' @param cmap the complement of the knockout map i.e. (1-map)
#' @param phi the S-gene model
#' 
#' @return the path count matrix and the propagation matrix
#' 
#' @author Sumana Srivatsa

getperturb.matrices<-function(cmap,phi,control){
  Sgenes <- colnames(phi)
  N <- length(Sgenes)
  phi2 <- matrix(NA,nrow=nrow(cmap),ncol=ncol(cmap)) # Initalize final matrix
  if(any(diag(phi) == 1)){ diag(phi) <- 0} # no loops
  PathMat <- matrix(NA,nrow=N,ncol=N) # Number of paths from a node(rows) to target node(columns)
  dimnames(PathMat) <- dimnames(phi)
  for(i in 1:N){
    for(j in 1:N){
      PathMat[i,j] <- PathNumber(Sgenes[i],Sgenes[j],phi)
    }
  }
  for(k in 1:nrow(cmap)){
    for(n in 1:ncol(cmap)){
      g <- which(PathMat[,n] != 0)
      phi2[k,n] <- 1-prod((cmap[k,g]^PathMat[g,n]))
    }
  }
  colnames(phi2) <- colnames(phi)
  rownames(phi2) <- rownames(cmap)
  
  if(control$selEGenes.method == "regularization"){
    Phi2 = cbind(phi2, double(nrow(phi2)))
    colnames(Phi2)[ncol(Phi2)] = "null"
    colnames(Phi2)[1:(ncol(phi2))]<-Sgenes
  }
  else 
    Phi2 = phi2
  return(list(PropMat=Phi2, PathMat=PathMat))
}


#' @title Propagation matrix for only a subset of nodes
#' 
#' @description a function to calculate the number of paths between rescore nodes and
#' compute the propagated perturbation probabilities in each experiment
#' 
#' @param phi the newly sampled S-gene model
#' @param cmap the complement of the knockout map i.e. (1-map)
#' @param phi2 the propagation matrix corresponding to the previous S-gene model
#' @param PathMat the path count matrix corresponding to the previous S-gene model
#' @param rescorenodes the nodes in the newly sampled S-gene model that are affected
#' by edge modification to the previous S-gene model
#' 
#' @return the path count matrix and the propagation matrix for the newly sampled S-gene model
#' 
#' @author Sumana Srivatsa

probmatrix.rescorenodes <- function(phi,phi2,cmap,PathMat,rescorenodes){
  Sgenes <- colnames(phi)
  N <- length(Sgenes)
  if(any(diag(phi) == 1)){ diag(phi) <- 0} # no loops
  for(i in 1:N){
    for(j in rescorenodes){
      PathMat[i,j] <- PathNumber(Sgenes[i],Sgenes[j],phi)
    }
  }
  for(k in 1:nrow(cmap)){
    for(n in rescorenodes){
      g <- which(PathMat[,n] != 0)
      phi2[k,n] <- 1-prod((cmap[k,g]^PathMat[g,n]))
    }
  }
  
  return(list(PropMat = phi2,PathMat = PathMat))
}

# Functions for DAG scoring 

#' @title Log likelihood score for a S-gene model
#' 
#' @description a function to calculate the log likelihood score for a S-gene model
#' 
#' @param phi2 the propagation matrix corresponding to the S-gene model
#' @param D1 Effects observed in the data
#' @param D0 Effects not observed in the data
#' @param control the control parameters
#' 
#' @return list of log likelihood and the marginalised likelihood matrix 
#' 
#' @author Sumana Srivatsa

DAGscore <- function(phi2,D1,D0,control,alpha,beta){
  para=NULL
  
  L <- matrix(NA,nrow = nrow(D1),ncol=ncol(phi2))
  for(l in 1:nrow(D1)){
    L[l,] <- apply(D1[l,]*(alpha*(1-phi2) + (1-beta)*phi2) +
                     D0[l,]* (beta*phi2 + (1-alpha)*(1-phi2)),2,prod)
  }
  if(!is.null(control$Pe))
    LP <- L*control$Pe
  else
    LP <- L	
  LLperGene = log(rowSums(LP))
  s  <- sum(LLperGene)
  return(list(mLL=s,LMat=L))
}


#' @title Updating the log likelihood using rescore nodes
#' 
#' @description a function to calculate the log likelihood score for newly sampled 
#' S-gene model by updating the old score using the rescore nodes
#' 
#' @param phi2 the propagation matrix corresponding to the S-gene model
#' @param D1 Effects observed in the data
#' @param D0 Effects not observed in the data
#' @param control the control parameters
#' @param L  marginalised likelihood matrix for the previous S-gene model
#' @param rescorenodes the nodes in the newly sampled S-gene model that are affected
#' by edge modification to the previous S-gene model
#' 
#' @return list of log likelihood and the marginalised likelihood matrix for the
#'  newly sampled S-gene model
#' 
#' @author Sumana Srivatsa
# B. calculate the Log likelihood score for a network only for the rescoring nodes

DAG.rescore <- function(phi2,D1,D0,control,L,rescorenodes,alpha,beta){
  
  for(l in 1:nrow(D1)){
    for(node in rescorenodes){
      L[l,node] <- prod(D1[l,]*(alpha*(1-phi2[,node]) + (1-beta)*phi2[,node]) +
                          D0[l,]* (beta*phi2[,node] + (1-alpha)*(1-phi2[,node])))
    }
    if(!is.null(control$Pe))
      LP <- L*control$Pe
    else
      LP <- L	
    LLperGene = log(rowSums(LP))
    s  <- sum(LLperGene)
  }
  return(list(mLL=s,LMat=L))
}

rec.sample <- function(current_value,proposed_value,Sigma,i){
  if(abs(proposed_value) > 2){
    proposed_value <- current_value + mvrnorm(n = 1, rep(0, 2), Sigma)[i]
    return(rec.sample(current_value,proposed_value,Sigma,i))
  }else
    return(proposed_value)
}
