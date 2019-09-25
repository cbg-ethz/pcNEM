# Function for learning the closest DAG to consensus network after using bootstrap
#
# Author: Sumana Srivatsa
###############################################################################

#' @title Learn the closest DAG to consensus network after using bootstrap
#'
#' @description a function which performs adaptive simulated annealing to find
#' the nearest DAG to the consensus directed graph after running bootstrap. 
#' NOTE : This function is different from the AdaSimulatedAnnealing which infers 
#' the MLE network.
#' The code for adaptive simulated annealing
#' has been developed using the code from the paper titled 'Partition MCMC for
#' Inference on Acyclic Digraphs' by Jack Kuipers & Giusi Moffa which is further
#' based on the code from the Dortmund course programmed by Miriam Lohr.
#' The code has been modified to match it to the NEMs framework.
#' The scoring function has been changed to compute the likelihood according to our model.
#' Further, the update steps have been modified.
#' (https://www.statistik.tu-dortmund.de/fileadmin/user_upload/Lehrstuehle/Genetik/GN10/structureMCMC.txt).
#'
#' @param n number of S-genes
#' @param phi the starting S-gene model
#' @param D the data set. The rows refer to the E-genes and the columns correspond
#' to the knockdown experiment.
#' @param control the control parameters
#' @param verbose print output
#'
#' @return a list containing the graphs searched, their likelihoods, acceptance rates, temperatures, distances
#' after every 'stepsave' number of steps, the best graph and it's corresponding likelihood and distance,
#' the minimum number of steps required to reach the maximum likelihood, and the
#' transformation exponent for making the ideal acceptance rate symmetrical
#'
#' @author Sumana Srivatsa


SearchNearestDAG <- function(n,phi,D,control,consensus){
  
  moveprobs <- control$moveprobs
  Sigma <- control$sigma
  iterations <- control$iterations
  stepsave <- control$stepsave
  fan.in <- (n-1)
  revallowed <- control$revallowed
  if(is.null(control$AcceptRate)){ # ideal acceptance rate
    AcceptRate <- 1/n
  }else{
    AcceptRate <- control$AcceptRate
  }
  Temp <- control$Temp
  AdaptRate <- control$AdaptRate
  
  L1 <- list() # stores the adjecency matrices
  L2 <- list() # stores the temperatures
  L3 <- list() # stores the acceptance rates after 'stepsave' iterations
  L4 <- list() # stores the distances after 'stepsave' iterations
  
  transformAR <- log(0.5)/log(AcceptRate) # to make the ideal acceptance rate symmetric
  
  # Initializing parameters
  Sgenes <- colnames(phi)
  cmap <- 1-control$map # complementary probabilities
  dimnames(cmap) <- dimnames(control$map)
  nexp <- rownames(cmap) # number of experiments
  D1 = sapply(nexp, function(s) rowSums(D[, colnames(D) == s, drop = FALSE])) # Signals present in the data
  D0 = sapply(nexp, function(s) sum(colnames(D) == s)) - D1 # Signals absent in the data
  
  # For marginalizing the effects we use uniform priors across all Sgenes
  if (is.null(control$Pe)) {
    control$Pe <- matrix(1/n, nrow = nrow(D1), ncol = n)
    colnames(control$Pe) <- Sgenes
  }
  if (control$selEGenes.method == "regularization" && ncol(control$Pe) == n) {
    control$Pe = cbind(control$Pe, double(nrow(D1)))
    control$Pe[, ncol(control$Pe)] = control$delta/n
    control$Pe = control$Pe/rowSums(control$Pe)
  }
  
  initPertMats <- getperturb.matrices(cmap,phi,control) # To get initial propagation matrix and path count matrix
  PropMat <- initPertMats$PropMat # Propagation matrix
  PathMat <- initPertMats$PathMat # Path count matrix
  Alpha <- control$para[1]
  Beta <- control$para[2]
  DAGParams <- DAGscore(PropMat,D1,D0,control,Alpha,Beta)
  DAGLMat <- DAGParams$LMat
  currentDistance <- shd(as(phi,'graphNEL'),as(consensus,"graphNEL"))
  
  L1 <- c(L1,list(phi)) # starting adjacency matrix
  L2 <- c(L2,Temp) # starting temperature
  L3 <- c(L3,0) # starting acceptance rate
  L4 <- c(L4,currentDistance)
    
  bestScore <- DAGParams$mLL
  bestGraph <- phi
  bestPropMat <- PropMat
  minSteps <- 0
  bestDistance <- currentDistance
  
  naccepts_DAG <- 0
  
  # first ancestor matrix
  ancest1 <- ancestor(phi)
  
  ####### ... the number of neighbour graphs/proposal probability for the FIRST graph
  ### 1.) number of neighbour graphs obtained by edge deletions
  num_deletion <- sum(phi)
  
  emptymatrix <- matrix(numeric(n*n),nrow=n)
  fullmatrix <- matrix(rep(1,n*n),nrow=n)
  Ematrix <- diag(1,n,n)
  
  ### 2.) number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add <- which(fullmatrix - Ematrix - phi - ancest1 >0)
  add <- emptymatrix
  add[inter_add] <- 1
  add[,which(colSums(phi)>fan.in-1)] <- 0
  num_addition <- sum(add)
  
  ### 3.) number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
  inter_rev <- which(phi - t(t(phi)%*% ancest1)==1)
  re <- emptymatrix
  re[inter_rev] <- 1
  re[which(colSums(phi)>fan.in-1),] <- 0 # this has to be this way around
  num_reversal <- sum(re)
  
  ##### total number of neighbour graphs:
  currentnbhoodnorevs <- sum(num_deletion,num_addition)+1
  currentnbhood <- currentnbhoodnorevs+num_reversal
  
  for (z in 1:iterations){
    
    ### sample one of the three single edge operations
    if(revallowed==1){
      operation<-sample.int(4,1,prob=c(num_reversal,num_deletion,num_addition,1)) # sample the type of move including staying still
    }else{
      operation<-sample.int(3,1,prob=c(num_deletion,num_addition,1))+1 # sample the type of move including staying still
    }
    
    # 1 is edge reversal, 2 is deletion and 3 is additon. 4 represents choosing the current DAG
    if(operation<4){ # if we don't stay still
      
      #### shifting of the phi matrix
      new_phi <- phi
      
      # creating a matrix with dimensions of the phi matrix and all entries zero except for the entry of the chosen edge
      help_matrix <- emptymatrix
      
      if (operation==1){    # if edge reversal was sampled
        new_edge <- propersample(which(re==1))      # sample one of the existing edges where a reversal leads to a valid graph
        new_phi[new_edge] <- 0             # delete it
        help_matrix[new_edge] <- 1               # an only a "1" at the entry of the new (reversed) edge
        new_phi <- new_phi + t(help_matrix) # sum the deleted matrix and the "help-matrix"
      }
      if (operation==2){              # if edge deletion was sampled
        new_edge <- propersample(which(phi>0)) # sample one of the existing edges
        new_phi[new_edge] <- 0            # and delete it
        help_matrix[new_edge] <- 1
      }
      if (operation==3){     # if edge addition was sampled
        new_edge <- propersample(which(add==1)) # sample one of the existing edges where a addition leads to a valid graph
        new_phi[new_edge] <- 1             # and add it
        help_matrix[new_edge] <- 1
      }
      
      ### Updating the ancestor matrix
      
      # numbers of the nodes that belong to the shifted egde
      parent <- which(rowSums(help_matrix)==1)
      child <- which(colSums(help_matrix)==1)
      
      ### updating the ancestor matrix (after edge reversal)
      ## edge deletion
      ancestor_new <- ancest1
      if (operation==1){
        ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants
        top_name <- des_top_order(new_phi, ancest1, child)
        for (d in top_name){
          for(g in which(new_phi[,d]==1)) {
            ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
          }
        }
        
        anc_parent <- which(ancestor_new[child,]==1)          # ancestors of the new parent
        des_child <- which(ancestor_new[,parent]==1)          # descendants of the child
        ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
      }
      
      ### updating the ancestor matrix (after edge deletion)
      if (operation==2){
        ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           #
        top_name <- des_top_order(new_phi, ancest1, child)
        for (d in top_name){
          for(g in which(new_phi[,d]==1)) {
            ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
          }
        }
      }
      
      # updating the ancestor matrix (after edge addition)
      if (operation==3){
        anc_parent <- which(ancest1[parent,]==1)             # ancestors of the new parent
        des_child <- which(ancest1[,child]==1)               # descendants of the child
        ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
      }
      
      ####### ... the number of neighbour graphs/proposal probability for the proposed graph
      ### 1.) number of neighbour graphs obtained by edge deletions
      num_deletion_new <- sum(new_phi)
      
      ### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
      inter_add.new <- which(fullmatrix - Ematrix - new_phi - ancestor_new >0)
      add.new <- emptymatrix
      add.new[inter_add.new] <- 1
      add.new[,which(colSums(new_phi)>fan.in-1)] <- 0
      num_addition_new <- sum(add.new)
      
      ### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
      inter_rev.new<- which(new_phi - t(t(new_phi)%*% ancestor_new)==1)
      re.new <- emptymatrix
      re.new[inter_rev.new] <- 1
      re.new[which(colSums(new_phi)>fan.in-1),] <- 0 # this has to be this way around!
      num_reversal_new <- sum(re.new)
      
      ##### total number of neighbour graphs:
      proposednbhoodnorevs<-sum(num_deletion_new, num_addition_new) + 1
      proposednbhood <- proposednbhoodnorevs + num_reversal_new
      
      
      if (operation==1){                      # if single edge operation was an edge reversal
        rescorenodes<-unique(c(which(ancestor_new[,child]==1),child,parent))
      }
      if (operation==2){                      # if single edge operation was an edge deletion
        rescorenodes <- which(ancest1[,parent]==1)
      }
      if (operation==3){                      # if single edge operation was an edge deletion
        rescorenodes <- which(ancestor_new[,parent]==1)
      }
      
      proposedDistance <- shd(as(new_phi,'graphNEL'),as(consensus,"graphNEL"))
      
      if(revallowed==1){
        scoreratio<-exp((proposedDistance-currentDistance)/Temp)*(currentnbhood/proposednbhood) #acceptance probability
      } else{
        scoreratio<-exp((proposedDistance-currentDistance)/Temp)*(currentnbhoodnorevs/proposednbhoodnorevs) #acceptance probability
      }
      
      #if(is.na(proposedDAGlogscore)){browser()}
      if(proposedDistance <= currentDistance || runif(1) < scoreratio){ #Move accepted then set the current order and scores to the proposal
        
        phi <- new_phi
        currentDistance <- proposedDistance
        ancest1 <- ancestor_new
        currentnbhood <- proposednbhood
        currentnbhoodnorevs <- proposednbhoodnorevs
        num_deletion <- num_deletion_new
        num_addition <- num_addition_new
        num_reversal <- num_reversal_new
        add <- add.new
        re <- re.new
        naccepts_DAG <- naccepts_DAG + 1
        
        if(bestDistance > currentDistance){
          proposedPertMats <- probmatrix.rescorenodes(phi,PropMat,cmap,PathMat,rescorenodes)
          PropMat <- proposedPertMats$PropMat
          PathMat <- proposedPertMats$PathMat
          proposedDAGlogParams <- DAG.rescore(PropMat,D1,D0,control,DAGLMat,rescorenodes,Alpha,Beta)
          DAGLMat <- proposedDAGlogParams$LMat
          
          bestScore <- proposedDAGlogParams$mLL
          bestGraph <- phi
          minSteps <- z
          bestPropMat <- PropMat
          bestDistance <- currentDistance
        }
        
        if(bestDistance == currentDistance){
          proposedPertMats <- probmatrix.rescorenodes(phi,PropMat,cmap,PathMat,rescorenodes)
          proposedPropMat <- proposedPertMats$PropMat
          proposedPathMat <- proposedPertMats$PathMat
          proposedDAGlogParams <- DAG.rescore(proposedPropMat,D1,D0,control,DAGLMat,rescorenodes,Alpha,Beta)
          proposedDAGlogscore <- proposedDAGlogParams$mLL
          proposedDAGLMat <- proposedDAGlogParams$LMat
          
          if(bestScore < proposedDAGlogscore){
            PropMat <- proposedPropMat
            PathMat <- proposedPathMat
            DAGLMat <- proposedDAGLMat
            
            bestScore <- proposedDAGlogscore
            bestGraph <- phi
            minSteps <- z
            bestPropMat <- PropMat
            bestDistance <- currentDistance
          }
        }
      } # end of acceptance condition
    }  # end of staying still loop
    if(naccepts_DAG == stepsave){
      currentAR <- naccepts_DAG/stepsave  # current acceptance rate
      Temp <- Temp * exp((0.5-currentAR^transformAR) * AdaptRate) # adapting the temperature according to the acceptance rate
      # Updating all
      L1 <- c(L1,list(phi))
      L2 <- c(L2,Temp)
      L3 <- c(L3,currentAR)
      L4 <- c(L4,currentDistance)
      naccepts_DAG = 0
    }
  }
  return(list(graphs=L1,Temp=L2,AcceptRate=L3,Distances=L4,
              maxLLscore=bestScore,maxDAG=bestGraph,bestDistance=bestDistance,
              minIter=minSteps,transformAR=transformAR))
}
