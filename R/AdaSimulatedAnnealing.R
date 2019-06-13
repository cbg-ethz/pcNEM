# Function for inferring networks using adaptive simulated annealing
#
# Author: Sumana Srivatsa
###############################################################################

#' @title Adaptive simulated annealing for probabilistic combinatorial knockdowns
#'
#' @description a function which performs adaptive simulated annealing to find
#' the network with maximum likelihood. The code for adaptive simulated annealing
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
#' @return a list containing the graphs searched, their likelihoods, acceptance rates, temperatures, typeI and
#' typeII errors after every 'stepsave' number of steps, the best graph and it's corresponding likelihood,
#' the best estimates of errors, the minimum number of steps required to reach the maximum likelihood, and the
#' transformation exponent for making the ideal acceptance rate symmetrical
#'
#' @author Sumana Srivatsa
#'

#### Add citation here ####

# Adaptive simulated annealing to find the maximum likelihood for probabilistic combinatorial knockdowns
#
# The code has been modified to match it to the NEMs framework. The scoring function has been changed to compute the likelihood according to our model. Further, the update steps have been modified.


##############################
# SIMULATED ANNEALING CODE #
##############################

nem.AdaSA <- function(n,phi,D,control){
    
    # initializations
    chosenmove     <- 1 # starting with structure move
    noiseEst       <- control$noiseEst # if noise parameters need to be inferred
    moveprobs      <- control$moveprobs # moving between the two spaces
    moveprobsNoise <- control$moveprobsNoise # moving between alpha and beta noise para
    Sigma          <- control$sigma # initial noise para covariance matrix
    iterations     <- control$iterations 
    stepsave       <- control$stepsave
    temper         <- control$temper # choose tempering or not
    
    # moves paras
    fan.in     <- (n-1) 
    revallowed <- control$revallowed
    
    # ideal acceptance rate
    if(is.null(control$AcceptRate)){ 
        AcceptRate <- 1/n
    }else{
        AcceptRate <- control$AcceptRate
    }
    
    # start temp and adaptation rate
    Temp      <- control$Temp
    AdaptRate <- control$AdaptRate
    
    L1 <- list() # stores the adjecency matrices
    L2 <- list() # stores the log likelihood score of the DAGs
    L3 <- list() # stores the temperatures
    L4 <- list() # stores the acceptance rates after 'stepsave' iterations
    L5 <- list() # stores the alpha values
    L6 <- list() # stores the beta values
    
    transformAR <- log(0.5)/log(AcceptRate) # to make the ideal acceptance rate symmetrical i.e. AcceptRate^transformAR = 0.5
    print(paste0("The exponent for the acceptance rate is ",transformAR))
    
    # Initializing nem parameters
    Sgenes         <- colnames(phi)
    cmap           <- 1-control$map # complementary probabilities
    dimnames(cmap) <- dimnames(control$map)
    nexp           <- rownames(cmap) # number of experiments
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
    
    
    # Starting SA
    initPertMats       <- getperturb.matrices(cmap,phi,control) # To get initial propagation matrix and path count matrix
    PropMat            <- initPertMats$PropMat # Propagation matrix
    PathMat            <- initPertMats$PathMat # Path count matrix
    currentAlpha       <- control$para[1]
    currentBeta        <- control$para[2]
    currentDAGParams   <- DAGscore(PropMat,D1,D0,control,currentAlpha,currentBeta)
    currentDAGlogscore <- currentDAGParams$mLL
    currentDAGLMat     <- currentDAGParams$LMat
    
    L1 <- c(L1,list(phi)) # first adjacency matrix
    L2 <- c(L2,currentDAGlogscore) # first DAG score
    L3 <- c(L3,Temp) # starting temperature
    L4 <- c(L4,0) # starting acceptance rate
    L5 <- c(L5,currentAlpha) # starting alpha
    L6 <- c(L6,currentBeta) # starting beta
    
    # Starting noise parameter chains
    AlphaChain <- NULL
    BetaChain  <- NULL
    
    # assigning starting DAG as best DAG
    bestScore   <- currentDAGlogscore
    bestGraph   <- phi
    bestPropMat <- PropMat
    bestAlpha   <- currentAlpha
    bestBeta    <- currentBeta
    minSteps    <- 0
    
    # initialize number of steps for each parameter to be inferred
    nsteps_alpha <- 0
    nsteps_beta  <- 0
    nsteps_DAG   <- 0
    naccepts_DAG <- 0
    #nsampled_alpha <- 0
    #nsampled_beta <- 0
    #adaptation_sigma <- 2*stepsave
    
    # first ancestor matrix
    ancest1 <- ancestor(phi)
    
    ####### ... the number of neighbour graphs/proposal probability for the FIRST graph
    ### 1.) number of neighbour graphs obtained by edge deletions
    num_deletion <- sum(phi)
    
    emptymatrix  <- matrix(numeric(n*n),nrow=n)
    fullmatrix   <- matrix(rep(1,n*n),nrow=n)
    Ematrix      <- diag(1,n,n)
    
    ### 2.) number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
    inter_add      <- which(fullmatrix - Ematrix - phi - ancest1 >0)
    add            <- emptymatrix
    add[inter_add] <- 1
    add[,which(colSums(phi)>fan.in-1)] <- 0
    num_addition   <- sum(add)
    
    ### 3.) number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
    inter_rev     <- which(phi - t(t(phi)%*% ancest1)==1)
    re            <- emptymatrix
    re[inter_rev] <- 1
    re[which(colSums(phi)>fan.in-1),] <- 0 # this has to be this way around
    num_reversal  <- sum(re)
    
    ##### total number of neighbour graphs:
    currentnbhoodnorevs <- sum(num_deletion,num_addition)+1
    currentnbhood       <- currentnbhoodnorevs+num_reversal
    
    if(iterations > 1){ # Iterations can be set to 0 to only infer noise parameters
      # looping through all iterations
      for (z in 1:iterations){
        if(noiseEst == TRUE){ # if we also want to infer noise para then we sample the space to explore
          chosenmove <- sample.int(2,1,prob=moveprobs) 
        }
        switch(as.character(chosenmove),
               "1"={ # searching for DAG
                 
                 nsteps_DAG <- nsteps_DAG + 1 # updating the iterations in DAG space
                 
                 ### sample one of the three single edge operations
                 if(revallowed==1){
                   operation <- sample.int(4,1,prob=c(num_reversal,num_deletion,num_addition,1)) # sample the type of move including reversal and staying still
                 }else{
                   operation <- sample.int(3,1,prob=c(num_deletion,num_addition,1))+1 # sample the type of move including staying still but without edge reversal
                 }
                 
                 # 1 is edge reversal, 2 is deletion and 3 is additon. 4 represents choosing the current DAG
                 if(operation < 4){ # for all other moves except staying in the current DAG
                   
                   #### shifting of the phi matrix
                   new_phi <- phi
                   
                   # creating a matrix with dimensions of the phi matrix and all entries zero except for the entry of the chosen edge
                   help_matrix <- emptymatrix
                   
                   if (operation==1){    # if edge reversal was sampled
                     new_edge          <- propersample(which(re==1)) # sample one of the existing edges where a reversal leads to a valid graph
                     new_phi[new_edge] <- 0   # delete it
                     help_matrix[new_edge] <- 1  # an only a "1" at the entry of the new (reversed) edge
                     new_phi  <- new_phi + t(help_matrix) # sum the deleted matrix and the "help-matrix"
                   }
                   if (operation==2){  # if edge deletion was sampled
                     new_edge          <- propersample(which(phi>0)) # sample one of the existing edges
                     new_phi[new_edge] <- 0  # and delete it
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
                   child  <- which(colSums(help_matrix)==1)
                   
                   ### updating the ancestor matrix (after edge reversal)
                   ## edge deletion
                   ancestor_new <- ancest1
                   if (operation == 1){
                     ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants
                     top_name <- des_top_order(new_phi, ancest1, child)
                     for (d in top_name){
                       for(g in which(new_phi[,d]==1)) {
                         ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
                       }
                     }
                     
                     anc_parent <- which(ancestor_new[child,]==1)          # ancestors of the new parent
                     des_child  <- which(ancestor_new[,parent]==1)          # descendants of the child
                     ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
                   }
                   
                   ### updating the ancestor matrix (after edge deletion)
                   if (operation == 2){
                     ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           #
                     top_name <- des_top_order(new_phi, ancest1, child)
                     for (d in top_name){
                       for(g in which(new_phi[,d]==1)) {
                         ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
                       }
                     }
                   }
                   
                   # updating the ancestor matrix (after edge addition)
                   if (operation == 3){
                     anc_parent <- which(ancest1[parent,]==1)             # ancestors of the new parent
                     des_child  <- which(ancest1[,child]==1)               # descendants of the child
                     ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
                   }
                   
                   ####### ... the number of neighbour graphs/proposal probability for the proposed graph
                   ### 1.) number of neighbour graphs obtained by edge deletions
                   num_deletion_new <- sum(new_phi)
                   
                   ### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
                   inter_add.new <- which(fullmatrix - Ematrix - new_phi - ancestor_new >0)
                   add.new       <- emptymatrix
                   add.new[inter_add.new] <- 1
                   add.new[,which(colSums(new_phi)>fan.in-1)] <- 0
                   num_addition_new <- sum(add.new)
                   
                   ### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
                   inter_rev.new <- which(new_phi - t(t(new_phi)%*% ancestor_new)==1)
                   re.new        <- emptymatrix
                   re.new[inter_rev.new] <- 1
                   re.new[which(colSums(new_phi)>fan.in-1),] <- 0 # this has to be this way around!
                   num_reversal_new <- sum(re.new)
                   
                   ##### total number of neighbour graphs:
                   proposednbhoodnorevs <- sum(num_deletion_new, num_addition_new) + 1
                   proposednbhood       <- proposednbhoodnorevs + num_reversal_new
                   
                   
                   if (operation==1){                      # if single edge operation was an edge reversal
                     rescorenodes <- unique(c(which(ancestor_new[,child]==1),child,parent))
                   }
                   if (operation==2){                      # if single edge operation was an edge deletion
                     rescorenodes <- which(ancest1[,parent]==1)
                   }
                   if (operation==3){                      # if single edge operation was an edge deletion
                     rescorenodes <- which(ancestor_new[,parent]==1)
                   }
                   
                   # Updating the nem propagation matrix, path count matrix and log likelihood for the proposed DAG
                   proposedPertMats <- probmatrix.rescorenodes(new_phi,PropMat,cmap,PathMat,rescorenodes)
                   proposedPropMat  <- proposedPertMats$PropMat
                   proposedPathMat  <- proposedPertMats$PathMat
                   proposedDAGlogParams <- DAG.rescore(proposedPropMat,D1,D0,control,currentDAGLMat,rescorenodes,currentAlpha,currentBeta)
                   proposedDAGlogscore  <- proposedDAGlogParams$mLL
                   proposedDAGLMat      <- proposedDAGlogParams$LMat
                   
                   
                   if(revallowed==1){
                     scoreratio <- exp((proposedDAGlogscore-currentDAGlogscore)/Temp)*(currentnbhood/proposednbhood) #acceptance probability
                   } else{
                     scoreratio <- exp((proposedDAGlogscore-currentDAGlogscore)/Temp)*(currentnbhoodnorevs/proposednbhoodnorevs) #acceptance probability
                   }
                   
                   #if(is.na(proposedDAGlogscore)){browser()}
                   # updating current DAG when proposal is accepted
                   if(currentDAGlogscore < proposedDAGlogscore || runif(1) < scoreratio){ #Move accepted then set the current order and scores to the proposal
                     PropMat <- proposedPropMat
                     PathMat <- proposedPathMat
                     phi     <- new_phi
                     currentDAGlogscore <- proposedDAGlogscore
                     currentDAGLMat     <- proposedDAGLMat
                     # updating best DAG and corresponding log likelihood 
                     if(bestScore < currentDAGlogscore){
                       bestScore   <- currentDAGlogscore
                       bestGraph   <- phi
                       minSteps    <- z
                       bestPropMat <- PropMat
                     }
                     # updating ancestor and neighbourhood
                     ancest1 <- ancestor_new
                     currentnbhood <-proposednbhood
                     currentnbhoodnorevs <- proposednbhoodnorevs
                     # updating prob of sampling operations 1,2, or 3
                     num_deletion <- num_deletion_new 
                     num_addition <- num_addition_new
                     num_reversal <- num_reversal_new
                     add <- add.new
                     re  <- re.new
                     # updating number of accepted DAGs
                     naccepts_DAG <- naccepts_DAG + 1
                   }
                 }  # end of staying still condition
                 
                 if(temper == FALSE){ #using adaptive steps and constant rate of cooling
                   if(naccepts_DAG == stepsave){
                     currentAR <- naccepts_DAG/stepsave  # current acceptance rate
                     Temp      <- Temp * exp((-0.5) * AdaptRate) 
                     # Updating all lists
                     L1 <- c(L1,list(phi))
                     L2 <- c(L2,currentDAGlogscore)
                     L3 <- c(L3,Temp)
                     L4 <- c(L4,currentAR)
                     naccepts_DAG <- 0 # restarting number of accepted moves
                   }
                 }else{
                   if(nsteps_DAG == stepsave){ #using the formula in the paper
                     currentAR <- naccepts_DAG/stepsave  # current acceptance rate
                     Temp      <- Temp * exp((0.5-currentAR^transformAR) * AdaptRate) # adapting the temperature according to the acceptance rate
                     # Updating all lists
                     L1 <- c(L1,list(phi))
                     L2 <- c(L2,currentDAGlogscore)
                     L3 <- c(L3,Temp)
                     L4 <- c(L4,currentAR)
                     naccepts_DAG <- 0 # restarting number of accepted moves
                     nsteps_DAG   <- 0 # restarting number of steps
                   }
                 }
               },
               "2"={
                 chosenNoise = sample.int(2,1,prob=moveprobsNoise)
                 #print(chosenNoise)
                 #print(c(nsteps_alpha,nsteps_beta))
                 
                 if(chosenNoise == 1){ # update alpha
                   
                   proposedAlpha = currentAlpha + mvrnorm(n = 1, rep(0, 2), Sigma)[1]
                   #nsampled_alpha <- nsampled_alpha + 1
                   
                   # Ensuring that alpha stays in a range of 0 and 0.5
                   if(proposedAlpha < 0){
                     proposedAlpha = abs(proposedAlpha)
                   }
                   if(proposedAlpha > 1 && proposedAlpha < 2){
                     proposedAlpha = proposedAlpha - 2*(proposedAlpha-1)
                   }
                   if(proposedAlpha > 2){
                     proposedAlpha = proposedAlpha %% 1
                     #proposedAlpha = rec.sample(currentAlpha,proposedAlpha,Sigma,1)
                   }
                   
                   #sampled_alphachain <- c(sampled_alphachain,proposedAlpha)
                   #print(proposedAlpha)
                   
                   proposedDAGlogParams <- DAGscore(PropMat,D1,D0,control,proposedAlpha,currentBeta)
                   proposedDAGlogscore  <- proposedDAGlogParams$mLL
                   proposedDAGLMat      <- proposedDAGlogParams$LMat
                   
                   #acceptance probability of alpha
                   scoreratio <- exp((proposedDAGlogscore-currentDAGlogscore)/Temp) 
                   
                   #if(is.na(proposedDAGlogscore)){browser()}
                   # updating if alpha is accepted
                   if(currentDAGlogscore < proposedDAGlogscore || runif(1) < scoreratio){ #Move accepted then set the current order and scores to the proposal
                     currentDAGlogscore <-proposedDAGlogscore
                     currentDAGLMat     <- proposedDAGLMat
                     currentAlpha       <- proposedAlpha
                     #replacing best alpha with proposed if likelihood is better than current best
                     if(bestScore < currentDAGlogscore){
                       bestScore <- currentDAGlogscore
                       bestAlpha <- proposedAlpha
                       minSteps  <- z
                     }
                     AlphaChain   <- c(AlphaChain,currentAlpha)
                     nsteps_alpha <- nsteps_alpha + 1
                   }else{
                     AlphaChain   <- c(AlphaChain,currentAlpha)
                     nsteps_alpha <- nsteps_alpha + 1
                   }
                   #print(nsteps_alpha)
                 }else if(chosenNoise == 2){ # update beta
                   
                   
                   proposedBeta = currentBeta + mvrnorm(n = 1, rep(0, 2), Sigma)[2]
                   #nsampled_beta <- nsampled_beta + 1
                   
                   # Ensuring that alpha stays in a range of 0 and 0.5
                   if(proposedBeta < 0){
                     proposedBeta = abs(proposedBeta)
                   }
                   
                   if(proposedBeta > 1  && proposedBeta < 2){
                     proposedBeta = proposedBeta - 2*(proposedBeta-1)
                   }
                   if(proposedBeta > 2){
                     proposedBeta = proposedBeta %% 1
                     #proposedBeta = rec.sample(currentBeta,proposedBeta,Sigma,2)
                   }
                   
                   #sampled_betachain <- c(sampled_betachain,proposedBeta)
                   #print(proposedBeta)
                   proposedDAGlogParams <- DAGscore(PropMat,D1,D0,control,currentAlpha, proposedBeta)
                   proposedDAGlogscore  <- proposedDAGlogParams$mLL
                   proposedDAGLMat      <- proposedDAGlogParams$LMat
                   
                   #acceptance probability
                   scoreratio <- exp((proposedDAGlogscore-currentDAGlogscore)/Temp) 
                   
                   #if(is.na(proposedDAGlogscore)){browser()}
                   # accepting proposed beta 
                   if(currentDAGlogscore < proposedDAGlogscore || runif(1) < scoreratio){ #Move accepted then set the current order and scores to the proposal
                     currentDAGlogscore <-proposedDAGlogscore
                     currentDAGLMat     <- proposedDAGLMat
                     currentBeta        <- proposedBeta
                     if(bestScore < currentDAGlogscore){
                       bestScore <- currentDAGlogscore
                       bestBeta  <- proposedBeta
                       minSteps  <- z
                     }
                     BetaChain <- c(BetaChain,currentBeta)
                     nsteps_beta <- nsteps_beta + 1
                   }else{
                     BetaChain <- c(BetaChain,currentBeta)
                     nsteps_beta <- nsteps_beta + 1
                   }
                 }else stop("\nem:AdaSA> sampling error")
                 
                 # Once nsteps of alpha and beta are equal to stepsave, update the covariance matrix and lists 
                 if(nsteps_alpha >= stepsave && nsteps_beta >= stepsave){
                   
                   Sigma <- cov(cbind(AlphaChain[1:stepsave],BetaChain[1:stepsave]),
                                cbind(AlphaChain[1:stepsave],BetaChain[1:stepsave]))
                   
                   L5 <- c(L5,AlphaChain[stepsave])
                   L6 <- c(L6,BetaChain[stepsave])
                   
                   #print(Sigma)
                   #print(nsteps_alpha)
                   #print(nsteps_beta)
                   
                   
                   # Reinitializing the chains and steps 
                   if(nsteps_alpha > stepsave){ # Since the choice between alpha and beta is random
                     #L5 <- c(L5,AlphaChain[stepsave])
                     AlphaChain   <- AlphaChain[(stepsave+1):length(AlphaChain)]
                     nsteps_alpha <- nsteps_alpha - stepsave
                   }else{
                     #L5 <- c(L5,AlphaChain[stepsave])
                     AlphaChain   <- NULL
                     nsteps_alpha <- 0
                   }
                   
                   if(nsteps_beta > stepsave){ # Since the choice between alpha and beta is random
                     #L6 <- c(L6,BetaChain[stepsave])
                     BetaChain   <- BetaChain[(stepsave+1):length(BetaChain)]
                     nsteps_beta <- nsteps_beta - stepsave
                   }else{
                     #L6 <- c(L6,BetaChain[stepsave])
                     BetaChain   <- NULL
                     nsteps_beta <- 0
                   }
                 }
               },
               {# if neither is chosen, we have a problem
                 print('The move sampling has failed!')
               })
      }
    }
    # Again ensuring that alpha and beta stay in the range
    if(bestAlpha > 0.5){
      bestAlpha = 1-bestAlpha
    }
    if(bestBeta > 0.5){
      bestBeta = 1-bestBeta
    }
    # For transitively closed networks
    if(control$trans.close == TRUE){
      bestGraph <- as(transitive.closure(bestGraph),'matrix')
      diag(bestGraph) <- 0
      bestGraph <- bestGraph[Sgenes,Sgenes]
      bestPropMat <- getperturb.matrices(cmap,bestGraph,control)$PropMat
      bestScore <- DAGscore(bestPropMat,D1,D0,control,bestAlpha,bestBeta)
    }
    
    # Function to optimise the best noise para further using a solver
    optimNoise <- function(x){
        #-(DAGscore(bestPropMat,D1,D0,control,(exp(x[1]))/(1+exp(x[1])),(exp(x[2]))/(1+exp(x[2])))$mLL)
        -(DAGscore(bestPropMat,D1,D0,control,x[1],x[2])$mLL)
    }

    if(noiseEst == TRUE){
      #NoiseEst <- bobyqa(c(log(bestAlpha/(1-bestAlpha)),log(bestBeta/(1-bestBeta))), optimNoise, lower = c(1e-20, 1e-20), upper = c(0.5, 0.5))
      NoiseEst <- bobyqa(c(bestAlpha, bestBeta), optimNoise, lower = c(1e-15, 1e-15), upper = c(0.5, 0.5))
      
      return(list(graphs=L1,LLscores=L2,Temp=L3,AcceptRate=L4,AlphaVals=L5, 
                  BetaVals=L6,maxLLscore=bestScore,maxDAG=bestGraph,
                  typeI=NoiseEst$par[1],typeII=NoiseEst$par[2],minIter=minSteps,transformAR=transformAR))
    }else{
      return(list(graphs=L1,LLscores=L2,Temp=L3,AcceptRate=L4,AlphaVals=L5, 
                  BetaVals=L6,maxLLscore=bestScore,maxDAG=bestGraph,
                  typeI=bestAlpha,typeII=bestBeta,minIter=minSteps,transformAR=transformAR))
    }
    
    
}

