# pcNEM
Probabilistic combinatorial nested effects model (pc-NEM), is a variant of NEM, which probabilistically models the perturbations 
for network reconstruction from combinatorial gene knockdown data. The model and the inference algorithm are implemented as part of the R/Bioconductor package `nem`. 

### Installation with devtools ###

```
install.packages("devtools") 
library(devtools) 
install_github("cbg-ethz/pcNEM")
```
### Running pc-NEM ### 
Small toy example with 6 S-genes, 6 knockdown experiments, and 90 E-genes. First, we sample a random network and then generate a perturbation map with on-target probabilities set to 0.8 and off-target probabilities set to 0.2. Next, simulate data using the network and the perturbation map and setting FPR = 0.05 and FNR = 0.01. Then use pc-nem to learn the network and noise parameters. 
```
library(pcnem)
library(Rgraphviz)

set.seed(42)

# Sample a toy network with N S-genes and M E-genes
N      <- 6
M      <- 90
Sgenes <- paste0("S", 1:N)
Phi    <- pcnem:::sampleRndNetwork(Sgenes = Sgenes, trans.close = TRUE)

# Generate knockout map with one experiment for each Sgene
KOmap           <- matrix(0, nrow = N, ncol = N )
rownames(KOmap) <- paste0("K", 1:N)
colnames(KOmap) <- Sgenes
# On target probabilities set to 0.8 in this example
diag(KOmap)     <- 0.8
# Off-target probabilities set to 0.2 in this example
KOmap[sample(1:N^2, 4)] <- 0.2

# Generating data
alpha <- 0.05
beta  <- 0.01
D     <- pcnem:::sampleData.pcnem(Phi = Phi, M = M, map = KOmap, typeI.err = alpha, typeII.err = beta)$D

# Setting all the control parameters and runnin pc-nem
control            <- set.default.parameters(unique(colnames(D)),type="mLL",pcombi = TRUE, trans.close=FALSE)
control$map        <- as.matrix(KOmap)
control$iterations <- 10000
control$temper     <- TRUE
control$AcceptRate <- 0.1
pcnem_mle          <- nem(D, inference = "AdaSimAnneal", control = control, verbose = FALSE)

# True network likelihood
control$para <- c(alpha,beta)
true_mle     <- nem(D,inference="search",control=control,verbose=FALSE, models = list(Phi))

# Plotting true and inferred network
par(mfrow = c(1,2))
plot(as(Phi, "graphNEL"), main = "True network")
plot(pcnem_mle$graph, main = "MLE graph")

# Estimated noise
cat("Estimated type I error:", pcnem_mle$typeIEst, "\n")
cat("Estimated type II error:", pcnem_mle$typeIIEst)

```
The FNR is slightly underestimated due to finite sampling effect from small number of experiments and effects in this example. 

### pc-NEM parameters ###
You can tune several hyperparameters for your case using the `set.default.parameters()` function.

`pcombi` :  Logical parameter set to TRUE for pc-NEM. Default set to FALSE. <br/>
`temper` :  Binary parameter to choose between two variant implementations of adaptive simulated annealing. FALSE corresponds to adaptation of temperatures at varying intervals but cooled at a fixed rate. TRUE corresponds to the scheme described in the paper. Both are very similar in performance. Default is set to FALSE.<br/>
`iterations` :  Number of iterations for adaptive simulated annealing (ASA). Deafult is 2e4.<br/>
`stepsave` : The length of intervals to adapt tempertature and noise. Default is 1e2.<br/>
`revallowed` : Binary parameter for including reversal moves. Default allows reversals and is set to 1.<br/>
`AcceptRate` : The ideal acceptance rate for ASA.  <br/>
`Temp`  : Initial temperature. Default is 50. <br/>
`AdaptRate` : Rate of adaptation of temperature. Default is 0.3.<br/>
`noiseEst` : Binary parameter to include estimation of noise parameters. Default is TRUE. <br/>
`moveprobs` : Probability of moving between DAG space and noise space if  `noiseEst`  = TRUE. Default = `c(0.6,0.4)`.<br/>
`moveprobsNoise` : Probability of moving between alpa and beta space for noise estimation. Default = `c(0.5,0.5)`. <br/>
`sigma` : Initial covariance matrix for noise estimation.

## Contact ##
[Sumana Srivatsa](https://www.bsse.ethz.ch/cbg/group/people/person-detail.MjAyOTQw.TGlzdC81MTYsOTQ0ODM3Mzc2.html) <br/>
sumana.srivatsa (at) bsse.ethz.ch

