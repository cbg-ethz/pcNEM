# pcNEM
Probabilistic combinatorial nested effects model (pc-NEM), is a variant of NEM, which probabilistically models the perturbations 
for network reconstruction from combinatorial gene knockdown data. The model and the inference algorithm are implemented as part of the R/Bioconductor package `nem`. 

#### Installation with devtools ####

```
install.packages("devtools") 
library(devtools) 
install_github("cbg-ethz/pcNEM")
```
#### Running pc-NEM #### 
```
data("BartonellaRNAi2017")

set.seed(456)

# Setting all the control parameters
control <- set.default.parameters(unique(colnames(D)),type="mLL",pcombi = TRUE, trans.close=FALSE)
control$map <- as.matrix(KOmap)
pcnem_mle <- nem(D,inference="AdaSimAnneal",control=control,verbose=FALSE)
```
#### pc-NEM parameters #### 
You can tune several hyperparameters for your setting using the `set.default.parameters()` function.

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
