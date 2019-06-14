# pcNEM
Probabilistic combinatorial nested effects model (pc-NEM), is a variant of NEM, which probabilistically models the perturbations 
for network reconstruction from combinatorial gene knockdown data. The model and the inference algorithm are implemented as part of the R/Bioconductor package _nem_. 

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
Once you set default parameters using _set.default.parameters()_ function, you can tune several hyperparameters for your setting.

_pcombi_ :  Logical parameter set to TRUE for pc-NEM. Default set to FALSE. <br/>
_temper_ :  Two different implementations of adaptive simulated annealing. FALSE corresponds to adaptation of temperatures at varying intervals but cooled at a fixed rate. TRUE corresponds to the scheme described in paper. Both are very similar in performance. Default is set to FALSE.<br/>
_iterations_ :  Number of iterations for adaptive simulated annealing (ASA). Deafult is 2e4.<br/>
_stepsave_ : The length of interval to adapt tempertature and noise. Default is 1e2.<br/>
_revallowed_ : Binary parameter for including reversal moves. Default allows reversals and is set to 1.<br/>
_AcceptRate_ : The ideal acceptance rate for ASA.  <br/>
_Temp_  : Initial temperature. Default is 50. <br/>
_AdaptRate_ : Rate of adaptation of temperature. Default is 0.3.<br/>
_noiseEst_ : Binary parameter to include estimation of noise parameters. Default is TRUE. <br/>
_moveprobs_ : Probability of moving between DAG space and noise space if  _noiseEst_  = TRUE. Default = c(0.6,0.4).<br/>
_moveprobsNoise_ : Probability of moving between alpa and beta space for noise estimation. Default = c(0.5,0.5). <br/>
_sigma_ : Initial covariance matrix for noise estimation. Default = diag(x=1,2)/10)
