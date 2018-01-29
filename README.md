# pcNEM
Probabilistic combinatorial nested effects model (pc-NEM), is a variant of NEM, which probabilistically models the perturbations 
for network reconstruction from combinatorial gene knockdown data. The model and the inference algorithm are implemented as part of the R/Bioconductor package _nem_. This package will soon be merged with the existing _nem_ package.

#### Installation with devtools ####

```
install.packages("devtools") 
library(devtools) 
install_github("cbg-ethz/pcNEM")
```
#### Running pc-NEM #### 
```
data("BartonellaRNAi")

# Maximum likelihood estimate
set.seed(456)

# Setting all the control parameters
control <- set.default.parameters(unique(colnames(D)),type="mLL",pcombi = TRUE, trans.close=FALSE)
control$map <- as.matrix(KOmap)
pcnem_mle <- nem(D,inference="AdaSimAnneal",control=control,verbose=FALSE)
```
