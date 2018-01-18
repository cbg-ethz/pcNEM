###################################
local.model.prior <- function(size,n,bias){
  
  nrE <- n*(n - 1)/2
  p <- size/nrE
  model.prior <- c(1-p,p/(bias+2),p/(bias+2),bias*p/(bias+2))
  names(model.prior) <- c("..","->","<-","<->")
  if (!(all(model.prior>0) & sum(model.prior)==1)) stop("That's not a distribution - maybe 'size' is too big!")
  return(model.prior)
}
