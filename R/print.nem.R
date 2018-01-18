print.nem <- function(x, ...) {

  # general
  cat("Object of class ",class(x),"\n")
  cat("\n")
  
  # slots
  cat("$graph:  phenotypic hierarchy (graphNEL object) with",ncol(x$graph),"genes\n")  
  cat("Inference scheme: ",x$control$type,"\n")
  if(length(x$mLL) == 1)
  	cat("log posterior (marginal) likelihood $mLL:", x$mLL, "\n")  
  if(x$control$type == "mLL")
  	cat("Error probabilities alpha and beta:", x$control$para,"\n")
  if(x$control$type == "FULLmLL")
  	cat("Hyperparameters for error probability distributions:", x$control$hyperpara, "\n")
  cat("network structure regularization parameter $lambda (default: 0):",x$control$lambda ,"\n")  
  cat("Prior weight $delta for assigning E-genes to virtual S-gene 'null' (default: 1):",x$control$delta ,"\n")
  if(!is(x, "score.list")){
	  cat(length(x$selected), " selected E-genes:\n")
	  for(i in 1:length(x$mappos)){
		cat("-->", names(x$mappos)[i], ":", length(x$mappos[[i]]), " attached E-genes\n")
	  }
	  cat("\nNOTE: One E-gene can be attached to multiple S-genes\n")
	  cat("\n")
  }     
}
  
print.nem.greedy = function(x, ...){
	print.nem(x, ...)
}

print.ModuleNetwork= function(x, ...){
	print.nem(x, ...)
}

print.pairwise = function(x, ...){
	print.nem(x, ...)
	cat("$scores: posterior distributions of local models\n")
	cat("\n")
	
	# summary
	cat("Summary of MAP estimates:\n") 
	tmp         <- table(apply(x$scores[,1:4],1,which.max))
	summ        <- c(sum(tmp),tmp[1],tmp[2]+tmp[3],tmp[4])
	names(summ) <- c("all","..","->","<->")
	print(summ)
}

print.triples = function(x, ...){
	print.nem(x, ...)
}

print.nem.greedyMAP = function(x, ...){
	print.nem(x, ...)
}

print.nem.jackknife = function(x, ...){
	print.nem(x, ...)
}

print.nem.bootstrap = function(x, ...){
	print.nem(x, ...)
}

print.nem.consensus = function(x, ...){
	print.nem(x, ...)
}

print.nem.BN = function(x, ...){
	print.nem(x, ...)
}

print.mc.eminem = function(x, ...){
	print.nem(x, ...)
}

print.dynoNEM = function(x, ...){
	print.nem(x, ...)
}

print.score <- function(x, ...) {	
	cat("scores for ",length(x$mLL)," models\n")	
	best = which.max(x$mLL)
	cat("--> best model is number ", best,"\nInformation on this model:\n")			
	x$mLL = x$mLL[best]
	x$mappos = x$mappos[[best]]	
	print.nem(x, ...)
	#cat("\n")  
	#cat("plot this object to see the graph\n")
}

print.score.list <- function(x, ...) {	
	cat("scores for ",length(x$mLL)," models\n")	
	best = which.max(x$mLL)
	cat("--> best model is number ", best,"\nInformation on this model:\n")			
	x$mLL = x$mLL[best]
	print.nem(x, ...)	
	#cat("\n")  
	#cat("plot this object to see the graph\n")
}
