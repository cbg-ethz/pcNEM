nem.consensus <- function(D,thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE){
	if(inference == "dynoNEM")
		stop("nem.consensus is not applicable for dynoNEMs")
	overlapBoot = as(nem.bootstrap(D,thresh, nboot,inference,models,control, verbose)$graph,"matrix")
	overlapJack = as(nem.jackknife(D,thresh, inference,models,control, verbose)$graph,"matrix")
	consens = ((overlapBoot > thresh) & (overlapJack > thresh))*1
	res = nem(D,models=list(consens),inference="search",control, verbose=verbose)
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	res$para = res$para[[1]]
	res$control= control	
	class(res) <- "nem.consensus"
	res
}
