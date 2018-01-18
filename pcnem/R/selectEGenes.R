filterEGenes = function(Porig, D, Padj=NULL, ntop=100, fpr=0.05, adjmethod="bonferroni", cutoff=0.05){
	if(is.null(Padj))
		Padj = apply(Porig, 2, p.adjust, method="fdr")	
	ntop = min(nrow(Porig),ntop)
	n = ncol(Porig)
	I1 = apply(Porig,2,function(x) order(x)[1:ntop])
	I1 = unique(as.vector(I1))  # Länge ist abhängig von P!	
	print(paste("Selecting top",ntop," genes from each list -->",length(I1),"genes total"))				
	disc = (Padj[I1,] <= fpr)*1		
	nsig = colSums(Padj[I1,] < fpr)
	N = nrow(disc)
	patterns = unique(disc)		
	patterns = patterns[-which(apply(patterns,1,function(r) all(r == 0))),]	
	if(nrow(patterns) < 1)
		stop("No patterns found!")	
	idx = apply(patterns,1, function(p){
		cl = which(apply(disc,1, function(r) all(r == p)))
	})
	nobserved = sapply(idx, length)
	patterns = patterns[nobserved > 0,]
	idx = idx[nobserved > 0]
	nobserved = nobserved[nobserved > 0]
	cat("Testing ", nrow(patterns), " patterns\n")
	p.values = sapply(1:nrow(patterns), function(j){
		p = patterns[j,]		
		pexpected = max(0,prod((fpr*p*nsig + (1-fpr)*(1-p)*(N-nsig))/N))
		p.value = binom.test(nobserved[j], N, pexpected, alternative="greater")$p.value
		cat("pattern ", p, ": (#observed = ", nobserved[j], ", #expected = ", floor(pexpected*N), ", raw p-value = ", p.value,")\n")		
		p.value
	})
	cat("\n")
	p.values = p.adjust(p.values,method=adjmethod)
	if(!any(p.values < cutoff))
		stop("No significant patterns found!\n")
	patterns = patterns[p.values < cutoff,]	
	idx = idx[p.values < cutoff]
	nobserved = nobserved[p.values < cutoff]
	p.values = p.values[p.values < cutoff]
	I = I1[unlist(idx)]	
	cat(length(p.values), " significant patterns -->", length(I), "E-genes in total\n")
	D = D[I,]	
	list(selected=I, dat=D, patterns=patterns, nobserved=nobserved, p.values=p.values)
}

getRelevantEGenes <- function(Phi, D, control, nEgenes=min(10*nrow(Phi), nrow(D))){			
# 	if(type %in% c("CONTmLLRatio", "CONTmLLMAP")){		
# 		sc = score(list(Phi), D, type=type, para=para, hyperpara=hyperpara, Pe=Pe, Pm=Pm, lambda=lambda, delta=delta, verbose=FALSE, graphClass="matrix")						
# 	}
# 	else{						
# 		L <- sapply(1:nrow(D), function(idx){									
# 			score(list(Phi), D[idx,,drop=FALSE], type=type, para=para, hyperpara=hyperpara, Pe=Pe[idx,,drop=FALSE], Pm=Pm, lambda=lambda, verbose=FALSE, graphClass="matrix")$mLL	
# 		})			
		
		L = score(list(Phi), D, control, verbose=FALSE, graphClass="matrix")$LLperGene[[1]]	
		L = L/sum(L)		
		if(control$type %in% c("CONTmLLDens", "CONTmLLBayes","CONTmLLMAP","CONTmLLRatio"))
			nEgenes = max(10,length(which(L>0)))		
		sel <- unique(order(L,decreasing=TRUE)[1:nEgenes])
		control$Pe = control$Pe[sel,]
		sc = score(list(Phi), D[sel,], control, verbose=FALSE, graphClass="matrix")		
# 	}					
	list(selected=sel, mLL=sc$mLL[[1]], pos=sc$pos[[1]], mappos=sc$mappos[[1]], LLperGene=sc$LLperGene[[1]])
}
