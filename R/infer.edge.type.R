infer.edge.type = function(x, logFC, alpha=0.05, adj.method="BY", method=c("direct", "indirect")){	
	method = match.arg(method, c("direct", "indirect"))
	if(length(x$mLL) > 1 & class(x) != "dynoNEM"){
        	winner <- which.max(x$mLL)
		pos <- x$mappos[[winner]]
	}
	else{
		pos <- x$mappos
	}		
	g = x$graph	
	myedges = edges(x$graph)
	myedges = myedges[sapply(myedges,length) > 0]
	if(length(edgeDataDefaults(g)) == 0){
		edgeDataDefaults(g, "label") <- 1
		edgeDataDefaults(g, "weight") <- 1
	}
	Phi = as(g, "matrix")
	if(method == "direct"){
		for(i in rownames(Phi)){
			for(j in colnames(Phi)){
				if(i != j && Phi[i,j] != 0){
					if(j %in% rownames(logFC)){
						tmp = sign(logFC[i,j])
						if(tmp < 0) # A->B and B goes down => activation
							tmp = Phi[i, j] + 1
						else # A-|B and B goes up => inhibition
							tmp = -Phi[i, j]
						edgeData(g, from = i, to = j, attr = "weight") = tmp
						edgeData(g, from = i, to = j, attr = "label") = Phi[i, j]
					}
				}
			}
		}
	}
	else{
		p.values = c()
		k = 1
		for(eout in names(myedges)){				
			for(ein in myedges[[eout]]){						
				if(ein %in% rownames(logFC))
					eff = unique(c(pos[[ein]], ein))
				else
					eff = pos[[ein]]
				intype =  table(sign(logFC[eff,eout]))			
				if(length(intype) > 1)
					pval = binom.test(intype, length(eff))$p.value	
				else
					pval = 0								
				p.values = c(p.values, pval)	
				names(p.values)[k] = paste(eout, ein,sep="~")
				k = k + 1			
			}
		}			
		p.values = p.adjust(p.values, method=adj.method)	
		
		for(eout in names(myedges)){				
			for(ein in myedges[[eout]]){				
				if(ein %in% rownames(logFC))
					eff = unique(c(pos[[ein]], ein))
				else
					eff = pos[[ein]]
				intype =  table(sign(logFC[eff,eout]))				
				if(length(intype) > 1)
					type = -(intype[1] > intype[2])*1				
				else
					type = names(intype)
				pval = p.values[paste(eout, ein,sep="~")]
				
				if(pval < alpha){
					if(type > 0) # more are upregulated than downregulated					
						tmp = -Phi[eout,ein] # inhibition
					else
						tmp = Phi[eout,ein] + 1# activation
					edgeData(g, from = eout, to = ein, attr = "weight") = tmp
					edgeData(g, from = eout, to = ein, attr = "label") = Phi[eout, ein]
				} # otherwise not clear
			}
		}
	}	
	
	x$graph = g
	x
}

