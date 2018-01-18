plot.nem <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, D=NULL, draw.lines=FALSE, palette="BlueRed", ...) {
	
	if (!(what%in%c("graph","mLL","pos"))) stop("\nnem> invalid plotting type: plot either 'graph', 'mLL', or 'pos'")
	
	if (what=="graph"){		
		gR = x$graph	
		if(numEdges(gR) == 0)	
			stop("Graph contains no edges - nothing to draw!")
		M = as(gR, "matrix")
		toremove = which((abs(M) <= thresh) & (abs(M) > 0), arr.ind=TRUE)
		if(nrow(toremove) > 0){
			for(i in 1:nrow(toremove))
				gR = removeEdge(from=nodes(gR)[toremove[i,1]], to=nodes(gR)[toremove[i,2]], gR)
		}
		if(SCC){			
			gR = SCCgraph(gR)$graph
			M = as(gR, "matrix")	
		}
		if(numEdges(gR) == 0)
			edgeattr=list()			
		else{					
			if(transitiveReduction)
				M = transitive.reduction(M)
			eDDn <- names(edgeDataDefaults(gR))
			if(!"weight" %in% eDDn)
				edgeDataDefaults(gR, "weight") <- 1
			if(!"label" %in% eDDn)
				edgeDataDefaults(gR, "label") <- 1
			if(!"arrowhead" %in% eDDn)
				edgeDataDefaults(gR, "arrowhead") = "normal"
			if(!"style" %in% eDDn)
				edgeDataDefaults(gR, "style") = "bold" 
			nodes <- colnames(M)
			nodenames = vector("character", length(M[abs(M) > 0]))
			probs = double(length(nodenames))	    
			arr = character(length(probs))	   
			penwidth = rep("bold",length(probs))
			k = 1
			for (i in 1:ncol(M)) {
				for (j in 1:nrow(M)) {
					if (M[i, j] != 0) {		
						if(class(x) != "dynoNEM")
							probs[k] = signif(ifelse(abs(M[i,j]) > 1, abs(M[i,j])-1, abs(M[i,j])), 2)
						else
							probs[k] = signif(M[i,j], 2)
						edgeData(gR, from = nodes[i], to = nodes[j], attr = "style") = "bold"
						edgeData(gR, from = nodes[i], to = nodes[j], attr = "label") = probs[k]
						edgeData(gR, from = nodes[i], to = nodes[j], attr = "weight") = M[i,j]
						if((M[i,j] > 0) & (M[i,j] <= 1) || class(x) == "dynoNEM"){
							edgeData(gR, from = nodes[i], to = nodes[j], attr = "arrowhead") = "normal"
							arr[k] = "normal"
						}
						else if(M[i,j] > 1 & class(x) != "dynoNEM"){
							edgeData(gR, from = nodes[i], to = nodes[j], attr = "arrowhead") = "vee"
							arr[k] = "vee"
						}
						else{
							edgeData(gR, from = nodes[i], to = nodes[j], attr = "arrowhead") = "tee"
							arr[k] = "tee"
						}
						nodenames[k] <- paste(nodes[i], "~", nodes[j],
						sep = "")
						k = k + 1
					}
					else{
						if(nodes[i] %in% unlist(inEdges(nodes[j], gR)))
							gR = removeEdge(from=nodes[i], to=nodes[j], gR)
					}
				}   
			}         
			names(arr) = nodenames 
			names(probs) = nodenames
			names(penwidth) = nodenames	
			fontcol = arr
			fontcol[arr == "tee"] = "blue"
			fontcol[arr == "normal"] = "black"
			fontcol[arr == "vee"] = "red"
			if(plot.probs)
				edgeattr = list(label = probs, arrowhead = arr, fontcolor = fontcol, color=fontcol, style=penwidth)
			else
				edgeattr = list(arrowhead = arr, fontcolor = fontcol, color=fontcol, style=penwidth)	
		}		
		edge.type = ifelse(plot.probs | class(x) == "dynoNEM", "distinct", "combined")
		el = buildEdgeList(gR, recipEdges=edge.type, edgeAttrs=edgeattr) 
		nodeattr=list(color=rep("white",length(nodes(gR))), penwidth=rep(0, length(nodes(gR))), fontsize=rep(14,length(nodes(gR))))			
		names(nodeattr$color)=nodes(gR)
		names(nodeattr$penwidth)=nodes(gR)
		names(nodeattr$fontsize)=nodes(gR)
		args = list(...)			
		if("nodeAttrs" %in% names(args))
			nodeattr = c(nodeattr, args[[match("nodeAttrs", names(args))]])
		if("edgeAttrs" %in% names(args))
			edgeattr = c(edgeattr, args[[match("edgeAttrs", names(args))]])			
		main=NULL			
		if("main" %in% names(args))
			main = args[["main"]]					
		G = agopen(gR,name="test",edges=el, edgeAttrs=edgeattr, nodeAttrs=nodeattr, recipEdges=edge.type)		
		
		if (PDF) pdf(file=filename)   
		par(cex.main=2, cex=1) 		
		if(is.null(D))
			plot(G, main=main)				
		else{
			zlim = NULL
			if("zlim" %in% names(args))
				zlim = args[["zlim"]]			
			plotnem(D, G, x, SCC=SCC, main=main, zlim=zlim, draw.lines=draw.lines, palette=palette)		
		}

		if (PDF) dev.off()
		save(gR, file=paste(unlist(strsplit(filename,".pdf")),".rda",sep=""))	
		toDotR(gR, paste(unlist(strsplit(filename,".pdf")),".dot",sep=""))
    }

	if(what=="mLL"){
		if(PDF) pdf(file=filename)
		par(cex=1.3)
		ss <- sort(unique(x$mLL),decreasing=TRUE)[1:min(30,length(x$mLL))]
		plot(x=1:length(ss), y=ss, pch=19, main="Score distribution",
			xlab=paste(length(ss),"top ranked models"),
			ylab="Marginal log-likelihood", 
			ylim=c(ss[length(ss)]-10,ss[1]+10)
			)
		points(1,max(unique(x$mLL)),pch=21,cex=1.7,lwd=2)
		if(PDF) dev.off()
	}
	
	if(what=="pos"){    		
		if(length(x$mLL) > 1 & class(x) != "dynoNEM"){
			winner <- which.max(x$mLL)
			pos <- x$pos[[winner]]
			effects <- rownames(x$pos[[winner]])
		}
		else{
			pos <- x$pos
			effects <- rownames(x$pos)
		}	
		pos[is.na(pos)] = 0
		if(PDF) pdf(file=filename)
		par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.3,cex.main=1.7)        	
		image(x=1:ncol(pos),
		y=1:nrow(pos),
		z = t(pos),
		main = "Posterior effect positions",
		xlab="Perturbations",
		xaxt="n",
		ylab="Effect reporters",
		yaxt="n",
		col=gray(seq(.95,0,length=10))
		)
		abline(v=(1:(ncol(pos)-1))+.5)
		axis(1,1:ncol(pos),colnames(pos))        
		axis(2,1:length(effects),effects)
		if(PDF) dev.off()
	}

}

plot.nem.consensus <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}
  
plot.nem.bootstrap <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.nem.jackknife <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.nem.greedy <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.ModuleNetwork <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.score <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.pairwise <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.triples <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.nem.greedyMAP <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.nem.BN <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.mc.eminem = function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs, SCC, ...)
}

plot.dynoNEM = function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs=TRUE, SCC=FALSE, ...)
}

plot.score.list = function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...) {
	x$mappos = x$mappos[[1]]
	x$pos = x$pos[[1]]
	plot.nem(x, what, remove.singletons, PDF, filename, thresh, transitiveReduction, plot.probs=TRUE, SCC=FALSE, ...)
}
