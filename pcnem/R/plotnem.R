plotnem = function(D, G, x, SCC, main=NULL, zlim=NULL, draw.lines=FALSE, palette="BlueRed"){ 
  if(!"Rgraphviz" %in% loadedNamespaces())
    stop("Rgraphviz is required for plotting!")
	if(is(D, "list"))
		stop("Plotting of effects matrix is only implemented for one dataset, not for several.")
	if(length(x$mLL) > 1 & class(x) != "dynoNEM"){
		winner <- which.max(x$mLL)
		mappos <- x$mappos[[winner]]		
	}
	else{
		mappos <- x$mappos	
		if(length(mappos) == 1)
			mappos = mappos[[1]]	
	}						
	nf = split.screen(rbind(c(0,1,0.4,1),c(0,0.8,0,0.4), c(0.9,1,0,0.4)))
	erase.screen(2)
	screen(2)	
	mynodes = AgNode(G)
	nodenames = sapply(mynodes, name)		
	xy = getNodeXY(G)
	left = xy$x <= max(xy$x)*0.5
	right = xy$x > max(xy$x)*0.5	
	nodenames = c(nodenames[left][order(xy$x[left] + xy$y[left])], nodenames[right][order(xy$x[right] + xy$y[right],decreasing=TRUE)])		
	if(!is.null(zlim))
		ord = plotEffects(D, x, legend=FALSE, order=nodenames, orderSCC=SCC, zlim=zlim, palette=palette)	
	else
		ord = plotEffects(D, x, legend=FALSE, order=nodenames, orderSCC=SCC, palette=palette)	
	erase.screen(3)		
	screen(3)
	nrcolors = 200; half = 1+nrcolors/2 # nrcolors must be an even number
	if(palette == "BlueRed")
	 	colpal = c(brewer.pal(9,"Blues")[9:1],brewer.pal(9,"OrRd")[1:9])
	else if(palette == "Grey")
		colpal = brewer.pal(9, "Greys")
	else
		stop("Unknown palette!")
	allcolors = colorRampPalette(colpal)(nrcolors)
	if(is.null(zlim)){
		r = c(quantile(D[D<0],0.95), quantile(D[D>0],0.95))
		rangeall = c(-max(abs(r),na.rm=TRUE), max(abs(r),na.rm=TRUE))	
	}
	else
		rangeall = zlim
	color.legend(0.4,0.1,1,1,signif(seq(rangeall[1],rangeall[2],length.out=5),digits=1),rect.col=allcolors,gradient="y", cex=0.75)

	erase.screen(1)
	screen(1)
	par(mar=c(0,0,0,0), cex=1, cex.main=2)
	plot(G, main=main)
	if(draw.lines){
		c = 0.7912	# Problem: Woher bekomme ich den rechten Rand des unteren Bildes in Koordinaten des oberen?
		ma = c*par()$usr[2]			
		mi = getX(botLeft(boundBox(G)))		
		may = getY(upRight(boundBox(G)))			
		miy = min(xy$y)	 # weiteres Problem bei größeren Datensätzen: unteres Bild wird über oberes geschoben. Wo ist das tatsächliche Minimum des oberen?
		if(miy - max(getNodeHeight(G)) > 0)	
			miy = 0
		xrescale = function(x,a2=mi,b2=ma){		
			alpha = (b2-a2)/length(ord)
			beta = mi
			alpha*x + beta
		}
		for(i in 1:length(mynodes)){
			xy = getNodeXY(mynodes[[i]])				
			nodenames = strsplit(name(mynodes[[i]]),":")[[1]]				
			to = xrescale(match(unique(unlist(mappos[nodenames])),ord)-0.5)		
			to = sort(to)				
			for(j in 1:length(to)){			
				segments(xy$x, max(miy, xy$y-7), to[j], miy, col="grey", lty=1, lwd=1, type="b", pch=19, cex=0.3)			
			}		
		}		
		screen(1)
		plot(G, main=main)
	}	
	close.screen(c(1,2,3),all.screens=TRUE)					
}
