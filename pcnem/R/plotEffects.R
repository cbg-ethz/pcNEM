plotEffects <- function(D,nem,border=TRUE,legend=TRUE,order=NULL,orderSCC=TRUE,palette="BlueRed",...){

if(!(class(D) %in% c("matrix","data.frame"))) stop("First argument has to be the data matrix and second the nem object!")
# int        <- unique(colnames(x))
sccg       <- SCCgraph(nem$graph,name=TRUE)
topo.order <- tsort(sccg$graph)
if(is.null(order))
	myorder = topo.order
else{		
	if(!orderSCC)
		myorder = rev(unique(sccg$which.scc[match(order, names(sccg$which.scc))]))	
	else
		myorder = rev(names(sccg$scc[order]))
}

#----------------------------
# reorder COLUMNS            
#----------------------------

# order
# ord        <- unlist(sccg$scc[topo.order])
# col.order  <- unlist(lapply(ord,function(y) which(colnames(x) == y)))
# D          <- x[,col.order]
# Dcn        <- colnames(D)

#----------------------------
# reorder ROWS               
#----------------------------

# estimate effect positions
# old: take them directly from hierarchy:
#      if (class(nem)=="pairwise") mappos <- nem$mappos
#      if (class(nem)=="score")    mappos <- nem$mappos[[which.max(nem$mLL)]]
# new: estimate them from scc graph
# colnames(D) <- sccg$which.scc[Dcn]
# M <- as(sccg$graph,"matrix") + diag(length(sccg$scc)) 
# mappos <- score(D,models=M,type=nem$type,para=nem$para,hyperpara=nem$hyperpara,verbose=FALSE)$mappos[[1]]
# colnames(D) <- Dcn
if(class(nem) == "score"){
	mappos = nem$mappos[[which.max(nem$mLL)]]	
}
else{
	mappos = nem$mappos
	if(length(mappos) == 1)
		mappos = mappos[[1]]	
}
selected = nem$selected
#if(nem$control$type != "CONTmLLMAP"){
#	if(!is.null(rownames(D)))
#		null.genes = setdiff(rownames(D), selected)	
#	else
#		null.genes = setdiff(1:nrow(D), selected)
#	mappos[["null"]] = null.genes	
#}
#else{
	null.genes =  unique(unlist(mappos["null"], use.names=FALSE))
#}
if(length(null.genes) > 0){
	sccg$scc[["null"]] = "null"
	myorder = c("null", myorder)
}
# order 
v <- list()
nr <- list()
for (i in myorder){	
	w = unique(unlist(mappos[sccg$scc[[i]]], use.names=FALSE))	
	
	if (length(w)==0){ 
		v[[i]]  <- NA
		nr[[i]] <- 1 
	}
	if (length(w)==1){
		v[[i]] <- w
		nr[[i]] <- 1
	}  
	if (length(w) >1){
		d       <- dist(D[w,,drop=FALSE],method="manhattan")
		v[[i]]  <- w[hclust(d)$order]
		nr[[i]] <- length(w) 
	}    
}
v <- unique(rev(unlist(v)))
nr <- rev(unlist(nr))

D2 <- matrix(0,nrow=length(v),ncol=ncol(D))
D2[which(!is.na(v)),] <- D[v[which(!is.na(v))],]
dimnames(D2) <- list(rownames(D)[v],colnames(D))
colorder = unlist(sccg$scc[topo.order], use.names=FALSE)
D2 = D2[,colorder]
#----------------------------
# PLOT                       
#----------------------------
cs <- cumsum(nr)[-length(nr)]

nrcolors = 200; half = 1+nrcolors/2 # nrcolors must be an even number
if(palette == "BlueRed")
	colpal = c(brewer.pal(9,"Blues")[9:1],brewer.pal(9,"OrRd")[1:9])
else if(palette == "Grey")
	colpal = brewer.pal(9, "Greys")[9:1]
else
	stop("Unknown palette!")
allcolors = colorRampPalette(colpal)(nrcolors)

maxwrittenrows = 60    	
if(legend){
	nf = split.screen(rbind(c(0,0.75,0,1),c(0.75,1,0,1)))
	erase.screen(1)
}
if(nrow(D2) < maxwrittenrows)
	par(las=2,mgp=c(5.5,1,0),mar=c(5,1.235,0,0),cex.lab=1.7,cex.main=2,lwd=2, oma=c(8,0,0,0))
else
	par(las=2,mgp=c(5.5,1,0),mar=c(0,1.235,0,0),cex.lab=1.7,cex.main=2,lwd=2, oma=c(8,0,0,0))

args = list(...)
if("zlim" %in% names(args))
	rangeall = args[["zlim"]]
else{
	r = c(quantile(D2[D2<0],0.95), quantile(D2[D2>0],0.95))
	rangeall = c(-max(abs(r),na.rm=TRUE), max(abs(r),na.rm=TRUE))
}
D2[D2<rangeall[1]] = rangeall[1]
D2[D2>rangeall[2]] = rangeall[2]
a = (length(allcolors)-1)/(rangeall[2] - rangeall[1])
b = 1 - a*rangeall[1]
colors = allcolors[round(a*min(D2) + b):round(a*max(D2) + b)]

image(x=1:nrow(D2),y=1:ncol(D2),z=D2,xaxt="n",yaxt="n",xlab="",ylab="",col=colors, ...)
if(nrow(D2) < maxwrittenrows)
	axis(side=1,at=1:nrow(D2),labels=rownames(D2),las=2,cex.axis=0.75, tick=!is.na(v))
axis(side=4,at=1:ncol(D2),labels=colnames(D2),tick=FALSE,las=1, cex.axis=0.75) 
axis(side=1, at=c(0,cs) + nr/2 ,labels=rev(myorder),tick=FALSE, cex.axis=0.6, las=2, outer=TRUE)
box()

# image(x   = 1:ncol(D2),
#       y   = 1:nrow(D2),
#       z   = t(D2),
#       xlab= "",
#       ylab= "",
#       xaxt= "n",
#       yaxt= "n",
#       col = gray(seq(.95,0,length=10)),
#       ...
#       )
# axis(1,1:ncol(D2),colnames(D2))
# axis(4,1:nrow(D2),rownames(D2),tick=!is.na(v))
# axis(2, c(0,cs) + nr/2 , rev(topo.order),tick=FALSE)
if (border) abline(v=cs+.5,col="black",lwd=2)

if(legend){
	erase.screen(2)
	screen(2)	
	par(mar=c(0,0,0.1,1))
	if(palette == "BlueRed")
		color.legend(0.5,0.1,1,1,signif(seq(rangeall[1],rangeall[2],length.out=10),digits=1),rect.col=allcolors,gradient="y", cex=0.75)
	else if(palette == "Grey")
		color.legend(0.5,0.1,1,1,signif(seq(rangeall[1],rangeall[2],length.out=10),digits=1),rect.col=allcolors,gradient="y", cex=0.75)
	close.screen(c(1,2),all.screens=TRUE)
}
return(v)
}
