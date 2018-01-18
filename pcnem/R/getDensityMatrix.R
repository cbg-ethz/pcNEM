fitBUM <- function(cname,Porig,dirname,startab,startlam,tol){			
	cat("Fitting BUM model for S-gene ",cname,"\n")
	x = Porig[,cname]
	res = bum.EM(x,startab,startlam,tol)
	a = res$a
	lam = res$lambda
	cat("--> final values for alpha, beta: ",a,"\n")
	cat("--> final values for mixing coefficients: ",lam,"\n")		
	if(!is.null(dirname))
		pdf(file=file.path(dirname,paste(cname,"_histogram.pdf",collapse="")))
	else
		x11()
	bum.histogram(x,a,lam,breaks=100)	
	if(!is.null(dirname))	
		pdf(file=file.path(dirname,paste(cname,"_QQplot.pdf",collapse="")))
	else
		x11()
	qqbum(x,a,lam)
	if(is.null(dirname))
		readline("Press any key to continue...")	
	dev.off()
	dev.off()
	if(!is.null(dirname)){
		params = list(a=a,lam=lam)
		save(params, file=file.path(dirname,paste(cname,"_parameter.rda",collapse="")))
	}
	dens = bum.dalt(x,a,lam)	
}

getDensityMatrix = function(Porig, dirname=NULL, startab=c(0.3,10), startlam=c(0.6,0.1,0.3), tol=1e-4){
	D = sapply(colnames(Porig), fitBUM, Porig=Porig, dirname=dirname, startab=startab, startlam=startlam, tol=tol)	
	log(D)
}
