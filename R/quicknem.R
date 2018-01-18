quicknem <- function(D,type="CONTmLLDens",inference="nem.greedy",controls.name=NULL,contrasts=NULL,normalize=FALSE,cutoff=0.05,DIR="bum",plot=TRUE,bootstrap=0, ...) {	
	# preprocess data
	if(class(D)=="ExpressionSet") {
		print("Extracting data from expression set...")
		D <- exprs(D)
	}
	dat <- as.matrix(D)
	## do data normalization
	if(normalize) {
		print("Quantile normalization...")
		order <- colnames(dat)
		nas <-apply(dat,2,function(x) all(is.na(x)))
		dat <- cbind(normalizeQuantiles(dat[,!nas]),dat[,nas])[,order]
	}
	# run DEPPNs
	if(type=="depn") {
		sgenes <- colnames(dat)
		control <- set.default.parameters(sgenes, type=type,...)
		result <- nem(dat,control=control,...)
	} else {
		if(is.null(contrasts)&&is.null(controls.name)) {
			stop("Please provide a description of the controls or the desired contrasts.")
		}
		if(length(as.vector(controls.name))>1) {
			stop("Please give a single name for the control experiment or specify contrasts explicitely.")
		}
		# check whether the control arguments for the discrete NEM-version are present
		if(type %in% c("mLL","FULLmLL")) {
			args = list(...)
			args.na = names(args)
			if(!any(c("neg.control","pos.control") %in% args.na)) {
				stop("Please provide at least a negative or a positive control when using type=\"mLL\" or \"FULLmLL\"")
			}
		}
		
		## limma
		print("Peforming limma analysis...")
		targets <- colnames(dat)
		samples <- rownames(dat)
		targets.u <- unique(targets)
		design <- model.matrix(~ -1 + factor(targets,levels=targets.u))
		colnames(design) <- targets.u
		
		# all experiments versus the control columns
		# or defined contrasts
		if(is.null(contrasts)) {
			controls.ind <- grep(controls.name,colnames(dat))
			controls.na <- unique(colnames(dat)[controls.ind])
			sgenes <- targets.u[-match(controls.na,targets.u)]
			contrasts <- paste(sgenes,controls.na,sep="-")
		} else {
			sgenes <- contrasts
		}
		
		# make the contrasts matrix
		contr <- sapply(contrasts, makeContrasts, levels=design)
		rownames(contr) <- targets.u
		colnames(contr) <- contrasts
		
		# fit the model
		fit <- lmFit(dat,design=design)
		fit2 <- contrasts.fit(fit, contr)
		fitE <- eBayes(fit2)
		
		# extract results
		print(paste("Extracting the differential genes for p-value cutoff ", cutoff,sep=" "))		
		pvals <- apply(fitE$p.value, 2, p.adjust, method="BH")
		diff <- apply(pvals, 1, function(x,cutoff) any(x<cutoff), cutoff=cutoff)
		diffgenes <- names(diff[which(diff==TRUE)])
		pvals <- pvals[diffgenes,]
		pvals.noadj <- fitE$p.value[diffgenes,]
		logodds <- fitE$lods[diffgenes,]
		
		# discretize data
		if(type=="mLL" || type=="FULLmLL") {
			print("Discretizing data...")
			dat.discr <- dat[diffgenes,]
			D <- nem.discretize(dat.discr,args[["neg.control"]],args[["pos.control"]])$dat
			sgenes <- colnames(D)
		}
		# estimate p-value densities	
		if(type=="CONTmLLDens" || type=="CONTmLLBayes") {
			if(!file_test("-d",DIR)) {
				dir.create(DIR)				
			}
			print("Estimating BUM density matrix...")
			D <- getDensityMatrix(pvals.noadj,dirname=DIR)
		}
		# take log-odds ratios from limma
		if(type=="CONTmLLMAP" || type=="CONTmLLRatio"){
			D <- logodds
		}
		# take effect probabilities from limma
		if(type=="CONTmLL") {
			D <- pvals
		}
		# starting NEM calculation
		print("Starting NEM calculation...")
		control <- set.default.parameters(sgenes, type=type, ...)
		if(bootstrap!=0) {
			result <- nem.bootstrap(D,inference=inference,control=control,nboot=bootstrap)
		} else {
			result <- nem(D,inference=inference,control=control)			
		}

	}
	if(plot) {
		plot(result$graph)
	}
	return(result)
}

