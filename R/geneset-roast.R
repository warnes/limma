##  ROAST.R

setClass("Roast",
#  rotation gene set test
representation("list")
)

setMethod("show","Roast",
#  Di Wu, Gordon Smyth
#  14 May 2010.  Last modified 19 May 2010.
function(object) print(object$p.value)
)

roast <- function(y,...) UseMethod("roast")

roast.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,approx.zscore=TRUE,...)
# Rotation gene set testing for linear models
# Gordon Smyth and Di Wu
# Created 24 Apr 2008.  Last modified 22 March 2016.
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Check index
	if(is.list(index)) return(mroast(y=y,index=index,design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,array.weights=array.weights,weights=weights,block=block,correlation=correlation,var.prior=var.prior,df.prior=df.prior,trend.var=trend.var,nrot=nrot))

#	Extract components from y
	y <- getEAWP(y)
	ngenes <- nrow(y$exprs)
	n <- ncol(y$exprs)

#	Check index
	if(is.null(index)) index <- rep.int(TRUE,ngenes)

#	Check design
	if(is.null(design)) design <- y$design
	if(is.null(design))
		stop("design matrix not specified")
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

	p <- ncol(design)
	p0 <- p-1L
	d <- n-p

#	Check set.statistic
	set.statistic <- match.arg(set.statistic,c("mean","floormean","mean50","msq"))

#	Check var.prior and df.prior
	if(!is.null(var.prior) && var.prior<0) stop("var.prior must be non-negative")
	if(!is.null(df.prior) && df.prior<0) stop("df.prior must be non-negative")

#	Check array weights
	if(!is.null(array.weights)) {
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		if(any(array.weights <= 0)) stop("array.weights must be positive")
	}

#	Check observational weights
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		dimw <- dim(weights)
		if(dimw[1] != ngenes || dimw[2] != n) stop("weights must have same dimensions as y")
		if(any(weights <= 0)) stop("weights must be positive")
		if(!is.null(array.weights)) {
			weights <- .matvec(weights, array.weights)
			array.weights <- NULL
		}
	}

#	Reduce to numeric expression matrix
	y <- y$exprs

#	Divide out array weights
	if(!is.null(array.weights)) {
		sw <- sqrt(array.weights)
		design <- design*sw
		y <- .matvec(y,sw)
		array.weights <- NULL
	}

#	Divide out block correlation
	if(!is.null(block)) {
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)
		y <- t(backsolve(R, t(y), transpose = TRUE))
		dn <- dimnames(design)
		design <- backsolve(R, design, transpose = TRUE)
		dimnames(design) <- dn
 	}

#	Check contrast
	if(is.character(contrast)) {
		if(length(contrast)>1L) {
			warning("using only first entry for contrast")
			contrast <- contrast[1]
		}
		contrast <- which(contrast==colnames(design))
		if(length(contrast)==0L) stop("coef ",contrast," not found")
	}
	if(all(contrast == 0)) stop("contrast all zero")

#	Reform design matrix so that contrast is last coefficient
	if(length(contrast) == 1L) {
		contrast <- as.integer(contrast)
		if(contrast < p)
			X <- cbind(design[,-contrast,drop=FALSE],design[,contrast,drop=FALSE])
		else
			X <- design
	} else {
		if(length(contrast) != p) stop("length of contrast must match column dimension of design")
		X <- contrastAsCoef(design, contrast, first=FALSE)$design
	}

	qr <- qr(X)
	signc <- sign(qr$qr[p,p])

	if(is.null(var.prior) || is.null(df.prior)) {
#		Fit model to all genes
		if(is.null(weights)) {
			effects <- qr.qty(qr,t(y))
		} else {
			ws <- sqrt(weights)
			effects <- matrix(0,n,ngenes)
			signc <- rep.int(0,ngenes)
			for (g in 1:ngenes) {
				wX <- X*ws[g,]
				wy <- y[g,]*ws[g,]
				qrX <- qr(wX)
				signc[g] <- sign(qrX$qr[p,p])
				effects[,g] <- qr.qty(qrX,wy)
			}
			signc <- signc[index]
		}
#		Estimate global parameters s0 and d0
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(trend.var) covariate <- rowMeans(y,na.rm=TRUE) else covariate <- NULL
		sv <- squeezeVar(s2,df=d,covariate=covariate)
		d0 <- sv$df.prior
		s02 <- sv$var.prior
		if(trend.var) s02 <- s02[index]
		effects <- effects[,index,drop=FALSE]
		sd.post <- sqrt(sv$var.post[index])
	} else {
		d0 <- df.prior
		s02 <- var.prior
		if(length(s02)>1) {
			names(s02) <- rownames(y)
			s02 <- s02[index]
		}
		y <- y[index,,drop=FALSE]
		if(is.null(weights)) {
			effects <- qr.qty(qr,t(y))
		} else {
			ws <- sqrt(weights[index,,drop=FALSE])
			nset <- nrow(y)
			effects <- matrix(0,n,nset)
			signc <- rep.int(0,nset)
			for (g in 1:nset) {
				wX <- X*ws[g,]
				wy <- y[g,]*ws[g,]
				qrX <- qr(wX)
				signc[g] <- sign(qrX$qr[p,p])
				effects[,g] <- qr.qty(qrX,wy)
			}
		}
		s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)
		if(is.finite(d0))
			sd.post <- sqrt( (d0*s02+d*s2)/(d0+d) )
		else
			sd.post <- sqrt(s02)
	}

#	From here, all results are for set only
	nset <- ncol(effects)
	if(p0>0)
		Y <- effects[-(1:p0),,drop=FALSE]
	else
		Y <- effects
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post

	statobs <- p <- rep(0,4)
	names(statobs) <- names(p) <- c("down","up","upordown","mixed")
	statrot <- array(0,c(nrot,4),dimnames=list(NULL,names(p)))

#	Convert to z-scores
	modt <- zscoreT(modt,df=d0+d,approx=approx.zscore)

#	Active proportions
	if(!is.null(gene.weights)) {
		lgw <- length(gene.weights)
		if(lgw > nset && lgw==ngenes) {
			gene.weights <- gene.weights[index]
		} else {
			if(lgw != nset) stop("length of gene.weights disagrees with size of set")
		}
		s <- sign(gene.weights)
		ss <- sum(abs(s))
		r1 <- sum(s*modt > sqrt(2)) / ss
		r2 <- sum(s*modt < -sqrt(2)) / ss
	} else {
		r1 <- mean(modt > sqrt(2))
		r2 <- mean(modt < -sqrt(2))
	}

#	Random rotations
	R <- matrix(rnorm(nrot*(d+1)),nrot,d+1)
	R <- R/sqrt(rowSums(R^2))
	Br <- R %*% Y
	s2r <- (matrix(YY,nrot,nset,byrow=TRUE)-Br^2)/d
	if(is.finite(d0))
		sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
	else
		sdr.post <- sqrt(s02)
	modtr <- signc*Br/sdr.post
	modtr <- zscoreT(modtr,df=d0+d,approx=approx.zscore)

	switch(set.statistic,
	"mean" = { 
#		Observed statistics
		if(!is.null(gene.weights)) modt <- gene.weights*modt
		m <- mean(modt)
		statobs["down"] <- -m
		statobs["up"] <- m
		statobs["mixed"] <- mean(abs(modt))
#		Simulated statistics
		if(!is.null(gene.weights)) modtr <- t(gene.weights*t(modtr))
		m <- rowMeans(modtr)
		statrot[,"down"] <- -m
		statrot[,"up"] <- m
		statrot[,"mixed"] <- rowMeans(abs(modtr))
#		p-values
		p["down"] <- sum(statrot[,c("down","up")] > statobs["down"])
		p["up"] <- sum(statrot[,c("down","up")] > statobs["up"])
		p["upordown"] <- min(p[c("down","up")])
		p["mixed"] <- sum(statrot[,c("mixed")] > statobs["mixed"])
		p <- (p+1) / (c(2,2,1,1)*nrot + 1)
	},

	"floormean" = { 
#		Observed statistics
		chimed <- qchisq(0.5,df=1)
		amodt <- pmax(abs(modt),chimed)
		if(!is.null(gene.weights)) {
			amodt <- gene.weights*amodt
			modt <- gene.weights*modt
		}
		statobs["down"] <- mean(pmax(-modt,0))
		statobs["up"] <- mean(pmax(modt,0))
		statobs["upordown"] <- max(statobs[c("down","up")])
		statobs["mixed"] <- mean(amodt)
#		Simulated statistics
		amodtr <- pmax(abs(modtr),chimed)
		if(!is.null(gene.weights)) {
			amodtr <- t(gene.weights*t(amodtr))
			modtr <- t(gene.weights*t(modtr))
		}
		statrot[,"down"] <- rowMeans(pmax(-modtr,0))
		statrot[,"up"] <- rowMeans(pmax(modtr,0))
		i <- statrot[,"up"] > statrot[,"down"]
		statrot[i,"upordown"] <- statrot[i,"up"]
		statrot[!i,"upordown"] <- statrot[!i,"down"]
		statrot[,"mixed"] <- rowMeans(amodtr)
#		p-values
		p["down"] <- sum(statrot[,c("down","up")] > statobs["down"])
		p["up"] <- sum(statrot[,c("down","up")] > statobs["up"])
		p["upordown"] <- sum(statrot[,c("upordown")] > statobs["upordown"])
		p["mixed"] <- sum(statrot[,c("mixed")] > statobs["mixed"])
		p <- (p+1) / (c(2,2,1,1)*nrot + 1)
	},

	"mean50" = { 
		if(nset%%2L == 0L) {
			half1 <- nset %/% 2L
			half2 <- half1 + 1L
		} else {
			half1 <- half2 <- nset %/% 2L + 1L
		}
#		Observed statistics
		if(!is.null(gene.weights)) modt <- gene.weights*modt
		s <- sort(modt,partial=half2)
		statobs["down"] <- -mean(s[1:half1])
		statobs["up"] <- mean(s[half2:nset])
		statobs["upordown"] <- max(statobs[c("down","up")])
		s <- sort(abs(modt),partial=half2)
		statobs["mixed"] <- mean(s[half2:nset])
#		Simulated statistics
		if(!is.null(gene.weights)) modtr <- t(gene.weights*t(modtr))
		for (g in 1L:nrot) {
			s <- sort(modtr[g,],partial=half2)
			statrot[g,"down"] <- -mean(s[1:half1])
			statrot[g,"up"] <- mean(s[half2:nset])
			statrot[g,"upordown"] <- max(statrot[g,c("down","up")])
			s <- sort(abs(modtr[g,]),partial=half2)
			statrot[g,"mixed"] <- mean(s[half2:nset])
		}
#		p-values
		p["down"] <- sum(statrot[,c("down","up")] > statobs["down"])
		p["up"] <- sum(statrot[,c("down","up")] > statobs["up"])
		p["upordown"] <- sum(statrot[,c("upordown")] > statobs["upordown"])
		p["mixed"] <- sum(statrot[,c("mixed")] > statobs["mixed"])
		p <- (p+1) / (c(2,2,1,1)*nrot + 1)
	},

	"msq" = {
#		Observed statistics
		modt2 <- modt^2
		if(!is.null(gene.weights)) {
			modt2 <- abs(gene.weights)*modt2
			modt <- gene.weights*modt
		}
		statobs["down"] <- sum(modt2[modt < 0])/nset
		statobs["up"] <- sum(modt2[modt > 0])/nset
		statobs["upordown"] <- max(statobs[c("down","up")])
		statobs["mixed"] <- mean(modt2)
#		Simulated statistics   
		if(!is.null(gene.weights)) {
			gene.weights <- sqrt(abs(gene.weights))
			modtr <- t(gene.weights*t(modtr))
		}
		statrot[,"down"] <- rowMeans(pmax(-modtr,0)^2)
		statrot[,"up"] <- rowMeans(pmax(modtr,0)^2)
		i <- statrot[,"up"] > statrot[,"down"]
		statrot[i,"upordown"] <- statrot[i,"up"]
		statrot[!i,"upordown"] <- statrot[!i,"down"]
		statrot[,"mixed"] <- rowMeans(modtr^2)
#		p-values
		p["down"] <- sum(statrot[,c("down","up")] > statobs["down"])
		p["up"] <- sum(statrot[,c("down","up")] > statobs["up"])
		p["upordown"] <- sum(statrot[,c("upordown")] > statobs["upordown"])
		p["mixed"] <- sum(statrot[,c("mixed")] > statobs["mixed"])
		p <- (p+1) / (c(2,2,1,1)*nrot + 1)
	})

#	Output
	out <- data.frame(c(r2,r1,max(r1,r2),r1+r2),p)
	dimnames(out) <- list(c("Down","Up","UpOrDown","Mixed"),c("Active.Prop","P.Value"))
	new("Roast",list(p.value=out,var.prior=s02,df.prior=d0,ngenes.in.set=nset))
}

mroast <- function(y,...) UseMethod("mroast")

mroast.default <- function(y,index=NULL,design=NULL,contrast=ncol(design),set.statistic="mean",gene.weights=NULL,array.weights=NULL,weights=NULL,block=NULL,correlation,var.prior=NULL,df.prior=NULL,trend.var=FALSE,nrot=999,approx.zscore=TRUE,adjust.method="BH",midp=TRUE,sort="directional",...)
#  Rotation gene set testing with multiple sets
#  Gordon Smyth and Di Wu
#  Created 28 Jan 2010. Last revised 22 Dec 2015.
{
#	Extract components from y
	y <- getEAWP(y)
	ngenes <- nrow(y$exprs)
	n <- ncol(y$exprs)

#	Check index
	if(is.null(index)) index <- rep(TRUE,ngenes)
	if(!is.list(index)) index <- list(set = index)
	nsets <- length(index)
	if(nsets==0) stop("index is empty")
	if(is.null(names(index))) names(index) <- paste("set",1:nsets,sep="")

#	Check design matrix
	if(is.null(design)) design <- y$design
	if(is.null(design))
		stop("design matrix not specified")
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

#	Check gene.weights
	if(!is.null(gene.weights)) if(length(gene.weights) != ngenes) stop("gene.weights must have length equal to nrow(y)")

#	Check array.weights
	if(!is.null(array.weights)) {
		if(length(array.weights) != n) stop("array.weights wrong length")
		if(any(array.weights <= 0)) stop("array.weights must be positive")
	}

#	Check weights
	if(is.null(weights) && is.null(array.weights)) weights <- y$weights
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		dimw <- dim(weights)
		if(dimw[1] != ngenes || dimw[2] != n) stop("weights must have same dimensions as y")
		if(any(weights <= 0)) stop("weights must be positive")
		if(!is.null(array.weights)) {
			weights <- .matvec(weights,array.weights)
			array.weights <- NULL
		}
	}

#	Reduce to numeric expression matrix
	y <- y$exprs

#	Divide out array.weights
	if(!is.null(array.weights)) {
		sw <- sqrt(array.weights)
		design <- design*sw
		y <- .matvec(y,sw)
		array.weights <- NULL
	}

#	Divide out block correlation
	if(!is.null(block)) {
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)
		y <- t(backsolve(R, t(y), transpose = TRUE))
		design <- backsolve(R, design, transpose = TRUE)
		block <- NULL
 	}

#	Estimate var.prior and df.prior if not preset
	if(is.null(var.prior) || is.null(df.prior)) {
		fit <- lmFit(y,design=design,weights=weights)
		if(trend.var) {
			covariate <- fit$Amean
			if(is.null(covariate)) covariate <- rowMeans(y)
		} else {
			covariate=NULL
		}
		sv <- squeezeVar(fit$sigma^2,df=fit$df.residual,covariate=covariate)
		var.prior <- sv$var.prior
		df.prior <- sv$df.prior
	}

	pv <- adjpv <- active <- array(0,c(nsets,4),dimnames=list(names(index),c("Down","Up","UpOrDown","Mixed")))
	NGenes <- rep(0,nsets)
	if(nsets<1) return(pv)
	for(i in 1:nsets) {
		out <- roast(y=y,index=index[[i]],design=design,contrast=contrast,set.statistic=set.statistic,gene.weights=gene.weights,weights=weights,var.prior=var.prior,df.prior=df.prior,nrot=nrot,approx.zscore=approx.zscore,...)
		pv[i,] <- out$p.value$P.Value
		active[i,] <- out$p.value$Active.Prop
		NGenes[i] <- out$ngenes.in.set
	}

#	Use mid-p-values or ordinary p-values?
#	pv2 <- pv
#	if(midp) pv2 <- pv2-1/2/(nrot+1)
#	adjpv[,"Down"] <- p.adjust(pv2[,"Down"], method=adjust.method)
#	adjpv[,"Up"] <- p.adjust(pv2[,"Up"], method=adjust.method)
#	adjpv[,"Mixed"] <- p.adjust(pv2[,"Mixed"], method=adjust.method)
#	list(P.Value=pv, Adj.P.Value=adjpv, Active.Proportion=active)

#	New-style output
	Up <- pv[,"Up"] < pv[,"Down"]
	Direction <- rep.int("Down",nsets); Direction[Up] <- "Up"
	TwoSidedP2 <- pv[,"UpOrDown"]
	MixedP2 <- pv[,"Mixed"]
	if(midp) {
		TwoSidedP2 <- TwoSidedP2 - 1/2/(nrot+1)
		MixedP2 <- MixedP2 - 1/2/(nrot+1)
	}

	tab <- data.frame(
		NGenes=NGenes,
		PropDown=active[,"Down"],
		PropUp=active[,"Up"],
		Direction=Direction,
		PValue=pv[,"UpOrDown"],
		FDR=p.adjust(TwoSidedP2,method="BH"),
		PValue.Mixed=pv[,"Mixed"],
		FDR.Mixed=p.adjust(MixedP2,method="BH"),
		row.names=names(index),
		stringsAsFactors=FALSE
	)

	if(midp) {
		tab$FDR <- pmax(tab$FDR, pv[,"UpOrDown"])
		tab$FDR.Mixed <- pmax(tab$FDR.Mixed, pv[,"Mixed"])
	}

#	Sort by p-value
	sort <- match.arg(sort,c("directional","mixed","none"))
	if(sort=="none") return(tab)
	if(sort=="directional") {
		Prop <- pmax(tab$PropUp,tab$PropDown)
		o <- order(tab$PValue,-Prop,-tab$NGenes,tab$PValue.Mixed)
	} else {
		Prop <- tab$PropUp+tab$PropDown
		o <- order(tab$PValue.Mixed,-Prop,-tab$NGenes,tab$PValue)
	}
	tab[o,,drop=FALSE]
}

