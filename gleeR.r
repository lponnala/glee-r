
#' check the spectral count data
#'
#' dimensions must match with specification, and all values must be positive
#'
#' @export
data_ok = function(D, nA, nB) {
	# ensure that the input file has as many columns as claimed, and at least one row
	if (dim(D)[2] != (1 + nA + nB )) { return(FALSE) }
	if (dim(D)[1]<2) { return(FALSE) }
	
	# ensure that the spectral-count data contains only finite positive values
	dat = tryCatch({ matrix(as.numeric(unlist(D[,-1])), nrow=dim(D)[1]) }, error = function(e) { NA }, warning = function(w) { NA })
	if (!all(is.finite(dat)) || any(dat<0)) { return(FALSE) }
	
	return(TRUE)
}

#' fit the specified model
#'
#' in each condition separately
#'
#' @export
fit_model = function(A,B,fit_type) {
	conditions = c("A","B")
	
	xbar = cbind(rowMeans(A), rowMeans(B))
	stdev = cbind(apply(A, FUN=sd, MARGIN=1), apply(B, FUN=sd, MARGIN=1))
	
	m = vector("list", length(conditions))
	names(m) = conditions
	for (j in 1:length(conditions)) {
		ind_valid = is.finite(xbar[,j]) & xbar[,j]>0 & is.finite(stdev[,j]) & stdev[,j]>0
		xbar_log = log(xbar[ind_valid,j])
		stdev_log = log(stdev[ind_valid,j])
		
		if (fit_type=="linear") {
			model = tryCatch({ lm(stdev_log ~ xbar_log) }, error = function(e) { NA }, warning = function(w) { NA })
		} else if (fit_type=="cubic") {
			model = tryCatch({ lm(stdev_log ~ poly(xbar_log,3)) }, error = function(e) { NA }, warning = function(w) { NA })
		}
		if (class(model) != "lm") { stop("error fitting model: condition ", conditions[j]) }
		
		adjRsq = summary(model)$adj.r.squared
		coeffs = model$coefficients
		
		# pick 20 values to show model fitted values
		xbar_showfit = seq(min(xbar_log), max(xbar_log), length.out=20)
		stdev_showfit = predict(model, data.frame(xbar_log = xbar_showfit))
		
		# remove outliers, i.e. keep points that are within 5 "sigma" of stdev
		ind_show = (stdev_log>(mean(stdev_log)-5*sd(stdev_log))) & (stdev_log<(mean(stdev_log)+5*sd(stdev_log)))
		xbar_showpts = xbar_log[ind_show]
		stdev_showpts = stdev_log[ind_show]
		
		# collect into a list to be returned
		m[[j]] = list(model = model, adjRsq = adjRsq, coeffs = coeffs, xbar_showfit = xbar_showfit, stdev_showfit = stdev_showfit, xbar_showpts = xbar_showpts, stdev_showpts = stdev_showpts, xbar = xbar[,j], stdev = stdev[,j])
	}
	return(m)
}

#' plot model-fit
#'
#' mean-vs-stdev and signal-level
#'
#' @export
model_fit_plots = function(m, outfile) {
	png(filename=outfile, width=800, height=800, units="px")
	par(mfrow=c(2,2), mar=c(5,4,4,2)+0.1)
	for (j in 1:length(m)) {
		plot(m[[j]]$xbar_showpts, m[[j]]$stdev_showpts, type="p", pch=20, col="gray", xlab="log(mean)", ylab="log(stdev)", main=paste("Condition ",names(m)[j],sep=""))
		lines(m[[j]]$xbar_showfit, m[[j]]$stdev_showfit, lty=1, col="blue", lwd=2)
		legend("topleft", paste("adj.Rsq = ",round(m[[j]]$adjRsq,3),sep=""), bty="n")
		plot(sort(m[[j]]$xbar_showpts), type="p", pch=20, col="gray", xlab="protein", ylab="signal level (mean)", main=paste("Condition ",names(m)[j],sep=""))
	}
	dev.off()
}

#' calculate standard deviation
#'
#' using the log-values of xbar
#'
#' @export
calcs = function(model, x) {
	exp(predict(model, data.frame(xbar_log=log(x))))
}

#' calculate signal-to-noise and p-values
#'
#' using resampled distribution
#'
#' @importFrom dplyr %>%
#'
#' @export
calc_stn_pval = function(A, B, m, num_iter) {
	best_idx = as.numeric(which.max(sapply(m, FUN = function(x) { x$adjRsq })))
	best_model = m[[best_idx]]$model
	
	xbar = c()
	for (j in 1:length(m)) {
		xbar_j = m[[j]]$xbar
		# replace xbar zeroes with min-positive xbar from the same condition
		xbar_j[xbar_j==0] = min(xbar_j[xbar_j>0])
		xbar = cbind(xbar, xbar_j)
	}
	min_xbar_value = apply(xbar, MARGIN=2, FUN=min)

	# calculate model-based STN
	if (any(xbar[,1] <= 0) || any(xbar[,2] <= 0)) { stop("xbar has non-positive values") }
	model_stn = (xbar[,2] - xbar[,1]) / (calcs(best_model,xbar[,1]) + calcs(best_model,xbar[,2]))
	if (!all(is.finite(model_stn))) { stop("model_stn contains non-finite values") }
	
	# calculate null distribution of model_stn using the baseline (i.e. best fit) condition
	if (best_idx==1) {
		# condition A is the baseline
		pValue = getSTNdistrib(A = A, nA = dim(A)[2], nB = dim(B)[2], iter = num_iter, xbar0_replace = min_xbar_value[best_idx], best_model = best_model) %>% getPvals(model_stn, num_iter, .)
	} else {
		# condition B is the baseline, so assign it to A (since A is used as baseline while resampling)
		pValue = getSTNdistrib(A = B, nA = dim(B)[2], nB = dim(A)[2], iter = num_iter, xbar0_replace = min_xbar_value[best_idx], best_model = best_model) %>% getPvals(model_stn, num_iter, .)
	}
	if (!all(is.finite(pValue))) { stop("pValue contains non-finite values") }
	
	return(list(model_stn = model_stn, pValue = pValue))
}

#' estimate the STN distribution
#'
#' as a dplyr data_frame
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_rows
#'
#' @export
getSTNdistrib = function(A, nA, nB, iter, xbar0_replace, best_model) {
	# A: data matrix for the baseline (i.e. best fit) condition
	# nA: number of replicates in the baseline condition
	# nB: number of replicates in the other condition
	cat("generating resampled STN distribution ...", "\n", sep="")
	XbarStar = vector("list", length=nrow(A))
	for (i in 1:nrow(A)) {
		if (i %% 100 == 0) { cat("i = ", i, "\n", sep="") }
		Astar_i = matrix(sample(A[i,], size=nA*iter, replace=TRUE), nrow=iter, ncol=nA)
		Bstar_i = matrix(sample(A[i,], size=nB*iter, replace=TRUE), nrow=iter, ncol=nB)
		XbarStar[[i]] = dplyr::data_frame(xbarA = rowMeans(Astar_i), xbarB = rowMeans(Bstar_i)) %>%
										dplyr::mutate(xbarA = ifelse(xbarA==0,xbar0_replace,xbarA),
																	xbarB = ifelse(xbarB==0,xbar0_replace,xbarB),
																	model_stn_dist = (xbarB - xbarA) / (calcs(best_model,xbarB) + calcs(best_model,xbarA)))
	}
	XbarStar = do.call("bind_rows", XbarStar)
	if (XbarStar %>% dplyr::filter(!is.finite(model_stn_dist)) %>% nrow() > 0) { stop("model_stn_dist contains non-finite values") }
	return(XbarStar %>% dplyr::transmute(model_stn_dist))
}

#' calculate p-values
#'
#' for each protein, using the generated STN distribution
#'
#' @export
getPvals = function(model_stn, iter, model_stn_dist) {
	cat("calculating p-values ...", "\n", sep="")
	pValue = rep(NA, length(model_stn))
	F = ecdf(unlist(model_stn_dist, use.names=FALSE))
	for (i in 1:length(model_stn)) {
		p = F(model_stn[i])
		pValue[i] = 2*min(p,1-p)
	}
	if (!all(is.finite(pValue))) { stop("pValue contains non-finite values") }
	return(pValue)
}

#' plot stn against p-value
#'
#' and stn for each protein
#'
#' @export
stn_pval_plots = function(stn_pval, outfile) {
	png(filename=outfile, width=800, height=400, units="px")
	par(mfrow=c(1,2), mar=c(5,4,4,2)+0.1)
	plot(stn_pval$model_stn, stn_pval$pValue, type="p", pch=20, col="gray", xlab="STN", ylab="p-value", main="P-value distribution")
	plot(sort(stn_pval$model_stn), type="p", pch=20, col="gray", xlab="protein", ylab="STN", main="STN")
	dev.off()
}

#' calculate diff-exp statss
#'
#' and return the results as a data.frame
#'
#' @export
diff_exp_table = function(stn_pval, Prot, num_digits) {
	ind = order(stn_pval$pValue)
	D = data.frame( Name = Prot[ind], STN = as.character(round(stn_pval$model_stn[ind],num_digits)), pVal = as.character(round(stn_pval$pValue[ind],num_digits)) )
	return(D)
}
