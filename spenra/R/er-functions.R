#' Estimate Predictive / Conditional Density
#' 
#' estimate.pred.dens estimates the predictive / conditional density of
#' order p associated with the scalar time series x using np.
#' @param x A continuous-valued scalar time series.
#' @param p The order of the conditional density.
#' @param kerntype The kernel type, one of {'guassian', 'uniform', 'epanechnikov'}.
#' @param bwtype The bandwidth type, one of {'fixed', 'generalized_nn', 'adaptive_nn'}.
#' @param bwmethod The method used to choose the bandwidths, one of {cv.ml, cv.ls, normal-reference}. Defaults to cv.ml to perform LOOCV on the specific entropy rate.
#' @param bwuse Specify a fixed bandwidth to use across all lags, or NA to estimate the bandwidths in a data-driven way.
#' @param nmulti The number of random initial conditions to use by the optimizer for choosing the bandwidths.
#' @keywords predictive density, conditional density, kernel density estimator
#' @export

estimate.pred.dens = function(x, p, kerntype = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', bwuse = NA, nmulti = 1){
	Xt = embed.ts(x, p+1)
	
	npcdensbw.out = np::npcdensbw(ydat = as.matrix(Xt[, p+1]), xdat = as.matrix(Xt[, 1:p]), tol = 0.1, ftol = 0.1, cxkertype = kerntype, cykertype = kerntype, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)
	npcdens.out = np::npcdens(npcdensbw.out)
	
	return(npcdens.out)
}

#' Estimate the Entropy Rate of a Scalar Time Series Using a Plug-In Leave-Window-Out Estimator
#' 
#' timeseries.entropyrate.lwo estimates the specific entropy rate using a plug-in, leave-window-out estimator.
#' @param x The scalar time series.
#' @param pred.dens The npcdens-based predictive density used to estimate the specific entropy rate.
#' @param half.window.length The half-length of the window about the present removed in estimating the specific entropy rate.
#' @export

timeseries.entropyrate.lwo = function(x, pred.dens, half.window.length = 1){
	p = length(pred.dens$xbw)
	
	Xt = embed.ts(x, p+1)
	
	entropyrate = rep(NA, length(x))
	
		for (loo.ind in 1:dim(Xt)[1]){
			loo.inds = max(c(1, loo.ind - half.window.length)):min(c(dim(Xt)[1], loo.ind + half.window.length))
			
			npcdens.out = npcdens(bws = pred.dens$bws, txdat = Xt[-loo.inds, -dim(Xt)[2]], tydat = Xt[-loo.inds, dim(Xt)[2]], exdat = matrix(Xt[loo.ind, -dim(Xt)[2]], ncol = p), eydat = Xt[loo.ind, dim(Xt)[2]])
			
			entropyrate[loo.ind + p] = -log(fitted(npcdens.out))
		}
	
	return(entropyrate = entropyrate)
}



choose.p.lwo = function(x, p.max, p.min = 1, half.window.length = 1, kerntype = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', bwuse = NA, nmulti = 1){
	ps = seq(p.min, p.max, 1)
	
	pred.dens.by.p = list()
	er.by.p = rep(NA)
	
	for (p.ind in ps){
		p = ps[p.ind]

		pred.dens = estimate.pred.dens(x, p, kerntype, bwtype, bwmethod, bwuse, nmulti)
		
		entropyrate = timeseries.entropyrate.lwo(x, pred.dens, half.window.length)
		
		er.by.p[p.ind] = mean(entropyrate[(p.max + 1):length(entropyrate)])
		
		pred.dens.by.p[[as.character(p.ind)]] = pred.dens
	}
	
	p.best = ps[which.min(er.by.p)]
	
	return(list(ps = ps, er.by.p = er.by.p, pred.dens.by.p = pred.dens.by.p, p.best = p.best, pred.dens = pred.dens.by.p[[as.character(p.best)]]))
}