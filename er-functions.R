library(forecast)
library(extrafont)
library(np)
library(pdc) # Needed for entropyHeuristic

estimate.pdfs = function(data, kernel.type = 'uniform', bw = NA, nmulti = 1){
  if (is.na(bw)){
	  bw.out.marg = npudensbw(dat = data[, -dim(data)[2]], bwmethod = 'cv.ml', nmulti = nmulti, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  bw.out.full  = npudensbw(dat = data[, ], bwmethod = 'cv.ml', nmulti = nmulti, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  
	  bws.marg = bw.out.marg$bw
	  bws.full = bw.out.full$bw
  }else{
	  bws.marg = rep(bw, dim(data)[2] - 1)
	  bws.full = rep(bw, dim(data)[2])
  }

  return(list(bw.out.marg = bw.out.marg, bw.out.full = bw.out.full))
}

entropy.rate.np.diff.loo = function(data, kernel.type = 'uniform', bw = NA){
    full.entropy.all = rep(NA, dim(data)[1])
    marg.entropy.all = rep(NA, dim(data)[1])

  if (is.na(bw)){
	  bw.out.marg = npudensbw(dat = data[, -dim(data)[2]], bwmethod = 'cv.ml', nmulti = 4, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  bw.out.full  = npudensbw(dat = data[, ], bwmethod = 'cv.ml', nmulti = 4, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  
	  bws.marg = bw.out.marg$bw
	  bws.full = bw.out.full$bw
  }else{
	  bws.marg = rep(bw, dim(data)[2] - 1)
	  bws.full = rep(bw, dim(data)[2])
  }

  for (loo.ind in 1:(dim(data)[1])){
    np.out.marg = npudens(tdat = data[-loo.ind, -dim(data)[2]], bws = bws.marg)
    np.out.full = npudens(tdat = data[-loo.ind, ], bws = bws.full)

    dens.out.marg = predict(np.out.marg, edat = t(as.matrix(data[loo.ind, -dim(data)[2]])))
    dens.out.full = predict(np.out.full, edat = t(as.matrix(data[loo.ind, ])))

    if (dens.out.full != 0){
      full.entropy.all[loo.ind] = -log(dens.out.full)
    }

    if (dens.out.marg != 0){
      marg.entropy.all[loo.ind] = -log(dens.out.marg)
    }
  }

 entropy.rate = mean(full.entropy.all) - mean(marg.entropy.all)

  return(list(bw.out.marg = bws.marg, bw.out.full = bws.full, entropy.rate = entropy.rate, marg.entropy.all = marg.entropy.all, full.entropy.all = full.entropy.all))
}

entropy.rate.np.diff.loi = function(data, kernel.type = 'uniform', bw = NA){
  full.entropy.all = rep(NA, dim(data)[1])
  marg.entropy.all = rep(NA, dim(data)[1])
	
	
	
  if (is.na(bw)){
	  bw.out.marg = npudensbw(dat = data[, -dim(data)[2]], bwmethod = 'cv.ml', nmulti = 4, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  bw.out.full  = npudensbw(dat = data[, ], bwmethod = 'cv.ml', nmulti = 4, tol = 0.1, ftol = 0.1, ckertype = kernel.type)
	  
	  bws.marg = bw.out.marg$bw
	  bws.full = bw.out.full$bw
  }else{
	  bws.marg = rep(bw, dim(data)[2] - 1)
	  bws.full = rep(bw, dim(data)[2])
  }

  np.out.marg = npudens(tdat = data[, -dim(data)[2]], bws = bws.marg, ckertype = kernel.type)
  np.out.full = npudens(tdat = data[, ], bws = bws.full, ckertype = kernel.type)

  for (loo.ind in 1:(dim(data)[1])){
    dens.out.marg = predict(np.out.marg, edat = t(as.matrix(data[loo.ind, -dim(data)[2]])))
    dens.out.full = predict(np.out.full, edat = t(as.matrix(data[loo.ind, ])))

    if (dens.out.full != 0){
      full.entropy.all[loo.ind] = -log(dens.out.full)
    }

    if (dens.out.marg != 0){
      marg.entropy.all[loo.ind] = -log(dens.out.marg)
    }
  }

  entropy.rate = mean(full.entropy.all) - mean(marg.entropy.all)

  return(list(bw.out.marg = bws.marg, bw.out.full = bws.full, entropy.rate = entropy.rate, marg.entropy.all = marg.entropy.all, full.entropy.all = full.entropy.all))
}

entropy.rate.np.ave = function(data){
  full.entropy = 0.
  
  bw.out.full  = npudensbw(dat = data[, ], bwmethod = 'cv.ml', nmulti = 4, tol = 0.1, ftol = 0.1)
  
  for (loo.ind in 1:(dim(data)[1])){
    np.out.full = npudens(tdat = data[-loo.ind, ], bws = bw.out.full$bw)
    
    dens.out.full = predict(np.out.full, edat = t(as.matrix(data[loo.ind, ])))
    
    if (dens.out.full != 0){
      full.entropy = full.entropy - log(dens.out.full)
    }
  }
  
  entropy.rate = (full.entropy)/(dim(data)[1]*dim(data)[2])
  
  return(list(bw.out.full = bw.out.full, entropy.rate = entropy.rate))
}

embed.ts = function(ts, L = 2){
  n = length(ts)
  ts.embedded = matrix(0, nrow = n - (L - 1), ncol = L)
  
  for (start.ind in 1:L){
    ts.embedded[, start.ind] = ts[((start.ind-1) + 1):(n - (L - 1) + (start.ind-1))]
  }
  
  return(ts.embedded)
}

mv.epa.kernel = function(u){
	p = length(u)
	
	us = sum(u*u)
	if (us <= 1){
		return(p*(p+2)*gamma(p/2)/(4*pi^(p/2)) * (1 - us))
	}else{
		return(0)
	}
}