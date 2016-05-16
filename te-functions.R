library(np)

estimate.predictive.densities = function(Xt, Yt, p, kern.type = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', nmulti = 1, tol = 0.1, ftol = 0.1){
	tmp.df = data.frame(futureY = Yt[, dim(Yt)[2]], lagsX = Xt[, -dim(Xt)[2]], lagsY = Yt[, -dim(Yt)[2]])

	if (p > 1){
	  regressors.joint = c(paste("lagsX.", 1:(p), sep=""), paste("lagsY.", 1:(p), sep=""))
	  regressors.marginal = c(paste("lagsY.", 1:(p), sep=""))
	} else{
	  regressors.joint = c('lagsX', 'lagsY')
	  regressors.marginal = c('lagsY')
	}

	formula.for.np.joint = paste("futureY ~ ", paste(regressors.joint, collapse= "+"))
	formula.for.np.marginal = paste("futureY ~ ", paste(regressors.marginal, collapse= "+"))

	npcdensbw.out.joint = eval(parse(text = paste("npcdensbw(formula = ", formula.for.np.joint,", data = tmp.df, tol = tol, ftol = ftol, cxkertype = kern.type, cykertype = kern.type, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)")))
	npcdensbw.out.marginal = eval(parse(text = paste("npcdensbw(formula = ", formula.for.np.marginal,", data = tmp.df, tol = tol, ftol = ftol, cxkertype = kern.type, cykertype = kern.type, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)")))
  
	fitted.dens.joint = npcdens(npcdensbw.out.joint)
	fitted.dens.marginal = npcdens(npcdensbw.out.marginal)
	
	return(list(fitted.dens.marginal = fitted.dens.marginal, fitted.dens.joint = fitted.dens.joint))
}

estimate.predictive.density.marg = function(Yt, p, bw.use = NA, kern.type = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', nmulti = 1, tol = 0.1, ftol = 0.1){
	tmp.df = data.frame(futureY = Yt[, dim(Yt)[2]], lagsY = Yt[, -dim(Yt)[2]])

	if (p > 1){
	  regressors.marginal = c(paste("lagsY.", 1:(p), sep=""))
	} else{
	  regressors.marginal = c('lagsY')
	}
	
	formula.for.np.marginal = paste("futureY ~ ", paste(regressors.marginal, collapse= "+"))

	
  	
	if (is.na(bw.use)){
		npcdensbw.out.marginal = eval(parse(text = paste("npcdensbw(formula = ", formula.for.np.marginal,", data = tmp.df, tol = tol, ftol = ftol, cxkertype = kern.type, cykertype = kern.type, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)")))
		
		fitted.dens.marginal = npcdens(npcdensbw.out.marginal)
	}else{
		npcdensbw.out.marginal = rep(bw.use, dim(Yt)[2])
		fitted.dens.marginal = eval(parse(text = paste("npcdens(bws = npcdensbw.out.marginal, formula = ", formula.for.np.marginal,", data = tmp.df)")))
	}
	
	return(list(fitted.dens.marginal = fitted.dens.marginal, bw.marginal = npcdensbw.out.marginal))
}

estimate.predictive.density.joint = function(Xt, Yt, p, kern.type = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', nmulti = 1, tol = 0.1, ftol = 0.1){
	cat(sprintf('Using tols of %f / %f.', tol, ftol))
	
	tmp.df = data.frame(futureY = Yt[, dim(Yt)[2]], lagsX = Xt[, -dim(Xt)[2]], lagsY = Yt[, -dim(Yt)[2]])

	if (p > 1){
	  regressors.joint = c(paste("lagsX.", 1:(p), sep=""), paste("lagsY.", 1:(p), sep=""))
	} else{
	  regressors.joint = c('lagsX', 'lagsY')
	}

	formula.for.np.joint = paste("futureY ~ ", paste(regressors.joint, collapse= "+"))

	npcdensbw.out.joint = eval(parse(text = paste("npcdensbw(formula = ", formula.for.np.joint,", data = tmp.df, tol = tol, ftol = ftol, cxkertype = kern.type, cykertype = kern.type, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)")))
  
	fitted.dens.joint = npcdens(npcdensbw.out.joint)
	
	return(list(fitted.dens.joint = fitted.dens.joint, bw.joint = npcdensbw.out.joint))
}


filter.timeseries = function(Xt, Yt, fitted.dens.obj, p, dt.future = 5, future.offset = 50){
	future.eval = seq(min(Yt[, p+1])-future.offset, max(Yt[, p+1])+future.offset, by = dt.future)

	entropies.np = rep(NA, dim(Yt)[1]+p)
	variances.np = rep(NA, dim(Yt)[1]+p)
	means.np     = rep(NA, dim(Yt)[1]+p)
	modes.np     = rep(NA, dim(Yt)[1]+p)

	Dt = dt.future
	
		for (past.ind in 1:dim(Yt)[1]){
			
			if (is.na(sum(Xt))){
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsY = Xt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}else{
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsX.', p), 1:p, rep(' = Xt[past.ind, ', p), 1:p, rep('], ', p), rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}
  
		  extract.grid <- eval(parse(text = paste("expand.grid(futureY = future.eval,", expand.grid.opt, ")")))
  
		  fhat = predict(fitted.dens.obj, newdata = extract.grid)
  
		  entropy.conditional = -sum(fhat*log(fhat)*Dt, na.rm = TRUE)
  
		  entropies.np[past.ind+p] = entropy.conditional
  
		  mean.np = sum(fhat*future.eval*Dt)
		  var.np = sum(fhat*(future.eval - mean.np)^2*Dt)
  
		  means.np[past.ind+p] = mean.np
		  variances.np[past.ind+p] = var.np
  
		  if (compute.modes){
		    n.f = length(fhat)-1
		    deriv.signs = sign(diff(fhat))
		    modes.hat = deriv.signs[1:(n.f - 1)] - deriv.signs[2:n.f]
    
		    arg.mode = which(fhat == max(fhat[modes.hat == 2]))
    
		    modes.np[past.ind+p] = future.eval[arg.mode]
		  }
		}
	
	return(list(means.np = means.np, variances.np = variances.np, entropies.np = entropies.np))
}

filter.timeseries.er = function(Xt, Yt, fitted.dens.obj, p, dt.future = 5, future.offset = 50){
	future.eval = seq(min(Yt[, p+1])-future.offset, max(Yt[, p+1])+future.offset, by = dt.future)

	entropies.np = rep(NA, dim(Yt)[1]+p)

	Dt = dt.future
	
		for (past.ind in 1:dim(Yt)[1]){
			
			if (is.na(sum(Xt))){
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsY = Xt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}else{
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsX.', p), 1:p, rep(' = Xt[past.ind, ', p), 1:p, rep('], ', p), rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}
  
		  extract.grid <- eval(parse(text = paste("expand.grid(futureY = future.eval,", expand.grid.opt, ")")))
  
		  fhat = predict(fitted.dens.obj, newdata = extract.grid)
  
		  entropy.conditional = -sum(fhat*log(fhat)*Dt, na.rm = TRUE)
  
		  entropies.np[past.ind+p] = entropy.conditional
		}
	
	return(list(entropies.np = entropies.np))
}

filter.timeseries.er.integrate = function(Xt, Yt, fitted.dens.obj, p, dt.future = 5, future.offset = 50){
	future.eval = seq(min(Yt[, p+1])-future.offset, max(Yt[, p+1])+future.offset, by = dt.future)

	entropies.np = rep(NA, dim(Yt)[1]+p)

	Dt = dt.future
	
		for (past.ind in 1:dim(Yt)[1]){
			
			if (is.na(sum(Xt))){
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsY = Xt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}else{
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsX.', p), 1:p, rep(' = Xt[past.ind, ', p), 1:p, rep('], ', p), rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}
  
		  extract.grid <- eval(parse(text = paste("expand.grid(futureY = future.eval,", expand.grid.opt, ")")))
  
		  fhat = predict(fitted.dens.obj, newdata = extract.grid)
  
		  entropy.conditional = -sum(fhat*log(fhat)*Dt, na.rm = TRUE)
  
		  entropies.np[past.ind+p] = entropy.conditional
		}
	
	return(list(entropies.np = entropies.np, fhat = fhat, future.eval = future.eval))
}

filter.timeseries.er.loo = function(Xt, Yt, fitted.dens.obj, p){
	entropies.np = rep(NA, dim(Yt)[1]+p)
	
		for (loo.ind in 1:dim(Yt)[1]){
			if (is.na(sum(Xt))){
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = Yt[-loo.ind, -dim(Yt)[2]], tydat = Yt[-loo.ind, dim(Yt)[2]], exdat = matrix(Yt[loo.ind, -dim(Yt)[2]], ncol = p), eydat = Yt[loo.ind, dim(Yt)[2]])
			}else{
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = cbind(Xt[-loo.ind, -dim(Xt)[2]], Yt[-loo.ind, -dim(Yt)[2]]), tydat = Yt[-loo.ind, dim(Yt)[2]], exdat = matrix(c(Xt[loo.ind, -dim(Xt)[2]], Yt[loo.ind, -dim(Yt)[2]]), ncol = p + p), eydat = Yt[loo.ind, dim(Yt)[2]])
			}
			
			entropies.np[loo.ind + p] = -log(fitted(npcdens.out))
		}
	
	return(list(entropies.np = entropies.np))
}

filter.timeseries.er.lwo = function(Xt, Yt, fitted.dens.obj, p, window.length = 2){
	entropies.np = rep(NA, dim(Yt)[1]+p)
	
		for (loo.ind in 1:dim(Yt)[1]){
			loo.inds = max(c(1, loo.ind - window.length/2)):min(c(dim(Yt)[1], loo.ind + window.length/2))
			
			if (is.na(sum(Xt))){
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = Yt[-loo.inds, -dim(Yt)[2]], tydat = Yt[-loo.inds, dim(Yt)[2]], exdat = matrix(Yt[loo.ind, -dim(Yt)[2]], ncol = p), eydat = Yt[loo.ind, dim(Yt)[2]])
			}else{
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = cbind(Xt[-loo.inds, -dim(Xt)[2]], Yt[-loo.inds, -dim(Yt)[2]]), tydat = Yt[-loo.inds, dim(Yt)[2]], exdat = matrix(c(Xt[loo.ind, -dim(Xt)[2]], Yt[loo.ind, -dim(Yt)[2]]), ncol = p + p), eydat = Yt[loo.ind, dim(Yt)[2]])
			}
			
			entropies.np[loo.ind + p] = -log(fitted(npcdens.out))
		}
	
	return(list(entropies.np = entropies.np, npcdens.out.last = npcdens.out))
}

filter.timeseries.er.loo.integrate = function(Xt, Yt, fitted.dens.obj, p, dt.future = 1, future.offset = 50){
	future.eval = seq(min(Yt[, p+1])-future.offset, max(Yt[, p+1])+future.offset, by = dt.future)
	
	entropies.np = rep(NA, dim(Yt)[1]+p)
	
		for (past.ind in 1:dim(Yt)[1]){
			if (is.na(sum(Xt))){
	  		  if (p == 1){
	  		    # 0/0
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = Yt[-past.ind, -dim(Yt)[2]], tydat = Yt[-past.ind, dim(Yt)[2]], exdat = matrix(Yt[past.ind, -dim(Yt)[2]], ncol = p), eydat = Yt[past.ind, dim(Yt)[2]])
				
				expand.grid.opt = paste0('lagsY = Yt[past.ind, 1]')
	  		  }else{
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = Yt[-past.ind, -dim(Yt)[2]], tydat = Yt[-past.ind, dim(Yt)[2]], exdat = matrix(Yt[past.ind, -dim(Yt)[2]], ncol = p), eydat = Yt[past.ind, dim(Yt)[2]])
	  		    expand.grid.opt = paste0(rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}else{
	  		  if (p == 1){
	  		    # 0/0
				
				npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = cbind(Xt[-past.ind, -dim(Xt)[2]], Yt[-past.ind, -dim(Yt)[2]]), tydat = Yt[-past.ind, dim(Yt)[2]], exdat = matrix(c(Xt[past.ind, -dim(Xt)[2]], Yt[past.ind, -dim(Yt)[2]]), ncol = p + p), eydat = Yt[past.ind, dim(Yt)[2]])
				
				expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]') # was: expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]')
	  		  }else{
				 npcdens.out = npcdens(bws = fitted.dens.obj$bws, txdat = cbind(Xt[-past.ind, -dim(Xt)[2]], Yt[-past.ind, -dim(Yt)[2]]), tydat = Yt[-past.ind, dim(Yt)[2]], exdat = matrix(c(Xt[past.ind, -dim(Xt)[2]], Yt[past.ind, -dim(Yt)[2]]), ncol = p + p), eydat = Yt[past.ind, dim(Yt)[2]])
	  		    expand.grid.opt = paste0(rep('lagsX.', p), 1:p, rep(' = Xt[past.ind, ', p), 1:p, rep('], ', p), rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}
  
		  extract.grid <- eval(parse(text = paste("expand.grid(futureY = future.eval,", expand.grid.opt, ")")))
  
		  fhat = predict(npcdens.out, newdata = extract.grid)
  
		  entropy.conditional = -sum(fhat*log(fhat)*dt.future, na.rm = TRUE)
  
		  entropies.np[past.ind+p] = entropy.conditional
		}
	
	return(list(entropies.np = entropies.np, fhat = fhat, future.eval = future.eval))
}

get.predictive.density = function(Xt, Yt, time.ind, fitted.dens.obj, p, dt.future = 1, future.offset = 50){
	future.eval = seq(min(Yt[, p+1])-future.offset, max(Yt[, p+1])+future.offset, by = dt.future)

	Dt = dt.future
	
	past.ind = time.ind
			
			if (is.na(sum(Xt))){
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsY = Yt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}else{
	  		  if (p == 1){
	  		    expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]') # was: expand.grid.opt = paste0('lagsX = Xt[past.ind, 1], lagsY = Yt[past.ind, 1]')
	  		  }else{
	  		    expand.grid.opt = paste0(rep('lagsX.', p), 1:p, rep(' = Xt[past.ind, ', p), 1:p, rep('], ', p), rep('lagsY.', p), 1:p, rep(' = Yt[past.ind, ', p), 1:p, rep(']', p), collapse = ',')
	  		  }
			}
  
		  extract.grid <- eval(parse(text = paste("expand.grid(futureY = future.eval,", expand.grid.opt, ")")))
  
		  fhat = predict(fitted.dens.obj, newdata = extract.grid)
	
		return(list(future.eval = future.eval, fhat = fhat))
}