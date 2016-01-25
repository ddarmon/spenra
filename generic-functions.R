uniformize = function(x){
	x.transformed = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
	
	return(x.transformed)
}

gaussian.kernel = function(u){
	return(exp(-u^2/2)/sqrt(2*pi))
}

embed.ts = function(ts, L = 2){
  n = length(ts)
  ts.embedded = matrix(0, nrow = n - (L - 1), ncol = L)
  
  for (start.ind in 1:L){
    ts.embedded[, start.ind] = ts[((start.ind-1) + 1):(n - (L - 1) + (start.ind-1))]
  }
  
  return(ts.embedded)
}