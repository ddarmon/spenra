#' Transform to Uniform Scale
#' 
#' 
#' Transform a vector with support on [a, b] to have
#' support on [0, 1].
#' @param x A vector of values.
#' @export

uniformize = function(x){
	x.transformed = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
	
	return(x.transformed)
}

#' Gaussian Kernel of Scalar
#' 
#' Returns the standard Gaussian kernel evaluated at u.
#' @param u Evaluation points
#' @keywords kernel
#' @export

gaussian.kernel = function(u){
	return(exp(-u^2/2)/sqrt(2*pi))
}

#' Embed a Scalar Time Series
#' 
#' Embeds a scalar time series in a regression-like matrix.
#' @param ts The scalar time series.
#' @param L The block length used to embed.
#' @keywords time series, embedding, regression
#' @export

embed.ts = function(ts, L = 2){
  n = length(ts)
  ts.embedded = matrix(0, nrow = n - (L - 1), ncol = L)
  
  for (start.ind in 1:L){
    ts.embedded[, start.ind] = ts[((start.ind-1) + 1):(n - (L - 1) + (start.ind-1))]
  }
  
  return(ts.embedded)
}