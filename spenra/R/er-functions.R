#' Sentence-like Title for Function
#' 
#' Short description of what function does.
#' @param p1 Description of p1
#' @keywords keyword1, keyword2, keyword3
#' @export
#' function_name()

estimate.cond.dens = function(x, p, kerntype = 'gaussian', bwtype = 'fixed', bwmethod = 'cv.ml', bwuse = NA, nmulti = 1){
	cat(sprintf('I am here.'))
	Xt = embed.ts(x, p+1)
	
	show(Xt)
	
	npcdens.out = np::npcdensbw(ydat = as.matrix(Xt[, p+1]), xdat = as.matrix(Xt[, 1:p]), tol = 0.1, ftol = 0.1, cxkertype = kerntype, cykertype = kerntype, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)
	
	return(npcdens.out)
}