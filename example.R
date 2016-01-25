setwd('/Users/complexity/Documents/R/hrv-analysis')

library(statmod)
library(extrafont)
loadfonts()

plot.new()

source('/Users/complexity/Documents/R/hrv-analysis/ig-functions.R')
source('/Users/complexity/Documents/R/hrv-analysis/er-functions.R')
source('/Users/complexity/Documents/R/hrv-analysis/hrv-load-functions.R')

source('/Users/complexity/Documents/R/generic-functions.R')

# data.type = 'HDIG'
# data.type = 'rapp'
# data.type = 'tirp'
# data.type = 'tirp2'
# data.type = 'rachel'
# data.type = 'lorenz'
# data.type = 'physionet'
# data.type = 'sine'
# data.type = 'condstat'
# data.type = 'mixp'
# data.type = 'lar'
# data.type = 'rossler'
# data.type = 'lorenz+rossler'
# data.type = 'lorenz+rossler+lorenz'
data.type = 'nb'
# data.type = 'nb-mRSA'

entropy.type = 'er'
# entropy.type = 'ner'

front.back.cutoff = 50
# full.width = 120
full.width = 60

# ckertype = 'gaussian'
ckertype = 'uniform'

# for (tirp.number in c(1, 2, 3, 6, 7, 8)){
# for (tirp.number in c(1, 2, 3, 6, 7, 8)){
for (tirp.number in c(2)){
# for (tirp.number in c(1)){
# for (tirp.number in c(1, 4, 5)){

session = 1

# window.length = 100
window.length = 0

if ((data.type == 'tirp') | (data.type == 'nb') | (data.type == 'nb-mRSA')){
  figure.suffix = sprintf('-%g-%g', tirp.number, session)
}else{
  figure.suffix = ''
}

load.ts.out = load.ts(data.type, tirp.number = tirp.number, session = session)

Xt.all = load.ts.out$Xt
fname = load.ts.out$fname

beat.times.all = cumsum(Xt.all)

if ((data.type == 'nb') & (tirp.number == 6) & (session == 1)){
  Xt.all[1148] = Xt.all[1147]
}

segments = list()

if ((data.type == 'nb') | (data.type == 'nb-mRSA')){
  if (session == 1){
    marker.data = read.csv(sprintf('/Users/complexity/Documents/R/hrv-analysis/data/NB0%g_Respiration.txt', tirp.number), header = FALSE, sep = "")
  }else{
    marker.data = read.csv(sprintf('/Users/complexity/Documents/R/hrv-analysis/data/NB0%g_V2_Respiration.txt', tirp.number), header = FALSE, sep = "")
  }
  
  change.times = marker.data$V1[which(marker.data$V3 == 1)]
  
  if (session == 1)  {
    interval.to.times = read.csv(sprintf('/Users/complexity/Documents/R/hrv-analysis/data/NB0%g_RRinterval.txt', tirp.number), header = FALSE, sep = "")
  }else{
    interval.to.times = read.csv(sprintf('/Users/complexity/Documents/R/hrv-analysis/data/NB0%g_v2_RRinterval.txt', tirp.number), header = FALSE, sep = "")
  }
  
  
  for (segment.ind in 1:(length(change.times)-1)){
    segments[[segment.ind]] = c(which(interval.to.times[, 1] - change.times[segment.ind] > 0)[1], which(interval.to.times[, 1] - change.times[segment.ind + 1] > 0)[1] - 1)
  }
  
  segments[[length(change.times)]] = c(which(interval.to.times[, 1] - change.times[length(change.times)] > 0)[1], length(Xt.all))
} else if (data.type == 'lorenz+rossler+lorenz'){
  segments[[1]] = c(1, 500)
  segments[[2]] = c(501, 1000)
  segments[[3]] = c(1001, 1500)
}

if ((data.type == 'nb') & (tirp.number == 1) & (session == 2)){
  segments[[7]] = segments[[8]]
  segments[[8]] = c()
}

if ((data.type == 'nb') & (tirp.number == 4) & (session == 2)){
  segments[[1]][1] = 175
}

if ((data.type == 'nb') & (tirp.number == 6) & (session == 1)){
  segments[[4]][1] = 1150
}

beat.times.for.er = beat.times.all[1:(segments[[1]][1] - 1)]

for (segment.ind in 1:length(segments)){
  beat.times.for.er = c(beat.times.for.er, beat.times.all[segments[[segment.ind]][1]:segments[[segment.ind]][2]])
}

# smooth.mean.var = 30
for (smooth.mean.var in c(full.width/2)){
# for (smooth.mean.var in seq(10, 120, by = 10)){
# for (smooth.mean.var in seq(30, 50, by = 10)){
# for (smooth.mean.var in c(30)){
# for (smooth.mean.var in seq(10, 30, by = 10)){
# for (smooth.mean.var in c(10)){
# for (smooth.mean.var in c(20)){

  kde.by.segments = list()
  er.by.segments = list()
  sd.by.segments = list()
  Xt.by.segments = list()
  logvar.by.segments = list()
  mean.by.segments = list()
  smooth.er.by.segments = list()
  beat.times.by.segments = list()
  
apen.by.segments = rep(0, length(segments))

ps.savename = sprintf('/Users/complexity/Documents/R/hrv-analysis/ps-bysegment-saved-lwo%g-adaptive-nn/%s.dat', window.length, strsplit(tail(strsplit(fname, '/')[[1]], n = 1)[1], split = '\\.')[[1]][1])

ps.by.segment = read.csv(ps.savename)$x

design.points = 1:length(Xt.all)

beat.mean.all = fitted(npreg(bws = c(smooth.mean.var), txdat = design.points, tydat = Xt.all))

beat.var.all = fitted(npreg(bws = c(smooth.mean.var), txdat = design.points, tydat = (Xt.all - beat.mean.all)^2))
beat.logvar.all = log(beat.var.all)

temporal.mean.all = fitted(npreg(bws = c(smooth.mean.var), txdat = beat.times.all, tydat = Xt.all, ckertype = ckertype))

temporal.var.all = fitted(npreg(bws = c(2*smooth.mean.var), txdat = beat.times.all, tydat = (Xt.all - temporal.mean.all)^2, ckertype = ckertype))
temporal.logvar.all = log(temporal.var.all)

design.points.all = 1:(segments[[1]][1] - 1)
er.all = c(rep(NA, segments[[1]][1] - 1))

for (segment.ind in 1:length(segments)){
  inds.for.segment = segments[[segment.ind]]
  
  design.points.all = c(design.points.all, inds.for.segment[1]:inds.for.segment[2])
  
  beat.times.by.segments[[segment.ind]] = beat.times.all[inds.for.segment[1]:inds.for.segment[2]]
  
  Xt = Xt.all[inds.for.segment[1]:inds.for.segment[2]]
  
  # The bounds for integration.
  
  lower.bound = min(Xt) - IQR(Xt)
  upper.bound = max(Xt) + IQR(Xt)
  
  Xt.by.segments[[segment.ind]] = Xt
  
  sd.by.segments[[segment.ind]] = sd(Xt)
  
  # design.points = 1:length(Xt)
  design.points = beat.times.all[inds.for.segment[1]:inds.for.segment[2]]
  
  temporal.mean = fitted(npreg(bws = c(smooth.mean.var), txdat = beat.times.by.segments[[segment.ind]], tydat = Xt, exdat = beat.times.by.segments[[segment.ind]], ckertype = ckertype))
  
#   temporal.logvar = fitted(npreg(bws = c(smooth.mean.var), txdat = design.points, tydat = log((Xt - temporal.mean)^2)))
#   temporal.var = exp(temporal.logvar)
  
  temporal.var = fitted(npreg(bws = c(2*smooth.mean.var), txdat = beat.times.by.segments[[segment.ind]], tydat = (Xt - temporal.mean)^2, exdat = beat.times.by.segments[[segment.ind]], ckertype = ckertype))
  temporal.logvar = log(temporal.var)
  
  if (front.back.cutoff > 0){
    temporal.logvar[1:front.back.cutoff] = NA
    temporal.logvar[(length(temporal.logvar) - (front.back.cutoff-1)):(length(temporal.logvar))] = NA
  }
  
  mean.by.segments[[segment.ind]] = temporal.mean
  
  logvar.by.segments[[segment.ind]] = temporal.logvar
  
  p = ps.by.segment[segment.ind]
  # p = 5
  
  # kern.type = 'epanechnikov'
  kern.type = 'gaussian'
  # kern.type = 'uniform'
  
  # bwtype = 'fixed'
  bwtype = 'adaptive_nn'
  # bwtype = 'generalized_nn'
  
  bwmethod = 'cv.ml'
  # bwmethod = 'normal-reference'
  
  bw.use = NA
  # bw.use = 0.2*sd(Xt)
  
  nmulti = 1
  # nmulti = 4
  
  # smooth.df = floor(0.9*length(Xt))
  smooth.df = length(Xt)
  
  include.time = FALSE
  compute.modes = TRUE
  to.sim = FALSE
  
  ts.embed = embed.ts(Xt, L = p+1)
  
  p.apen = p
  
  ts.embed.apen = embed.ts(Xt, L = p.apen+1)
  
  npcdens.apen = npcdens(bws = rep(0.2*sd(Xt), p.apen+1), txdat = ts.embed.apen[, 1:p.apen], tydat = ts.embed.apen[, p.apen+1], cxkertype = 'uniform', cykertype = 'uniform')
  
  if (entropy.type == 'er'){
    er.apen = -log(fitted(npcdens.apen))
  } else{
    er.apen = 0.5*log(2*pi*exp(1)*var(Xt))+log(fitted(npcdens.apen))
  }
  
  
  apen.by.segments[segment.ind] = mean(er.apen, na.rm = TRUE)
  
  x = ts.embed[, p+1]
  
  design.points = 1:length(Xt)
  
  # temporal.mean = smooth.spline(Xt, df = smooth.df)$y
  # 
  # temporal.var = smooth.spline((Xt - temporal.mean)^2, df = smooth.df)$y
  # temporal.var[temporal.var < 0] = 0
  # 
  # temporal.sd = sqrt(temporal.var)
  
  library(np)
  
  par(mfrow = c(1, 1))
  
  time.points = (p+1):length(Xt)
  
  if (include.time){
    tmp.df = data.frame(future = ts.embed[, dim(ts.embed)[2]], lags = ts.embed[, -dim(ts.embed)[2]], time.points = time.points)
  } else{
    tmp.df = data.frame(future = ts.embed[, dim(ts.embed)[2]], lags = ts.embed[, -dim(ts.embed)[2]])
  }
  
  if (p > 1){
    regressors = paste("lags.", 1:(p), sep="")
  } else{
    regressors = c('lags')
  }
  
  if (include.time){
    formula.for.np = paste("future ~ ", paste(c(regressors, 'time.points'), collapse= "+"))
  }else{
    formula.for.np = paste("future ~ ", paste(regressors, collapse= "+"))
  }
  
  dens.savename = sprintf('/Users/complexity/Documents/R/hrv-analysis/npcdens-saved-adaptive-nn/%s-seg%g-p%g.RData', strsplit(tail(strsplit(fname, '/')[[1]], n = 1)[1], split = '\\.')[[1]][1], segment.ind, p)
  
  dir.save = sprintf('/Users/complexity/Documents/R/hrv-analysis/er-saved-lwo%g', window.length)
  
  if(!file.exists(dir.save)){
    dir.create(dir.save)
  }
  
  er.savename = sprintf('/Users/complexity/Documents/R/hrv-analysis/er-saved-lwo%g-adaptive-nn/%s-seg%g-p%g.dat', window.length, strsplit(tail(strsplit(fname, '/')[[1]], n = 1)[1], split = '\\.')[[1]][1], segment.ind, p)
  
  if (file.exists(dens.savename)){
  # if(FALSE){
    load(dens.savename)
    
    fitted.dens = npcdens(npcdensbw.out)
  }else{
    if (is.na(bw.use)){
      npcdensbw.out = eval(parse(text = paste("npcdensbw(formula = ", formula.for.np,", data = tmp.df, tol = 0.1, ftol = 0.1, cxkertype = kern.type, cykertype = kern.type, bwtype = bwtype, bwmethod = bwmethod, nmulti = nmulti)")))
      
      fitted.dens = npcdens(npcdensbw.out)
      
      save(npcdensbw.out, file = dens.savename)
    }else{
      fitted.dens = eval(parse(text = paste("npcdens(formula = ", formula.for.np,", data = tmp.df, cxkertype = kern.type, cykertype = kern.type, bws = rep(bw.use, dim(ts.embed)[2]))")))
    }
  }
  
  plot(fitted.dens)
  
  cat('\n\n', fitted.dens$bws$xbw, fitted.dens$bws$ybw, '\n\n')
  
  kde.by.segments[[segment.ind]] = fitted.dens
  
  if (file.exists(er.savename)){
  # if (FALSE){
    entropies.np = read.csv(er.savename)$x
  } else{
    entropies.np = rep(NA, length(Xt))
    
    for (past.ind in (p+1):length(Xt)){
      if((past.ind %% 100 == 1)){
        cat(sprintf('\n\nFiltering on time point %g.\n\n', past.ind))
      }
      # time.ind = c(past.ind - p)
      
      time.ind = max(c(1, past.ind - p - window.length/2)):min(c(length(Xt) - p, past.ind - p + window.length/2))
      
      tmp.df.loo = tmp.df[-time.ind, ]
      
      if (p == 1){
        expand.grid.opt = paste0('lags = Xt[past.ind-', p, ']')
      }else{
        expand.grid.opt = paste0(rep('lags.', p), 1:p, rep(' = Xt[past.ind-', p), p:1, rep(']', p), collapse = ',')
      }
      
      # fitted.dens.loo = eval(parse(text = paste("npcdens(formula = ", formula.for.np,", data = tmp.df.loo, cxkertype = kern.type, cykertype = kern.type, bws = c(fitted.dens$bws$ybw, fitted.dens$bws$xbw))")))
      fitted.dens.loo = eval(parse(text = paste("npcdens(formula = ", formula.for.np,", data = tmp.df.loo, cxkertype = kern.type, cykertype = kern.type, bws = fitted.dens$bws)")))
      
      integrand = function(future.eval, Xt, fitted.dens, p, past.ind){
        if (p == 1){
          expand.grid.opt = paste0('lags = Xt[past.ind-', p, ']')
        }else{
          expand.grid.opt = paste0(rep('lags.', p), 1:p, rep(' = Xt[past.ind-', p), p:1, rep(']', p), collapse = ',')
        }
        
        extract.grid <- eval(parse(text = paste("expand.grid(future = future.eval, time.points = c(past.ind), ", expand.grid.opt, ")")))
        
        fhat = predict(fitted.dens, newdata = extract.grid)
        
        entropy = -log(fhat)*fhat
        
        entropy[is.na(entropy)] = 0
        # entropy[(entropy > 0) & (entropy < 1e-50)] = 0
      
        return(entropy)
      }
      
      fhat = function(future.eval, Xt, fitted.dens, p, past.ind){
        if (p == 1){
          expand.grid.opt = paste0('lags = Xt[past.ind-', p, ']')
        }else{
          expand.grid.opt = paste0(rep('lags.', p), 1:p, rep(' = Xt[past.ind-', p), p:1, rep(']', p), collapse = ',')
        }
        
        extract.grid <- eval(parse(text = paste("expand.grid(future = future.eval, time.points = c(past.ind), ", expand.grid.opt, ")")))
        
        fhat = predict(fitted.dens, newdata = extract.grid)
        
        # fhat[fhat < 1e-5] = 0
        
        return(list(fhat = fhat, ys = future.eval))
      }
      
      fhat.integrand = function(future.eval, Xt, fitted.dens, p, past.ind){
        if (p == 1){
          expand.grid.opt = paste0('lags = Xt[past.ind-', p, ']')
        }else{
          expand.grid.opt = paste0(rep('lags.', p), 1:p, rep(' = Xt[past.ind-', p), p:1, rep(']', p), collapse = ',')
        }
        
        extract.grid <- eval(parse(text = paste("expand.grid(future = future.eval, time.points = c(past.ind), ", expand.grid.opt, ")")))
        
        fhat = predict(fitted.dens, newdata = extract.grid)
        
        # fhat[fhat < 1e-5] = 0
        
        return(fhat)
      }
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #
      # TURN OFF BETWEEN
      #
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
#       fhat.out = fhat(seq(min(Xt) - IQR(Xt), max(Xt) + IQR(Xt), by = 0.001), Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind)
#     
#       fhat.integrand.out = integrate(fhat.integrand, min(Xt) - IQR(Xt), max(Xt) + IQR(Xt), rel.tol = 0.01, Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind)
#       
#       par(mfrow = c(1, 1))  
#       plot(fhat.out$ys, fhat.out$fhat, type = 'l', ylim = c(0, 40))
#       plot(fhat.out$ys, -log(fhat.out$fhat)*fhat.out$fhat, type = 'l')
#       abline(h = 0, col = 'red')
#       
#       integrate(integrand, 0, 1.4, rel.tol = 0.01, Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind)$value
#       integrate(integrand, min(Xt) - IQR(Xt), max(Xt) + IQR(Xt), rel.tol = 0.01, Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind)$value
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      #
      # TURN OFF BETWEEN
      #
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      # entropies.np[past.ind] = integrate(integrand, -Inf, Inf, rel.tol = 0.01, Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind)$value
      entropies.np[past.ind] = integrate(integrand, lower.bound, upper.bound, rel.tol = 0.01, Xt = Xt, fitted.dens = fitted.dens.loo, p = p, past.ind = past.ind, subdivisions = 1000)$value
    }
    
    write.csv(x = entropies.np, file = er.savename, quote = FALSE, row.names = FALSE)
  }
  
  er.all = c(er.all, entropies.np)
  
  smooth.er = c(rep(NA, p), fitted(npreg(bws = c(smooth.mean.var), txdat = design.points, tydat = entropies.np, ckertype = ckertype)))
  
  smooth.er.by.segments[[segment.ind]] = smooth.er
  
  if (entropy.type == 'er'){
    er.by.segments[[segment.ind]] = entropies.np
  } else{
    er.by.segments[[segment.ind]] = (0.5*log(2*pi*exp(1)*var(Xt)) - entropies.np)
  }
}


fig.fname = sprintf('%s%s-pCV-%s-by-time-segmented.pdf', data.type, figure.suffix, entropy.type)
pdf(fig.fname, family = "CMU Serif", width = 15)
par(mfrow = c(1, 1), mar = c(5, 5, 2, 1), cex.lab = 1.5, cex.axis = 1.5)

if (entropy.type == 'er'){
  if (data.type == 'nb'){
    ylims = c(-4, 0) # entropy rate for nb files
  }else{
    ylims = c(-1.5, 2) # entropy rate for mixed lorenz + rossler
  }
  
  plot(0, 0, xlim = c(1, length(Xt.all)), ylim = ylims, cex = 0, xlab = 'Beat Index', ylab = expression(paste(italic(h[t]), ' (nats / beat)')))
}else{
  if (data.type == 'nb'){
    ylims = c(-1, 3) # negentropy rate for nb files
  }else{
    ylims = c(-1, 4) # negentropy rate for mixed + rossler
  }
  
  plot(0, 0, xlim = c(1, length(Xt.all)), ylim = ylims, cex = 0, xlab = 'Beat Index', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
}

for (segment.ind in 1:length(segments)){
  lines(segments[[segment.ind]][1]:segments[[segment.ind]][2], er.by.segments[[segment.ind]])
  
  # segment.er = median(er.by.segments[[segment.ind]], na.rm = TRUE)
  segment.er = mean(er.by.segments[[segment.ind]], na.rm = TRUE)
  segment.sd = sd.by.segments[[segment.ind]]
  segment.apen = apen.by.segments[segment.ind]
  
  if (entropy.type == 'er'){
    segment.apen.unnormalized = segment.apen - log(0.2*segment.sd)
  }else{
    # segment.apen.unnormalized = -1*(segment.apen - 0.5*log(2*pi*exp(1)*var(Xt))) - log(0.2*segment.sd)
    segment.apen.unnormalized = 0.5*log(2*pi*exp(1)*segment.sd^2) - segment.apen - log(0.2*segment.sd)
  }
  
  
  segments(x0 = segments[[segment.ind]][1], x1 = segments[[segment.ind]][2], y0 = segment.er, col = 'red', lwd = 2)
  
  text(x = segments[[segment.ind]][1] + (segments[[segment.ind]][2] - segments[[segment.ind]][1])/2, y = ylims[1]+0.5, labels = sprintf(fmt = '%.2f', segment.er), col = 'red', cex = 1.5)
  text(x = segments[[segment.ind]][1] + (segments[[segment.ind]][2] - segments[[segment.ind]][1])/2, y = ylims[1]+0.25, labels = sprintf(fmt = '%.2f (%.2f)', segment.apen, segment.apen.unnormalized), col = 'darkgreen', cex = 1.5)
  text(x = segments[[segment.ind]][1] + (segments[[segment.ind]][2] - segments[[segment.ind]][1])/2, y = ylims[1], labels = sprintf(segment.sd, fmt = '%.3f'), col = 'blue', cex = 1.5)
}

lines(uniformize(Xt.all) + ylims[2]-1, type = 'l', lty = 1, lwd = 0.5)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# Segmenting estimation of both variance and entropy rate!

# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, max(beat.times.all)), main = sprintf('Smoothing = %g', smooth.mean.var))
# for (segment.ind in 1:length(segments)){
#   lines(beat.times.by.segments[[segment.ind]], 0.5*(logvar.by.segments[[segment.ind]] + log(2*pi*exp(1))), col = 'darkgreen')
#   lines(beat.times.by.segments[[segment.ind]], er.by.segments[[segment.ind]], col = 'red')
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 2, col = 'red', lwd = 1)
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 2, col = 'red', lwd = 1, lty = 2)
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 2, col = 'red', lwd = 1, lty = 2)  
# }
# lines(beat.times.all, uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)

fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-er-by-time-segmented-discontinuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, max(beat.times.all)), main = sprintf('Smoothing = %g, Segment Padding = %g', smooth.mean.var, front.back.cutoff), xlab = 'Time (s)', ylab = expression(paste(italic(h)[italic(t)], ' (nats / beat)')))
for (segment.ind in 1:length(segments)){
  lines(beat.times.by.segments[[segment.ind]], 0.5*(logvar.by.segments[[segment.ind]] + log(2*pi*exp(1))), col = 'darkgreen', lwd = 3)
  lines(beat.times.by.segments[[segment.ind]], smooth.er.by.segments[[segment.ind]], col = 'red', lwd = 3)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 1, col = 'red', lwd = 1)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)  
  abline(v = beat.times.all[segments[[segment.ind]][1]], col = 'red', lty = 3, lwd = 3)
}
lines(beat.times.all, 2*uniformize(Xt.all, a = 0, b = 1.4)+1, lwd = 0.5)
legend('bottomright', legend = c('Pointwise Gaussian Entropy', 'Entropy Rate'), col = c('darkgreen', 'red'), lwd = 3, cex = 2)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)), main = sprintf('Smoothing = %g', smooth.mean.var))
# for (segment.ind in 1:length(segments)){
#   lines(segments[[segment.ind]][1]:segments[[segment.ind]][2], 0.5*(logvar.by.segments[[segment.ind]] + log(2*pi*exp(1))) - er.by.segments[[segment.ind]])
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 2, col = 'red', lwd = 1)
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 2, col = 'red', lwd = 1, lty = 2)
#   lines(beat.times.by.segments[[segment.ind]], uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 2, col = 'red', lwd = 1, lty = 2)  
# }
# lines(uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)



fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-ner-by-time-segmented-discontinuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-1, 3), xlim = c(1, max(beat.times.all)), main = sprintf('Smoothing = %g, Segment Padding = %g', smooth.mean.var, front.back.cutoff), xlab = 'Time (s)', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
for (segment.ind in 1:length(segments)){
  lines(beat.times.by.segments[[segment.ind]], 0.5*(logvar.by.segments[[segment.ind]] + log(2*pi*exp(1))) - smooth.er.by.segments[[segment.ind]], lwd = 3)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 1, col = 'red', lwd = 1)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)  
  abline(v = beat.times.all[segments[[segment.ind]][1]], col = 'red', lty = 3, lwd = 3)
}
lines(beat.times.all, 2*uniformize(Xt.all, a = 0, b = 1.4)+1, lwd = 0.5)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# Not segmenting estimation of variance!

# plot.new()
# 
# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, max(beat.times.all)))
# for (segment.ind in 1:length(segments)){
#   lines(beat.times.by.segments[[segment.ind]], er.by.segments[[segment.ind]], col = 'red')
# }
# lines(beat.times.all, 0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen')
# lines(beat.times.all, uniformize(Xt.all)+2, lwd = 0.5)
# 
# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)))
# for (segment.ind in 1:length(segments)){
#   lines(segments[[segment.ind]][1]:segments[[segment.ind]][2], smooth.er.by.segments[[segment.ind]], col = 'red')
# }
# lines(0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen')
# lines(uniformize(Xt.all)+2, lwd = 0.5)
# 
# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)))
# for (segment.ind in 1:length(segments)){
#   lines(segments[[segment.ind]][1]:segments[[segment.ind]][2], 0.5*(temporal.logvar.all[segments[[segment.ind]][1]:segments[[segment.ind]][2]] + log(2*pi*exp(1))) - er.by.segments[[segment.ind]])
# }
# lines(uniformize(Xt.all)+2, lwd = 0.5)
# 
# plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)))
# for (segment.ind in 1:length(segments)){
#   lines(segments[[segment.ind]][1]:segments[[segment.ind]][2], 0.5*(temporal.logvar.all[segments[[segment.ind]][1]:segments[[segment.ind]][2]] + log(2*pi*exp(1))) - smooth.er.by.segments[[segment.ind]])
# }
# lines(uniformize(Xt.all)+2, lwd = 0.5)
# 
# # Plot by-beat, not segmenting estimation of variance or entropy rate!
# 
plot.new()

smooth.er.all = fitted(npreg(bws = c(smooth.mean.var), txdat = beat.times.for.er, tydat = er.all, exdat = beat.times.all, ckertype = ckertype))

par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)), xlab = 'Beat Index', ylab = expression(paste(italic(h)[italic(t)], ' (nats / beat)')))
lines(er.all, col = 'red')
lines(0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen')
lines(uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)
for(segment.ind in 1:length(segments)){
  abline(v = segments[[segment.ind]], col = 'red', lwd = 1, lty = 3)
}
legend('bottomright', legend = c('Gaussian Entropy', 'Entropy Rate'), col = c('darkgreen', 'red'), lwd = 3, cex = 2)

fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-er-by-beat-segmented-continuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, length(Xt.all)), xlab = 'Beat Index', ylab = expression(paste(italic(h)[italic(t)], ' (nats / beat)')))
lines(smooth.er.all, col = 'red', lwd = 3)
lines(0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen', lwd = 3)
for(segment.ind in 1:length(segments)){
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 1, col = 'red', lwd = 1)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)  
  abline(v = segments[[segment.ind]], col = 'red', lwd = 1, lty = 3)
}
legend('bottomright', legend = c('Gaussian Entropy', 'Entropy Rate'), col = c('darkgreen', 'red'), lwd = 3, cex = 2)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# plot(0, cex = 0, ylim = c(0, 3), xlim = c(1, length(Xt.all)), xlab = 'Beat Index', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
# lines(0.5*(temporal.logvar.all + log(2*pi*exp(1))) - er.all)
# lines(uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)
# for(segment.ind in 1:length(segments)){
#   abline(v = segments[[segment.ind]], col = 'red', lwd = 1, lty = 3)
# }

fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-ner-by-beat-segmented-continuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(0, 3), xlim = c(1, length(Xt.all)), xlab = 'Beat Index', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
lines(0.5*(temporal.logvar.all + log(2*pi*exp(1))) - smooth.er.all, lwd = 3)
lines(uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)
for(segment.ind in 1:length(segments)){
  abline(v = segments[[segment.ind]], col = 'red', lwd = 1, lty = 3)
}
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# Plot by-time, not segmenting estimation of variance or entropy rate!

plot.new()

par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, max(beat.times.all)), xlab = 'Time (s)', ylab = expression(paste(italic(h)[italic(t)], ' (nats / beat)')))
lines(beat.times.for.er, er.all, col = 'red')
lines(beat.times.all, 0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen')
lines(beat.times.all, uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)
for(segment.ind in 1:length(segments)){
  abline(v = beat.times.all[segments[[segment.ind]]], col = 'red', lwd = 1, lty = 3)
}
legend('bottomright', legend = c('Gaussian Entropy', 'Entropy Rate'), col = c('darkgreen', 'red'), lwd = 3, cex = 2)

fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-er-by-time-segmented-continuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-4, 3), xlim = c(1, max(beat.times.all)), xlab = 'Time (s)', ylab = expression(paste(italic(h)[italic(t)], ' (nats / beat)')))
lines(beat.times.all, smooth.er.all, col = 'red', lwd = 3)
lines(beat.times.all, 0.5*(temporal.logvar.all + log(2*pi*exp(1))), col = 'darkgreen', lwd = 3)
for(segment.ind in 1:length(segments)){
  abline(v = beat.times.all[segments[[segment.ind]]], col = 'red', lwd = 1, lty = 3)
}
lines(beat.times.all, 2*uniformize(Xt.all, a = 0, b = 1.4)+1, lwd = 0.5)
legend('bottomright', legend = c('Gaussian Entropy', 'Entropy Rate'), col = c('darkgreen', 'red'), lwd = 3, cex = 2)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

# plot(0, cex = 0, ylim = c(0, 3), xlim = c(1, max(beat.times.all)), xlab = 'Time (s)', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
# lines(beat.times.all, 0.5*(temporal.logvar.all + log(2*pi*exp(1))) - er.all)
# lines(beat.times.all, uniformize(Xt.all, a = 0, b = 1.4)+2, lwd = 0.5)
# for(segment.ind in 1:length(segments)){
#   abline(v = beat.times.all[segments[[segment.ind]]], col = 'red', lwd = 1, lty = 3)
# }

fig.fname = sprintf('/Users/complexity/Dropbox/papers/hrv-ActEn/figures/%s%s-pCV-ner-by-time-segmented-continuous-smooth%g.pdf', data.type, figure.suffix, smooth.mean.var)
pdf(fig.fname, family = "CMU Serif", width = 30)
par(mar=c(5,5,2,1), cex.lab = 2, cex.axis = 2)
plot(0, cex = 0, ylim = c(-1, 3), xlim = c(1, max(beat.times.all)), xlab = 'Time (s)', ylab = expression(paste(negent[italic(t)], ' (nats / beat)')))
lines(beat.times.all, 0.5*(temporal.logvar.all + log(2*pi*exp(1))) - smooth.er.all, lwd = 3)
for(segment.ind in 1:length(segments)){
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]], a = 0, b = 1.4) + 1, col = 'red', lwd = 1)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] + 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)
  lines(beat.times.by.segments[[segment.ind]], 2*uniformize(mean.by.segments[[segment.ind]] - 2*sqrt(exp(logvar.by.segments[[segment.ind]])), a = 0, b = 1.4) + 1, col = 'red', lwd = 1, lty = 2)  
  abline(v = beat.times.all[segments[[segment.ind]]], col = 'red', lwd = 1, lty = 3)
}
lines(beat.times.all, 2*uniformize(Xt.all, a = 0, b = 1.4)+1, lwd = 0.5)
dev.off()
embed_fonts(fig.fname, outfile = fig.fname)

fname.save = sprintf('/Users/complexity/Documents/R/hrv-analysis/smooth-er-saved/%s%s-pCV-ner-by-time-segmented-smooth%g.dat', data.type, figure.suffix, smooth.mean.var)

write.csv(x = data.frame(beat.times.all = beat.times.all, smooth.er.all = smooth.er.all, smooth.ner.all = 0.5*(temporal.logvar.all + log(2*pi*exp(1))) - smooth.er.all), file = fname.save, quote = FALSE, row.names = FALSE)

}
}

# library(fractal)
# 
# DFA.out = DFA(x = rnorm(100000), scale.ratio = 2, scale.min = 100)
# plot(DFA.out)

plot(uniformize(entropies.np), type = 'l', col = 'red')
lines(uniformize(Xt), col = 'blue')
