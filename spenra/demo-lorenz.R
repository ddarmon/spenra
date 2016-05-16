load('data/lorenz.RData')

# N = 1000 # Full time series: takes longer
N = 200 # Shorter time series: runs faster

half.window.length = 50

choose.p.out = choose.p.lwo(x = lorenz.ibi[1:N], p.max = 8, p.min = 1, half.window.length = half.window.length)

par(mfrow = c(1, 1))
plot(choose.p.out$ps, choose.p.out$er.by.p, pch = 16, cex = 2, xlab = expression(paste('Autoregressive Order ', italic(p))), ylab = substitute(paste(CV[l](italic(p), bold(k))), list(l = half.window.length)))
abline(v = choose.p.out$p.best)

cat(sprintf('LWOCV chose p = %g.\n\n', choose.p.out$p.best))

future.vals = seq(0, 5, by = 0.01)

par(mfrow = c(4, 1))
for (time.point in 250:273){
  specific.past = lorenz.ibi[1:(length(choose.p.out$pred.dens$xbw)) + time.point]

  fhat.out = evaluate.fhat(future.vals = future.vals, specific.past = specific.past, pred.dens = choose.p.out$pred.dens)

  plot(fhat.out$eval.points[, 1], fhat.out$fhat, type = 'l', ylim = c(0, 4), xlab = 'Future Value', ylab = 'Predictive Density')
}

cat(sprintf('Estimating the SPecific ENtropy RAte (spenra)...\n\n'))

estimate.spenra.out = estimate.spenra(lorenz.ibi[1:N], choose.p.out$pred.dens, half.window.length = half.window.length, integral.lowerbound = -Inf, integral.upperbound = Inf)

par(mar=c(5,7,2,1), cex.lab = 2, cex.axis = 2, mfrow = c(1, 1))
plot(estimate.spenra.out, type = 'l', xlab = 'Time (Interval Index)', ylab = expression(paste(widehat(italic(h))[italic(t)], ' (nats / event)')), ylim = c(-1, 0.5))
