load('data/condstat.RData')

# N = 1000 # Full time series: takes longer
N = 200 # Shorter time series: runs faster

half.window.length = 50

choose.p.out = choose.p.lwo(x = condstat[1:N], p.max = 4, p.min = 1, half.window.length = half.window.length)

plot(choose.p.out$ps, choose.p.out$er.by.p, pch = 16, cex = 2, xlab = expression(paste('Autoregressive Order ', italic(p))), ylab = substitute(paste(CV[l](italic(p), bold(k))), list(l = half.window.length)))
abline(v = choose.p.out$p.best)

cat(sprintf('LWOCV chose p = %g.\n\n', choose.p.out$p.best))

future.vals = seq(-20, 20, by = 0.1)

par(mfrow = c(4, 1))
for (time.point in 250:273){
  specific.past = condstat[1:(length(choose.p.out$pred.dens$xbw)) + time.point]

  fhat.out = evaluate.fhat(future.vals = future.vals, specific.past = specific.past, pred.dens = choose.p.out$pred.dens)

  plot(fhat.out$eval.points[, 1], fhat.out$fhat, type = 'l', ylim = c(0, 1), xlab = 'Future Value', ylab = 'Predictive Density')

  lines(c(-5, 5), y = uniformize(tail(specific.past, n = 2), a = -20, b = 20))
  points(c(-5, 5), y = uniformize(tail(specific.past, n = 2), a = -20, b = 20), pch = 16, cex = 4)
  abline(h = 0.5)
}

cat(sprintf('Estimating the SPecific ENtropy RAte (spenra)...\n\n'))

estimate.spenra.out = estimate.spenra(condstat[1:N], choose.p.out$pred.dens, half.window.length = half.window.length, integral.lowerbound = -Inf, integral.upperbound = Inf)

par(mar=c(5,7,2,1), cex.lab = 2, cex.axis = 2, mfrow = c(1, 1))
plot(estimate.spenra.out, type = 'l', ylim = c(1.7, 2.7), xlab = 'Time (au)', ylab = expression(paste(widehat(italic(h))[italic(t)], ' (nats / symbol)')))
abline(h = c(1.744021, 2.517551), col = 'blue', lty = 2, lwd = 3)