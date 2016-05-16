load('data/condstat.RData')

choose.p.out = choose.p.lwo(x = condstat[1:1000], p.max = 5, p.min = 1, half.window.length = 50)

par(mfrow = c(4, 1))
# for (time.point in 0:24){
# for (time.point in 175:200){
for (time.point in 250:273){
  specific.past = condstat[1:(length(choose.p.out$pred.dens$xbw)) + time.point]

  fhat.out = evaluate.fhat(future.vals = future.vals, specific.past = specific.past, pred.dens = choose.p.out$pred.dens)

  plot(fhat.out$eval.points[, 1], fhat.out$fhat, type = 'l', ylim = c(0, 1))

  lines(c(-5, 5), y = uniformize(specific.past, a = -20, b = 20))
  points(c(-5, 5), y = uniformize(specific.past, a = -20, b = 20), pch = 16, cex = 4)
  abline(h = 0.5)
}

estimate.spenra.out = estimate.spenra(condstat, choose.p.out$pred.dens, half.window.length = 50, integral.lowerbound = -Inf, integral.upperbound = Inf)

par(mfrow = c(1, 1))
plot(estimate.spenra.out, type = 'l', ylim = c(1.8, 2.7))
