##################################################
## Code to create Dushoff Flu A vs B vs weather plot for CI lecture
##################################################

## Fake data
xx <- rnorm(100,50, 10) 
yy <- rnorm(100,35, 7)
zz <- rnorm(100,10, 20)
dat <- cbind(xx,yy,zz)
xs <- seq(.5,2.5,1)

fgcol <- 'black' ## foreground color
for(ff in c(F,T)) {
  pdf(paste0('Flu A vs B vs Weather',' CIs'[ff], '.pdf'), w = 5, h = 5)
  par('ps' = 14, mgp = c(4,1,0), mar = c(4,6,2,2), ## point size, axis label locations , margins
      fg = fgcol, col.axis=fgcol, col.lab=fgcol) ## set foreground & other colors
  cols <- rainbow(3)
  plot(0,0,type = 'n', xlim = c(0,3), ylim = c( -100, 100),
       xlab = '', ylab = 'Attributable Deaths (per 10,000)', xaxt = 'n', bty = 'n',
       main = '')
  abline(h = 0, lwd = 2)
  points(xs, apply(dat,2,mean), pch = 19, col = cols)
  axis(1, at = xs, c('Flu A', 'Flu B', 'Weather'))
  if(ff) { ## add CIs?
    for(ii in 1:ncol(dat)) { ## for each CI
      arrows(xs[ii], quantile(dat[,ii],.025), xs[ii], quantile(dat[,ii],.975),  lwd = 2, col = cols[ii],
             angle = 90, length = .1, code = 3)
    }
  }
  dev.off() 
}
