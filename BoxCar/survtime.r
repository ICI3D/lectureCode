##################################################
## Set up various plots to show exponential assumptions and why we
## need boxcar models to get gamma-distributed waiting times.

## Graphical parameters
fgcol <- 'black' ## foreground color
hcol <- 'red' ## color of data histograms
ps <- 14 ## pointsize
xlim <- c(0,25)
ylim <- c(0,.3)
xlab <- "years since infection"
width <- 3.5 
height <- 4
bty <- "n"

##################################################
## calculate gamma hazards by using the CDF to get instantaneous hazards
step <- .2 ## to evaluate pdf
seqs <- seq(0,25,by = step)
gamsurv <- 100*(1-pgamma(seqs,4,1/2.5))
gamhaz <- (gamsurv[2:length(gamsurv)]-gamsurv[1:(length(gamsurv)-1)]) / (-step) / gamsurv[1:(length(gamsurv)-1)]

##################################################
## Hazards
pdf("haz (exponential).pdf", width = width, height = height)
par(mar=c(5,5,1,1), fg = fgcol, col.lab = fgcol, col.axis = fgcol, 'ps'=ps)
mpar <- par() ## save par for future plots
plot(seqs, rep(1/10, length(seqs)),xlab = xlab, ylab = "hazard of mortality", main="", type = "l", lwd = 2, xlim = xlim, ylim = ylim, cex.lab = cex.lab, bty=bty, col=fgcol)
dev.off() 
pdf("haz (gamma).pdf", width = width, height = height)
par(mar=mpar)
plot(seqs[-1],gamhaz, type = "l", xlab = xlab, ylab = "hazard of mortality", main="", lwd = 2, xlim = xlim, ylim = ylim, cex.lab = cex.lab, bty=bty, col=fgcol)
dev.off()

##################################################
## PDF of survival times w/ & w/o data
fdat <- rgamma(100, 4,1/2.5) ## fake gamma distributed data
breaks <- -1:500
ylab <- "probability of surv time"
xlab <- "surv time (years)"
ylim = c(0,.15)
for(dd in c(F,T)) {
  pdf(paste0("exp PDF",' data'[dd], ".pdf"), width = width, height = height)
  par(mar=mpar)
  curve(dexp(x,1/10), col = fgcol, lwd = 2, xlab = xlab,
        xlim = xlim, ylim = ylim, ylab = ylab, main = "", #add = dd,
        yaxt = "n", cex.lab = cex.lab, bty=bty)
    hist(fdat, freq=F, col = hcol, plot=dd, breaks = breaks, add = T)
    curve(dexp(x,1/10), add = T)
  axis(2, at = c(0,.05,.1,.15), labels = c("0",".05",".1",".15"),cex.lab = cex.lab)
  dev.off()
  pdf(paste0("gamma PDF",' data'[dd], ".pdf"), width = width, height = height)
  par(mar=mpar)
  curve(dgamma(x,4,1/2.5), col = fgcol, lwd = 2, xlab = xlab, xlim = xlim,
        ylim = ylim, ylab = ylab, main = "",yaxt = "n", cex.lab = cex.lab, bty=bty)
  hist(fdat, freq=F, col = hcol, plot=dd, breaks = breaks, add = T)
  curve(dgamma(x,4,1/2.5), add = T)
  axis(2, at = c(0,.05,.1,.15), labels = c("0",".05",".1",".15"))
  dev.off()
}

##################################################
## CDFs
ylab <- "probability of dying >t"
ylim <- c(0,1)
pdf("exp cdf.pdf", width = width, height = height)
par(mar=mpar)
curve(pexp(x,1/10), col = fgcol, add = F, lwd = 2, xlab = xlab, xlim = xlim, ylim = ylim,
      ylab = ylab, main = "", yaxt = "n", cex.lab = cex.lab, bty=bty)
axis(2, at = c(0,.5,1), labels = c("0","0.5","1"), cex.lab = cex.lab)
dev.off()
pdf("gamcdf.pdf", width = width, height = height)
par(mar=mpar)
curve(pgamma(x,4,1/2.5), col = fgcol, add = F, lwd = 2, xlab = xlab, xlim = xlim,
      ylim = ylim, ylab = ylab, main = "",yaxt = "n", cex.lab = cex.lab, bty=bty)
axis(2, at = c(0,.5,1), labels = c("0","0.5","1"), cex.lab = cex.lab)
dev.off()

##################################################
## 1 -CDFS
ylim <- c(0,1)
ylab <- "cum. prob. of survival to >t"
pdf("expSurv.pdf", width = width, height = height)
par(mar=mpar)
curve(1 - pexp(x,1/10), col = fgcol, add = F, lwd = 2, xlab = xlab, xlim = xlim, ylim = ylim,
      ylab = ylab, main = "", yaxt="n", cex.lab = cex.lab, bty=bty)
axis(2, at = c(0,.5,1), labels = c("0","0.5","1"), cex.lab = cex.lab)
dev.off()
pdf("gamSurv.pdf", width = width, height = height)
par(mar=mpar)
curve(1-pgamma(x,4,1/2.5), col = fgcol, add = F, lwd = 2, xlab = xlab, xlim = xlim,
      ylim = ylim, ylab = ylab, main = "", yaxt="n", cex.lab = cex.lab, bty=bty)
axis(2, at = c(0,.5,1), labels = c("0","0.5","1"), cex.lab = cex.lab)
dev.off()
