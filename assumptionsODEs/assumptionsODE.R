## (Hidden) Assumptions of Simple ODE Models
## Clinic on Meaningful Modeling of Epidemiological Data
## International Clinics on Infectious Disease Dynamics and Data Program
## African Institute for Mathematical Sciences, Muizenberg, South Africa
##
## Juliet R.C. Pulliam, 2014
## Some Rights Reserved
## CC BY-NC 3.0 (http://creativecommons.org/licenses/by-nc/3.0/)

dir <- "."
source("../pctVax.R")

pcol <- "#696464"

pdf('pvaccNone.pdf', w = 8, h=5)
# Percent needed to vaccinate plot
par(mfcol=c(1,1),mar=c(9,9,2,1),bty="n",lwd=2, 'ps' = 30, mgp = c(7,2,0),
		fg=pcol, col.axis=pcol, col.lab=pcol)
curve(1-1/x,xlab=expression(R[0]),ylab=expression(p[v]),ylim=c(0,1),xlim=c(0,20),lwd=6, las = 1)
dev.off()
