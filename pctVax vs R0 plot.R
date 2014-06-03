## pctVax.R
## JRCP 25 March 2014
##
## Unit 6 - Applied Ecology of Infectious Diseases
## ZOO4926/ZOO6927 Spring 2014, Department of Biology
## University of Florida, Gainesville, Florida, USA

# Disease list
# R0 estimates from Anderson & May 1982 Science:
dl <- list(measles=c(13.7,12.5,18),
	pertussis=c(14.3,17.1,12.2),
	chickenpox = c(8.5,9),
	diptheria=c(6.6,6.4),
	mumps=c(4.3,7.1),
	rubella=c(6,6.7),
	poliomyelitis=c(5.9,6.2))

## dir <- '~/Dropbox/BiomathsMedPH/Lectures/Code/lectureCode/'
setwd(dir)

pdf('pvacc.pdf', w = 8, h=5)
# Percent needed to vaccinate plot
if(!exists("pcol")){
	pcol <- 'white'
}
par(mfcol=c(1,1),mar=c(9,9,2,1),bty="n",lwd=2, 'ps' = 30, mgp = c(7,2,0),
    fg=pcol, col.axis=pcol, col.lab=pcol)
curve(1-1/x,xlab=expression(R[0]),ylab=expression(p[v]),ylim=c(0,1),xlim=c(0,20),lwd=6, las = 1)
sapply(1:length(dl),function(dis){
	xx <- dl[[dis]]
	yy <- 1-1/dl[[dis]]
	points(xx,yy,pch=16,cex=2.5,col=dis+1)
})
par('ps'=16)
legend('bottomright', names(dl),col=1:length(dl)+1, pch = 19, bty = 'n', cex = 1.5)
dev.off()
