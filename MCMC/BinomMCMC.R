library(car); library(parallel)
## ##################################################
## ## Plotting parameters
## fgcol <- 'white'
## brdcol <- 'black' ## barplot bar borders
## hcol <- gray(.5) ## hypothetical distribution colors
## xcol <- 'red' ## extreme value in tail color
## pcol <- 'purple' ## color of probability densities or exact values
## ycol <- 'yellow' ## color of emphasized text
## ps <- 24
## cx <- 1.2 ## ps/12 
## ymax <- .15

par('las'=1, bty = 'n', 'ps'=12) ## axes labels always written horizontal

size <- 100
truePrev <- .3
## Sample from this distribution once (use 28 for powerpoint slides).
sampPos <- 28                          #rbinom(1,100,.3)
sampPrev <- sampPos/size

gaussianProposal <- function(params, sd = .1)
    rnorm(length(params), mean = params, sd = sd)

unifProposal <- function(params, halfwidth = .1)
    runif(length(params), mean = params, min = params - halfwidth, max = params + halfwidth)
logPrior <- function(prevalence) ## uniform prior from 0 to 1 on prevalence
    dunif(prevalence, min = 0, max = 1, log = T)
logLikelihood <- function(prevalence, data = list(size=size, sampPos=sampPos))
    dbinom(data$sampPos, size = data$size, prob = prevalence, log = T)
logPosterior <- function(prevalence, data = list(size=size, sampPos=sampPos))
    logPrior(prevalence) + logLikelihood(prevalence, data)
posterior <- function(prevalence, data = list(size=size, sampPos=sampPos))
    exp(logPosterior(prevalence, data))

par(mfrow = c(2,2))
curve(logPrior, 0, 1)
curve(logLikelihood, 0, 1)
curve(logPosterior, 0,1)
curve(posterior, 0, 1)

## Sample on a logit-probability scale
runMCMC <- function(iterations, startvalue = runif(1, logit(.01), logit(.99)), 
                          proposer = gaussianProposal, verbose = 0){
    if(verbose > 0) browser()
    chain <- array(dim = c(iterations+1, length(startvalue)))
    chain[1,] <- startvalue
    for(ii in 1:iterations){
        proposal <- proposer(chain[ii,])
        MHratio <- exp(logPosterior(inv.logit(proposal)) - logPosterior(inv.logit(chain[ii,])))
        if(runif(1) < MHratio){
            chain[ii+1,] <- proposal
        }else{
            chain[ii+1,] <- chain[ii,]
        }
    }
    return(inv.logit(chain))
}

runMCMC(10)

nchains <- 4
its <- 1000
plot(0, type = 'n', xlab = 'iteration', ylab = 'prevalence', col = rainbow(nchains)[1], ylim = c(0,1), xlim = c(0, its))
for(ii in 1:nchains) lines(runMCMC(1000), col = rainbow(nchains)[ii])

defParList <- function() list(freq = F, xlab = '', ylab = '', xaxt = 'n',
                              breaks = seq(0, 1 , by = .01), ylim = c(0,100),
                              col = 'black', main = '') 

mcmcHist <- function(chains, parList=defParList(), proposer = gaussianProposal) {

propCol <- gray(.4)
par(mar = c(7,8,1,1))
    parList$x <- chains
    do.call(hist, parList)


    x0 <- chains[nrow(chains),]
    xseq <- seq(0,1, l = 1000)
    yseq <- dnorm(xseq, x0, sd = .05)
    xseqP <- c(xseq, rev(xseq))
    yseqP <- c(-yseq * 12/max(yseq), rep(0, length(yseq)))
par(xpd=NA)
    polygon(xseqP, yseqP, col = propCol, border=NA)
axis(1, at = seq(0,1, by = .1), line = 3)
mtext('prevalence', 1, 5)
mtext('proposal\ndistribution', 2, 2, srt=90, at = min(yseq)*.8, col = propCol)
mtext('MCMC sample\nfrequency', 2, 2, srt=90, at = parList$ylim[2]/2)

dbinom(
    
browser()
3
}

mcmcHist(runMCMC(100))

## ## Sample on a probability scale
## runMCMC <- function(iterations, startvalue = runif(1), proposer = gaussianProposal, logitTr=T, verbose = 0){
##     if(verbose > 0) browser()
##     chain <- array(dim = c(iterations+1, length(startvalue)))
##     chain[1,] <- startvalue
##     for(ii in 1:iterations){
##         proposal <- proposer(chain[ii,])
##         MHratio <- exp(logPosterior(proposal) - logPosterior(chain[ii,]))
##         if(runif(1) < MHratio){
##             chain[ii+1,] <- proposal
##         }else{
##             chain[ii+1,] <- chain[ii,]
##         }
##     }
##     return(chain)
## }
## logPosterior <- function(prevalence, data = list(size=size, sampPos=sampPos)) {
##     if(logPrior(prevalence[1]==-Inf)) {
##         return(-Inf) }else{
##             logPrior(prevalence) + logLikelihood(prevalence, data)
##         }
## }

