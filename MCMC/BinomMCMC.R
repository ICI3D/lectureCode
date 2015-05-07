library(boot); library(parallel); library(animation)
setwd('~/Documents/R Repos/lectureCode/MCMC/')
source('utilityFxns.R')
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


opar <- par('las'=1, bty = 'n', 'ps'=12) ## axes labels always written horizontal

size <- 100
truePrev <- .3
## Sample from this distribution once (use 28 for powerpoint slides).
sampPos <- 28                          #rbinom(1,100,.3)
sampPrev <- sampPos/size

gaussianProposal <- function(sd = .5, params = 1)
    list(rprop = function(params) rnorm(1, mean = params, sd = sd),
         dprop = function(x) dnorm(x, mean = params, sd = sd), sd = sd)
unifProposal <- function(halfwidth = .1, params = 1)
    list(rprop = function(params) runif(1, min = params - halfwidth, max = params + halfwidth),
         dprop = function(x) dunif(x, min = params - halfwidth, max = params + halfwidth), halfwidth = halfwidth)
logPrior <- function(prevalence) ## uniform prior from 0 to 1 on prevalence
    dunif(prevalence, min = 0, max = 1, log = T)
logLikelihood <- function(prevalence, data = list(size=size, sampPos=sampPos))
    dbinom(data$sampPos, size = data$size, prob = prevalence, log = T)
logLikePrior <- function(prevalence, data = list(size=size, sampPos=sampPos))
    logPrior(prevalence) + logLikelihood(prevalence, data)
## posterior <- function(prevalence, data = list(size=size, sampPos=sampPos))
##     exp(logPosterior(prevalence, data))

par(mfrow = c(2,2))
curve(logPrior, 0, 1)
curve(logLikelihood, 0, 1)
curve(logLikePrior, 0,1)
curve(posterior, 0, 1)

## Sample on a logit-probability scale
runMCMC <- function(iterations, startvalue = runif(1, logit(.01), logit(.99)),
                    plotter = NULL, plotNM = NULL, plotDIR = 'Videos',
                    proposer = gaussianProposal,verbose = 0){
    if(verbose > 0) browser()
    chain <- array(dim = c(iterations+1, length(startvalue)))
    chain[1,] <- startvalue
    for(ii in 1:iterations){
        ##browser()
        proposal <- proposer$rprop(chain[ii,])
        MHratio <- exp(logLikePrior(inv.logit(proposal)) - logLikePrior(inv.logit(chain[ii,])))
        if(runif(1) < MHratio){
            chain[ii+1,] <- proposal
        }else{
            chain[ii+1,] <- chain[ii,]
        }
        ## PLOTTING
        if(!is.null(plotter)) {
            if(!is.null(plotNM)) {
                if(!file.exists(file.path('Videos',plotNM))) dir.create(file.path('Videos',plotNM))
                jpeg(paste0('Videos/',plotNM, '/', plotNM, formatC(ii, 4, flag='0'), '.jpeg'), quality = 300)
            }
            par(opar)
            plotter(inv.logit(chain), proposal = inv.logit(proposal), proposer = proposer, verbose = 0)#ii > 10)
            if(!is.null(plotNM)) graphics.off()
        }
    }
    return(inv.logit(chain))
}

runMCMC(2000, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=.5), plotNM = 'test.5-')


defParList <- function() list(freq = T, xlab = '', ylab = '', xaxt = 'n',
                              breaks = seq(0, 1 , by = .01), ylim = c(0,100),
                              col = 'black', main = '') 

mcmcHist <- function(chains, parList=defParList(), proposer = gaussianProposal, proposal = NA, verbose = 0) {
    propCol <- 'red'
    chains <- chains[!is.na(chains[,1]),,drop=F]
    propColTr <- makeTransparent(propCol,50)
    par(mar = c(7,8,1,1))
    layout(matrix(1:2,2,1), h = c(2,1))
    parList <- within(parList, {
        x <- chains
        plot <- F
    })
    hh <- hist(chains, parList$breaks, plot = F)
    parList <- within(parList, {
        ylim <- c(0, ceiling(max(hh$counts)/100)*100)
        plot <- T
    })
    hh <- do.call(hist, parList)
    x0 <- chains[nrow(chains)-1,]
    xseq <- seq(.01,.99, l = 1000)
    yseq <- dnorm(logit(xseq), logit(x0), sd = proposer$sd) ##proposer$ddist(xseq, x0, sd = .1)
    xseqP <- c(xseq, rev(xseq))
    scl <- 50/max(yseq)
    dep <- 5
    yseqP <- c(-yseq * scl, rep(0, length(yseq)))
    par(xpd=NA)
    polygon(xseqP, yseqP-dep, col = propColTr, border=NA)
    if(verbose>0) browser()
    accepted <- chains[nrow(chains),]==proposal
    segments(x0, -dep, x0, min(yseqP)-dep, col = 'blue')
    segments(proposal, -dep, proposal, -dep-scl*dnorm(logit(proposal), logit(x0), proposer$sd), col = 'brown', lty = 2-accepted)
    ##axis(1, at = seq(0,1, by = .1), line = 3)
#    axis(1, at = seq(0,1, by = .1), pos = -dep, labels=F)
    mtext('proposal\ndistribution', 2, 2, srt=90, at = -(par('usr')[4] + par('usr')[3])/5 , col = propCol)
    mtext('MCMC sample\ndistribution', 2, 2, srt=90, at = parList$ylim[2]/2)
    par(mar = c(5,8,1,1))
    curve(posterior, 0, 1, ylab='', xlab='prevalence')
    segments(x0, 0, x0, posterior(x0), col = 'blue')
    segments(proposal, 0, proposal, posterior(proposal), col = 'brown', lty = 2-accepted)    
    axis(1, at = seq(0,1, by = .1))
    mtext('posterior', 2, 3, srt=90)#, at = parList$ylim[2]/2)
}


nchains <- 4
its <- 1000
par(mfrow=c(1,1))
plot(0, type = 'n', xlab = 'iteration', ylab = 'prevalence', col = rainbow(nchains)[1], ylim = c(0,1), xlim = c(0, its))
for(ii in 1:nchains) lines(runMCMC(1000), col = rainbow(nchains)[ii])

runMCMC(5, plotter=mcmcHist, verbose = 0)#, plotNM = 'test')

mcmcHist(runMCMC(100))

## ## Sample on a probability scale
## runMCMC <- function(iterations, startvalue = runif(1), proposer = gaussianProposal, logitTr=T, verbose = 0){
##     if(verbose > 0) browser()
##     chain <- array(dim = c(iterations+1, length(startvalue)))
##     chain[1,] <- startvalue
##     for(ii in 1:iterations){
##         proposal <- proposer(chain[ii,])
##         MHratio <- exp(logLikePrior(proposal) - logLikePrior(chain[ii,]))
##         if(runif(1) < MHratio){
##             chain[ii+1,] <- proposal
##         }else{
##             chain[ii+1,] <- chain[ii,]
##         }
##     }
##     return(chain)
## }
## logLikePrior <- function(prevalence, data = list(size=size, sampPos=sampPos)) {
##     if(logPrior(prevalence[1]==-Inf)) {
##         return(-Inf) }else{
##             logPrior(prevalence) + logLikelihood(prevalence, data)
##         }
## }

## start with just likelihood+ prior on bottom, then add sampling scheme above it
## show posterior at the end (after calculation norm constant) to show that we got it right

## show w/ informative prior
