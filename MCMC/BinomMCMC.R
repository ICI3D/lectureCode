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


mainCol <- 'white'
backCol <- 'black'
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
logUnifPrior <- function(prevalence) ## uniformative prior from 0 to 1 on prevalence
        dunif(prevalence, min = 0, max = 1, log = T) 
logBetaPrior <- function(prevalence, shape1 = 8, shape2 = 40) ## informative beta prior
        dbeta(prevalence, shape1=shape1, shape2=shape2, log = T)
logLikelihood <- function(prevalence, data = list(size=size, sampPos=sampPos))
    dbinom(data$sampPos, size = data$size, prob = prevalence, log = T)
logLikePrior <- function(prevalence, data = list(size=size, sampPos=sampPos), logPriorFxn = logUnifPrior, ...)
    logPriorFxn(prevalence, ...) + logLikelihood(prevalence, data)
Prior <- function(x, logPriorFxn=logUnifPrior, ...) exp(logPriorFxn(x, ...))
Likelihood <- function(x) exp(logLikelihood(x))
LikePrior <- function(x, logPriorFxn=logUnifPrior, ...) exp(logLikePrior(x, logPriorFxn = logPriorFxn, ...))

opar <- par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,0))
## par(mfrow = c(2,2), mar = c(5,6,1,0))
## curve(Prior(x, logBetaPrior), 0, 1)
## #curve(Prior(x, logUnifPrior), 0, 1)
## curve(Likelihood, 0, 1)
## curve(logLikePrior, 0,1)
## curve(LikePrior, 0, 1)

## Sample on a logit-probability scale
runMCMC <- function(iterations, startvalue = runif(1, logit(.01), logit(.99)),
                    plotter = NULL, plotNM = NULL, plotDIR = 'stills', logPriorFxn = logUnifPrior,
                    proposer = gaussianProposal,verbose = 0){
    if(verbose > 0) browser()
    chain <- array(dim = c(iterations+1, length(startvalue)))
    chain[1,] <- startvalue
    for(ii in 1:iterations){
        ##browser()
        proposal <- proposer$rprop(chain[ii,])
        MHratio <- exp(logLikePrior(inv.logit(proposal), logPriorFxn = logPriorFxn) - 
                       logLikePrior(inv.logit(chain[ii,]), logPriorFxn = logPriorFxn))
        if(runif(1) < MHratio){
            chain[ii+1,] <- proposal
        }else{
            chain[ii+1,] <- chain[ii,]
        }
        ## PLOTTING
        if(!is.null(plotter)) {
            if(!is.null(plotNM)) if(!file.exists(file.path('stills',plotNM))) dir.create(file.path('stills',plotNM))
            plotter(inv.logit(chain), proposal = inv.logit(proposal), proposer = proposer, plotNM=plotNM,
                    verbose = 0, logPriorFxn=logPriorFxn, ii = ii)
        }
    }
    return(inv.logit(chain))
}

lorange <- rgb(246,190,146, maxColorValue = 255)
lpurp <- rgb(178,160,198, maxColorValue = 255)
lblue <- rgb(185,221,231, maxColorValue = 255)
lgreen <- 'light green'#rgb(214,228,189, maxColorValue = 255)
defParList <- function() list(freq = T, xlab = '', ylab = '', xaxt = 'n', border = NA,
                              breaks = seq(0, 1 , by = .01), ylim = c(0,40),
                              col = lgreen, main = '') 

mcmcHist <- function(chains, parList=defParList(), proposer = gaussianProposal, proposal = NA, verbose = 0, lwd = 5, plotNM=NULL, 
                     propDistCol = 'yellow', propCol='brown', curCol = 'dodger blue', lmar = 19, logPriorFxn=logUnifPrior, ii=1) {
    chains <- chains[!is.na(chains[,1]),,drop=F]
    marLine <- 12
    for(bb in 1:2) { ## show 1st frame without proposal, then 2nd with proposal
        if(!is.null(plotNM)) png(paste0('stills/',plotNM, '/', plotNM, '-',
                                        formatC(ii, 4, flag='0'),'-',bb,'.png'), width = 800, height = 600)
        layout(matrix(1:4,4,1), h = c(1.5,1,1,1))
        par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,0))
            par(mar = c(9,lmar,1,1), 'ps'=18)
            parList <- within(parList, {
                if(bb==1) x <- chains[1:(nrow(chains)-1),] else x <- chains
                plot <- F
            })
            ## MCMC histogram
            hh <- hist(chains, parList$breaks, plot = F)
            parList <- within(parList, {
                ylim <- c(0, ceiling(max(hh$counts)/40)*40)
                plot <- T
            })
            hh <- do.call(hist, parList)
            x0 <- chains[nrow(chains)-1,]
            xseq <- seq(.01,.99, l = 1000)
            yseq <- dnorm(logit(xseq), logit(x0), sd = proposer$sd)
            xseqP <- c(xseq, rev(xseq))
            scl <- 25/max(yseq)*ceiling(max(hh$counts)/40)
            dep <- 5
            yseqP <- c(-yseq * scl, rep(0, length(yseq)))
            ## Proposal distribution
            par(xpd=NA)
            polygon(xseqP, yseqP-dep, col = makeTransparent(propDistCol, alpha = 100), border=NA)
            if(verbose>0) browser()
            accepted <- chains[nrow(chains),]==proposal
            if(accepted) ltyProp <- 1 else ltyProp <- 3
            segments(x0, -dep, x0, min(yseqP)-dep, col = curCol, lwd = lwd)
            if(bb>1) segments(proposal, -dep, proposal, -dep-scl*dnorm(logit(proposal), logit(x0), proposer$sd), 
                              col = propCol, lty = ltyProp, lwd = lwd)
            axis(1, at = seq(0,1, by = .1), pos = -dep, labels=F)
            mtext('proposal\ndistribution', 2, marLine, srt=90, at = -(par('usr')[4] + par('usr')[3])/5 , col = propDistCol, adj = .5)
            mtext('MCMC sample\ndistribution', 2, marLine, srt=90, at = parList$ylim[2]/2, col=lgreen, adj = .5)
            ## Likelihood X Prior
            par(mar = c(3,lmar,1,1), col.axis=mainCol)
            curve(LikePrior(x, logPriorFxn = logPriorFxn), 0, 1, ylab='', xlab='', xaxt='n', lwd = 3)
            segments(x0, 0, x0, LikePrior(x0, logPriorFxn = logPriorFxn), col = curCol, lwd = lwd)
            if(bb>1) segments(proposal, 0, proposal, LikePrior(proposal, logPriorFxn = logPriorFxn), col = propCol, lty = ltyProp, lwd = lwd)
            axis(1, at = seq(0,1, by = .1), labels=F)
            mtext('likelihood\nx\nprior', 2, marLine, srt=90, adj = .5)
            legend('right', leg =c('current', 'proposed (accepted)', 'proposed (rejected)'), lwd = lwd, col = c(curCol, propCol, propCol),
                   lty = c(1,1,3), bty = 'n', cex = 1.6, seg.len = 3)
            ## Likelihood
            curve(Likelihood, 0, 1, ylab='', xlab='', xaxt='n', col = lpurp, lwd = 3)
            axis(1, at = seq(0,1, by = .1), labels=F)
            mtext('likelihood', 2, marLine, adj = .5, srt=90, col = lpurp)
            ## Prior
            par(mar = c(6,lmar,1,1))
            curve(Prior(x, logPriorFxn = logPriorFxn), 0, 1, ylab='', lwd = 3, 
                  xlab='prevalence', xaxt='n', col.lab=mainCol,col.axis=mainCol, col=lorange)
            axis(1, at = seq(0,1, by = .1))
            mtext('prior', 2, marLine, adj = .5, srt=90, col=lorange) 
            if(!is.null(plotNM)) graphics.off()
        }
    }    

## Run 3 chains, save the first 10 frames of each
nits <- 3000
sdvec <- c(.05, .2, 2)
ff <- 'logBetaPrior'
ii <- 2
for(ff in c('logUnifPrior','logBetaPrior')) {
    for(ii in 1:length(sdvec)) {
        set.seed(4)
        nm <- paste0('movies/',ff,'-',sdvec[ii], '.mov')
        if(file.exists(nm)) file.remove(nm)
        saveVideo({
            ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
            runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), logPriorFxn = get(ff))
        }, video.name = nm, other.opts = "-b 1000k -pix_fmt yuv420p", ani.width = 800, ani.height = 600)
        set.seed(5)
        test <- runMCMC(7, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), startvalue = logit(.25),
                        plotNM = paste0(ff,'-',sdvec[ii]))
    }
}

runMCMC(1, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=.5))

runMCMC(10, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=.5))#, plotNM = 'test.5-')

nchains <- 4
its <- 1000
par(mfrow=c(1,1))
plot(0, type = 'n', xlab = 'iteration', ylab = 'prevalence', col = rainbow(nchains)[1], ylim = c(0,1), xlim = c(0, its))
for(ii in 1:nchains) lines(runMCMC(1000), col = rainbow(nchains)[ii])


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
