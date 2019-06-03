library(boot); library(parallel); library(animation); library(coda)
## setwd('~/Documents/R Repos/lectureCode/MCMC/')
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

size <- 100
truePrev <- .3
## Sample from this distribution once (use 28 for powerpoint slides).
sampPos <- 28                          #rbinom(1,100,.3)
sampPrev <- sampPos/size

gaussianProposal <- function(sd = .5, params = 1)
    list(rprop = function(params) rnorm(1, mean = params, sd = sd),
         dprop = function(x) dnorm(x, mean = params, sd = sd), sd = sd, type='norm')
unifProposal <- function(halfwidth = .1, params = 1)
    list(rprop = function(params) runif(1, min = params - halfwidth, max = params + halfwidth),
         dprop = function(x) dunif(x, min = params - halfwidth, max = params + halfwidth), halfwidth = halfwidth, type='unif')
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
dev.off()
## par(mfrow = c(2,2), mar = c(5,6,1,0))
## curve(Prior(x, logBetaPrior), 0, 1)
## #curve(Prior(x, logUnifPrior), 0, 1)
## curve(Likelihood, 0, 1)
## curve(logLikePrior, 0,1)
## curve(LikePrior, 0, 1)

## Sample on a logit-probability scale
runMCMC <- function(iterations, startvalue = runif(1, logit(.01), logit(.99)),
                    plotter = NULL, plotNM = NULL, plotDIR = 'stills', logPriorFxn = logUnifPrior,
                    proposer = gaussianProposal,verbose = 0, ...){
    if(verbose > 0) browser()
    chain <- array(dim = c(iterations+1, length(startvalue)))
    chain[1,] <- startvalue
    for(ii in 1:iterations){
        ##        browser()
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
                    verbose = verbose, logPriorFxn=logPriorFxn, ii = ii, ...)
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
                     ps = 25, propDistCol = 'yellow', propCol='brown', curCol = 'dodger blue', lmar = 25, logPriorFxn=logUnifPrior,
                     noProps = F, ii=1, iiShow = NULL, burn = 0) {
    if(is.null(iiShow) | ii %in% iiShow) {
            if(verbose > 0) browser()
        chains <- chains[!is.na(chains[,1]),,drop=F]
        marLine <- 14
        numPlot <- ifelse(ii==1 | !is.null(iiShow), 3, 2)
        for(bb in 1:numPlot) { ## show 1st frame without proposal, then 2nd with proposal
            if(!is.null(plotNM)) png(paste0('stills/',plotNM, '/', plotNM, '-',
                                            formatC(ii, 4, flag='0'),'-',bb,'.png'), width = 800*resScl, height = 600*resScl)
            layout(matrix(1:4,4,1), h = c(1.5,.8,.8,1))
            par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
                cex.axis=1.2, cex.lab=1.2, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,0))
            par(mar = c(13,lmar,1,1), 'ps'=ps)
            parList <- within(parList, {
                x <- chains[(burn+1):(nrow(chains)-1),]
                ##if(bb==1) x <- chains[1:(nrow(chains)-1),] else x <- chains
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
            if('norm'==proposer$type) {
                yseq <- dnorm(logit(xseq), logit(x0), sd = proposer$sd)
                propDens <- dnorm(logit(proposal), logit(x0), proposer$sd)
            }else{
                yseq <- dunif(logit(xseq), min = logit(x0) - proposer$halfwidth, max = logit(x0) + proposer$halfwidth)
                propDens <- dunif(logit(proposal), min = logit(x0) - proposer$halfwidth, max = logit(x0) + proposer$halfwidth)
            }
            xseqP <- c(xseq, rev(xseq))
            scl <- 18/max(yseq)*ceiling(max(hh$counts)/40)
            dep <- 5
            yseqP <- c(-yseq * scl, rep(0, length(yseq)))
            ## Proposal distribution
            if(!noProps | bb<3) {
                par(xpd=NA)
                polygon(xseqP, yseqP-dep, col = makeTransparent(propDistCol, alpha = 100), border=NA)
                if(verbose>0) browser()
                accepted <- chains[nrow(chains),]==proposal
                if(accepted) ltyProp <- 1 else ltyProp <- 3
                segments(x0, -dep, x0, min(yseqP)-dep, col = curCol, lwd = lwd)
                if(bb>1) segments(proposal, -dep, proposal, -dep-scl*propDens, 
                                  col = propCol, lty = ltyProp, lwd = lwd)
            }
            axis(1, at = seq(0,1, by = .1), pos = -dep, labels=F)
            mtext('proposal\ndistribution', 2, marLine, srt=90, at = -(par('usr')[4] + par('usr')[3])/2 , col = propDistCol, adj = .5)
            mtext('MCMC sample\ndistribution', 2, marLine, srt=90, at = parList$ylim[2]/2, col=lgreen, adj = .5)
            ## Likelihood X Prior
            par(mar = c(3,lmar,1,1), col.axis=mainCol)
            curve(LikePrior(x, logPriorFxn = logPriorFxn), 0, 1, ylab='', xlab='', xaxt='n', lwd = 3)
            if(!noProps | bb<3) {
                segments(x0, 0, x0, LikePrior(x0, logPriorFxn = logPriorFxn), col = curCol, lwd = lwd)
                if(bb>1) segments(proposal, 0, proposal, LikePrior(proposal, logPriorFxn = logPriorFxn), col = propCol, lty = ltyProp, lwd = lwd)
            }
            axis(1, at = seq(0,1, by = .1), labels=F)
            mtext('likelihood\nx\nprior', 2, marLine, srt=90, adj = .5)
            if(!noProps | bb <3)
                legend('right', leg =c('current', 'proposed (accepted)', 'proposed (rejected)'), lwd = lwd, col = c(curCol, propCol, propCol),
                       lty = c(1,1,3), bty = 'n', cex = 1.6, seg.len = 3)
            ## Likelihood
            curve(Likelihood, 0, 1, ylab='', xlab='', xaxt='n', col = lpurp, lwd = 3)
            axis(1, at = seq(0,1, by = .1), labels=F)
            mtext('likelihood', 2, marLine, adj = .5, srt=90, col = lpurp)
            ## Prior
            par(mar = c(8,lmar,1,1))
            curve(Prior(x, logPriorFxn = logPriorFxn), 0, 1, ylab='', lwd = 3, 
                  xlab='', xaxt='n', col.lab=mainCol,col.axis=mainCol, col=lorange)
            axis(1, at = seq(0,1, by = .1), mgp = c(3,3,0))
            mtext('prior', 2, marLine, adj = .5, srt=90, col=lorange) 
            mtext('prevalence', 1, 6, adj = .5)
            if(!is.null(plotNM)) graphics.off()
        }
    }
}    


mcmcHistTrace <- function(chains, parList=defParList(), proposer = gaussianProposal, proposal = NA, verbose = 0, lwd = 5, plotNM=NULL, 
                          ps = 25, propDistCol = 'yellow', propCol='brown', curCol = 'dodger blue', lmar = 25, logPriorFxn=logUnifPrior,
                          noProps = F, ii=1, iiShow = NULL, burn = 0) {
    if(is.null(iiShow) | ii %in% iiShow) {
        chains <- chains[!is.na(chains[,1]),,drop=F]
        marLine <- 14
        numPlot <- ifelse(ii==1, 3, 2)
        for(bb in 1:numPlot) { ## show 1st frame without proposal, then 2nd with proposal
            if(!is.null(plotNM)) png(paste0('stills/',plotNM, '/', plotNM, '-',
                                            formatC(ii, 4, flag='0'),'-',bb,'.png'), width = 800*resScl, height = 600*resScl)
            layout(matrix(1:3,3,1), h = c(1.5,1,1))
            par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
                cex.axis=1.2, cex.lab=1.2, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,0))
            par(mar = c(13,lmar,1,1), 'ps'=ps)
            parList <- within(parList, {
                x <- chains[(burn+1):(nrow(chains)-1),]
                ##if(bb==1) x <- chains[1:(nrow(chains)-1),] else x <- chains
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
            if('norm'==proposer$type) {
                yseq <- dnorm(logit(xseq), logit(x0), sd = proposer$sd)
                propDens <- dnorm(logit(proposal), logit(x0), proposer$sd)
            }else{
                yseq <- dunif(logit(xseq), min = logit(x0) - proposer$halfwidth, max = logit(x0) + proposer$halfwidth)
                propDens <- dunif(logit(proposal), min = logit(x0) - proposer$halfwidth, max = logit(x0) + proposer$halfwidth)
            }
            xseqP <- c(xseq, rev(xseq))
            scl <- 18/max(yseq)*ceiling(max(hh$counts)/40)
            dep <- 5
            yseqP <- c(-yseq * scl, rep(0, length(yseq)))
            ## Proposal distribution
            if(!noProps | bb<3) {
                par(xpd=NA)
                polygon(xseqP, yseqP-dep, col = makeTransparent(propDistCol, alpha = 100), border=NA)
                if(verbose>0) browser()
                accepted <- chains[nrow(chains),]==proposal
                if(accepted) ltyProp <- 1 else ltyProp <- 3
                segments(x0, -dep, x0, min(yseqP)-dep, col = curCol, lwd = lwd)
                if(bb>1) segments(proposal, -dep, proposal, -dep-scl*propDens, 
                                  col = propCol, lty = ltyProp, lwd = lwd)
            }
            axis(1, at = seq(0,1, by = .1), pos = -dep, labels=F)
            mtext('proposal\ndistribution', 2, marLine, srt=90, at = -(par('usr')[4] + par('usr')[3])/2 , col = propDistCol, adj = .5)
            mtext('MCMC sample\ndistribution', 2, marLine, srt=90, at = parList$ylim[2]/2, col=lgreen, adj = .5)
            ## Likelihood X Prior
            par(mar = c(8,lmar,1,1), col.axis=mainCol)
            curve(LikePrior(x, logPriorFxn = logPriorFxn), 0, 1, ylab='', xlab='', xaxt='n', lwd = 3)
            if(!noProps | bb<3) {
                segments(x0, 0, x0, LikePrior(x0, logPriorFxn = logPriorFxn), col = curCol, lwd = lwd)
                if(bb>1) segments(proposal, 0, proposal, LikePrior(proposal, logPriorFxn = logPriorFxn), 
                                  col = propCol, lty = ltyProp, lwd = lwd)
            }
            axis(1, at = seq(0,1, by = .1), mgp = c(3,3,0))
            mtext('likelihood\nx\nprior', 2, marLine, srt=90, adj = .5)
            if(!noProps | bb <3)
                legend('right', leg =c('current', 'proposed (accepted)', 'proposed (rejected)'), lwd = lwd, 
                       col = c(curCol, propCol, propCol),
                       lty = c(1,1,3), bty = 'n', cex = 1.6, seg.len = 3)
            mtext('prevalence', 1, 6, adj = .5)
            ## Trace
            xmax <- ceiling(nrow(chains)/100)*100
            plot(0, type = 'n', xlab = '', ylab = '', col = rainbow(1)[1], ylim = c(0,1), xlim = c(0, xmax), mgp = c(3,3,0))
            lines(chains, col = lgreen)
            mtext('prevalence\nMCMC chain', 2, marLine, adj = .5, srt=90)
            mtext('iteration', 1, 6, adj = .5)
            if(!is.null(plotNM)) graphics.off()
        }
    }
}    


## Trace
tracePlot <- function(chainList, upTo = 10^4, plotNM=NULL, marLine = 12, lmar=20, ps = 25, bb = 1, parList=defParList(), alpha = 150) {
    if(!is.null(plotNM)) if(!file.exists(file.path('stills',plotNM))) dir.create(file.path('stills',plotNM))
    if(!is.null(plotNM)) png(paste0('stills/',plotNM, '/', plotNM, '-', formatC(ii, 4, flag='0'),'-',bb,'.png'), width = 800*resScl, height = 600*resScl)
    layout(matrix(1:2,2,1), h = c(1.5,1))
    par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
        cex.axis=1.2, cex.lab=1.2, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,0))
    par(mar = c(8,lmar,1,1), 'ps'=ps)
    nchains <- length(chainList)
    ## subset chains
    for(ii in 1:nchains) chainList[[ii]] <- chainList[[ii]][1:min(upTo, nrow(chainList[[ii]])), ]
    mostIts <- max(unlist(lapply(chainList, length)))
    xmax <- ceiling(mostIts/100)*100
    ## Histograms (shaded)
    histLS <- list()
    for(ii in 1:nchains) histLS[[ii]] <- hist(chainList[[ii]], seq(0, 1 , by = .025), plot = F)
    maxCount <- max(unlist(lapply(histLS, function(x) x$counts)))
    ylim <- c(0, ceiling(maxCount/40)*40)    
    plot(0, type = 'n', xlab = '', ylab = '', ylim = ylim, xlim = c(0, 1), mgp = c(3,3,0), xaxt='n')
            axis(1, at = seq(0,1, by = .1), mgp = c(3,3,0))
    mtext('MCMC sample\ndistribution', 2, marLine, srt=90, at = 0, col=lgreen, adj = .5)
    histDens <- function(histObj, col) with(histObj, {
        xs <- rep(breaks,each=2)
        xs <- xs[-c(1,length(xs))]
        xs <- c(xs, rev(xs))
        ys <- rep(counts, each = 2)
        ys <- c(ys, rep(0, length(ys)))
        polygon(xs, ys, col = col, border = NA)
    })
    for(ii in 1:nchains) histDens(histLS[[ii]], col = makeTransparent(rainbow(nchains)[ii], alpha = alpha))
    mtext('prevalence', 1, 6, adj = .5)
    ## Trace Plots
    par(mar = c(8,lmar,1,1), col.axis=mainCol)
    plot(0, type = 'n', xlab = '', ylab = '', col = rainbow(1)[1], ylim = c(0,1), xlim = c(0, xmax), mgp = c(3,3,0))
    for(ii in 1:nchains) lines(chainList[[ii]], col = rainbow(nchains)[ii])
    mtext('prevalence\nMCMC chain', 2, marLine, adj = .5, srt=90)
    mtext('iteration', 1, 6, adj = .5)
    ## Gelman-Rubin
    mcCH <- chainList
    for(ii in 1:nchains) mcCH[[ii]] <- mcmc(mcCH[[ii]])
    mcCH <- as.mcmc(mcCH)
    gdag <- gelman.diag(mcCH)
    mtext(paste0('Gelman-Rubin Diagnostic = ', signif(gdag$psrf[,1],4)), 3, -3)
    if(!is.null(plotNM)) graphics.off()
}

## ## sample on a probability scale
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

## show posterior at the end (after calculation norm constant) to show that we got it right


