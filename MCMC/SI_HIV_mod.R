library(deSolve); library(ggplot2); library(MASS); library(sfsmisc); library(mnormt); library(ellipse); library(emdbook)
setwd('~/Documents/R Repos/lectureCode/MCMC/')
source('utilityFxns.R')

opar <- par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,2))
lorange <- rgb(246,190,146, maxColorValue = 255)
lpurp <- rgb(178,160,198, maxColorValue = 255)
lblue <- rgb(185,221,231, maxColorValue = 255)

disease_params <- function(Beta = 0.3, alpha = 1, progRt = 1/2.5,
                           birthRt = .03, deathRt = 1/60, ...)
    return(as.list(environment()))

disease_params()

init <- c(S=.99, I1=.0001, I2=0, I3=0, I4=0, CI = 0, CD = 0)
tseq <- seq(1970, 2015, by = 1/12)
Is <- paste0('I',1:4)

## SI ODE model
SImod <- function(tt, yy, parms) with(parms, {
    ## state variables
    S <- yy[1]  ## Susceptibles
    I1 <- yy[2] ## HIV stage 1
    I2 <- yy[3] ## HIV stage 2
    I3 <- yy[4] ## HIV stage 3
    I4 <- yy[5] ## HIV stage 4
    CI <- yy[7] ## cumulative incidence
    CD <- yy[8] ## cumulative mortality
    ## derived quantitties
    I <- I1+I2+I3+I4           ## total infecteds
    N <- I + S                 ## total population
    mort <- progRt * I4 / N ## HIV-related mortality
    FOI <- Beta * exp(-alpha * I/N) ## Force of infection
    ## state variable derivatives (ODE system)
    deriv <- rep(NA,5)  
    deriv[1]<-	birthRt*N - deathRt*S - FOI*S*I/N
    deriv[2]<-	FOI*S*I/N - progRt*I1 - deathRt*I1
    deriv[3]<-	progRt*I1 - progRt*I2 - deathRt*I2
    deriv[4]<-	progRt*I2 - progRt*I3 - deathRt*I3
    deriv[5]<-	progRt*I3 - progRt*I4 - deathRt*I4
    deriv[6]<-	FOI*S*I/N
    deriv[7]<-	progRt*I4
    return(list(deriv))
})

sampleEpidemic <- function(simDat, tseq = seq(1978, 2010, by = 2), numSamp = rep(80, length(tseq)), verbose=0) {
    if(verbose>0) browser()
    simDat$I <- rowSums(simDat[, Is])
    prev_at_sample_times <- simDat[simDat$time %in% tseq, 'I']
    numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
    lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = numSamp)
    uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = numSamp)    
    return(data.frame(time = tseq, numPos, numSamp, sampPrev =  numPos/numSamp,
                      lci = lci, uci = uci))
}

set.seed(4)
trueParms <- disease_params(Beta = .9, alpha = 8, progRt = 1/2.5)
simDat <- as.data.frame(lsoda(init, tseq, SImod, parms=trueParms))
simDat$I <- rowSums(simDat[, Is])
simDat$N <- rowSums(simDat[, c('S',Is)])
plot(simDat$time, simDat$I, xlab = '', ylab = 'prevalence', type = 'l', ylim = c(0,.4), col='red')

obsDat <- sampleEpidemic(simDat, verbose = 0)
points(obsDat$time, obsDat$sampPrev, col = 'red', pch = 16, cex = 2)
arrows(obsDat$time, obsDat$uci, obsDat$time, obsDat$lci, col = makeTransparent('red'), len = .025, angle = 90, code = 3)

## Log-Likelihood
llikelihood <- function(parms = disease_params(), obsDat, verbose = 0) {
    simDat <- as.data.frame(lsoda(init, tseq, SImod, parms=parms))
    simDat$I <- rowSums(simDat[, Is])
    if(verbose > 5) browser()
    lls <- dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$I[simDat$time %in% obsDat$time], log = T)
    return(sum(lls))
}
llikelihood(obsDat=obsDat)

## Log-Prior (assume uninformative)
lprior <- function(parms=disease_params()) with(parms, {
    lp <- 0
    return(lp)
})

llikePrior <- function(fit.params=NULL, ## parameters to fit
                       ref.params = disease_params(), ## reference parameters
                       obsDat, verbose = 0) { ## observed data
    parms <- within(ref.params, { ## subs fitting parameters into reference parameter vector
        for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
        rm(nm)
    })
    llikelihood(parms, obsDat=obsDat, verbose = verbose) + lprior(parms)
}
llikePrior(obsDat=obsDat)

logParms <- function(fit.params) {
    fit.params <- log(fit.params)
    names(fit.params) <- paste0('log',names(fit.params))
    return(fit.params)
}
unlogParms <- function(fit.params) {
    fit.params <- exp(fit.params)
    names(fit.params) <- sub('log','', names(fit.params))
    return(fit.params)
}
unlogParms(logParms(c(alpha = 3, Beta=.3)))

initBounds <- data.frame(rbind( ## for initial conditions
                               c(log(.2),log(2)) ## beta
                               ,c(log(1), log(30)) ## alpha
                               ,c(log(1),log(1/10)))) ## progRt
colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('logBeta','logalpha','logprogRt')
class(initBounds[,2]) <- class(initBounds[,1]) <- 'numeric'
initBounds

initRand <- function(fit.params) {
    fit.params <- logParms(fit.params)
    tempnm <- names(fit.params)
    for(nm in tempnm) fit.params[nm] <- runif(1, min = initBounds[rownames(initBounds)==nm, 'lower'], 
                                              max =  initBounds[row.names(initBounds)==nm, 'upper'])
    return(unlogParms(fit.params))
}
initRand(c(alpha = 3, Beta = 1))


mcmcSampler <- function(current.params, ref.params=disease_params(), obsDat, seed = 1,
                        proposer = sequential.proposer(sdProps=sdProps),
                        adaptiveMCMC = F, startAdapt = 300,
                        plotterLoops = NULL, ## repeat plot of some iterations more to lengthen their display
                        showing = showingFXN(),
                        plotter = plotterTS, randInit = T, niter = 100, nburn = 0, adptBurn = 200,
                        verbose=0, plotNM=NULL, tell = 100) {
    if(verbose>2) browser()
    if(randInit) current.params <- initRand(current.params)
    nfitted <- length(current.params)
    vv <- 2 ## mcmc iteration
    accept <- 0
    curVal <- llikePrior(current.params, ref.params = ref.params, obsDat=obsDat, verbose = verbose)
    out <- matrix(NA, nr = niter, nc=length(current.params)+1)
    out[1,] <- c(current.params, nll = -curVal)
    colnames(out) <- c(names(current.params), 'nll')
    max_index <- dim(out)[1]
    last.it <- 0
    ## Store original covariance matrix
    if(proposer$type=='block') originalCovar <- get('covar', envir = environment(proposer$fxn)) 
    while(vv <= max_index) {
        if ((verbose > 1) || (verbose && (vv%%tell == 0))) print(paste("on iteration",vv,"of",last.it + niter + 1))
        ## Adaptive MCMC
        ## adapt covariance every 50 iterations
        if(adaptiveMCMC & proposer$type=='block' & vv > startAdapt & vv %% 50 == 0) {
            adptBurn <- min((startAdapt-50), adptBurn)
            save(list = ls(all.names = TRUE), file='dbgAdapt.Rdata')
            ## load(file='dbgAdapt.Rdata')
            adaptedCovar <- 2.38^2 / nfitted * cov.wt(log(out[adptBurn:(vv-1),1:nfitted]))$cov ## will converge for vv large
            adaptedCovar <- adaptedCovar*.95 + originalCovar*.05 ## 95% adapted & 5% original
            rownames(adaptedCovar) <- colnames(adaptedCovar) <- names(current.params)
            assign('covar', adaptedCovar, envir = environment(proposer$fxn))
        }
        proposal0 <- proposer$fxn(logParms(current.params))
        onpar <- proposal0$onpar
        propt <- proposal0$type
        proposal <- unlogParms(proposal0$proposal)
        propVal <- llikePrior(proposal, ref.params = ref.params, obsDat=obsDat, verbose = verbose)
        lmh <- propVal - curVal ## likelihood ratio = log likelihood difference
        if (is.na(lmh)) { ## if NA, print informative info but don't accept it
            print(list(lmh=lmh, proposal=exp(proposal), vv=vv, seed=seed))
        } else { ## if it's not NA then do acception/rejection algorithm
            if (verbose > 1) print( c(lmh=lmh, propVal=propVal) )
            ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
            if ( (lmh >= 0) | (runif(1,0,1) <= exp(lmh)) ) {
                current.params <- proposal
                if (vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
                curVal <- propVal
            }
        }
        out[vv, ] <- c(current.params, nll=curVal)
        vv <- vv+1
        aratio <- accept/((vv-nburn))
        if(!is.null(plotter)) {
                par(opar)
                plotter(out, vv, ref.params=ref.params, obsDat=obsDat, proposer = proposer, proposal = proposal,
                        onpar=onpar,proptype=propt, aratio = aratio, showing = showing, plotterLoops=plotterLoops)
        }
    }
    if(!is.null(plotNM)) graphics.off()
    colnames(out) <- c(names(current.params), 'nll')
    return(list(out = out[1:nrow(out)>(nburn+1),], aratio = aratio, current.params = current.params,
                ref.params=ref.params))
}


showingFXN <- function(shPost = T, shHPD = F, shPropKern=T, shHist = T,
                       shProp=T, shDat=T, shCurr=T, shAratio=T)
    return(as.list(environment()))

plotterParmDens <- function(out, vv, ref.params=disease_params(), plotNM=NULL, obsDat,
                            verbose=0, proposer = NULL, onpar = onpar,
                            showing = showingFXN(),
                            proptype=proptype, aratio = NULL,
                            proposal = NA, propDistCol = 'yellow', propCol='brown', curCol = 'dodger blue',
                            obsCol = gray(.7), plotterLoops = NULL,
                            every = 200, burn = 100,
                            marLine = 8, lmar=23, ps = 25, xlim = c(1,50), ylim = c(.2,2), log = 'xy',
                            bump = 5, nlevs = 40, ts.lwd = 7,
                            yparnm = expression(beta), xparnm=expression(alpha)) {
    with(showing, {
        out <- out[!is.na(out[,1]), colnames(out) !='nll', drop=F]
        newParms <- out[nrow(out) , colnames(out) !='nll']
        lastParms <- fit.params <- out[nrow(out)-1,]
        out <- out[1:(nrow(out)-1),,drop=F]
        parnms <- colnames(out)
        if(verbose > 7) browser()
        
        for(bb in 1:2) {
            ## Number of times to repeat each plot

            numLoop <- ifelse(is.null(plotterLoops), 1, plotterLoops$repeats[min(which(vv <  plotterLoops$breaks))])
            for(ll in 1:numLoop) {

                if(vv==90) save(list = ls(all.names = TRUE), file='dbg.Rdata')
                if(vv==210) save(list = ls(all.names = TRUE), file='dbgE210.Rdata')
                if(vv==211) save(list = ls(all.names = TRUE), file='dbgE211.Rdata')
                ##load(file='dbg.Rdata')
                if(!is.null(plotNM)) {
                    dirnm <- file.path('movies',plotNM)
                    if(!file.exists(dirnm)) dir.create(dirnm)
                    png(file.path('movies',plotNM,paste0('BivariateSeq-',vv,'-',bb,'.png')),
                        width = 700*resScl, height = 700*resScl)
                }
                ##          load(file='dbgE210.Rdata')
                ## png(paste0('movies/test.png'), width = 700*resScl, height = 700*resScl) #
                layout(matrix(c(3,1,4,5,2,4), 3, 2), w = c(1,.5), h = c(.6,1,.8))
                par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
                    cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,2), 'ps'=ps)
                par(mar=c(8,lmar,1,1))
                ## rescale axes every so many iterations
                if(vv + 1 > every & bb==1) {
                    axisNum <- floor(vv/every)
                    ##if(vv+1==every) browser()
                    tout <- out[burn:min(nrow(out), axisNum*every),]
                    xlim <- quantile(tout[,1], c(.02,.98)) * c(.7, 1/.7)
                    ylim <- quantile(tout[,2], c(.02,.98)) * c(.9, 1/.9)
                    nlevs <- round(nlevs*2.5)
                }
                plot(1,1, type = 'n', xlim = xlim, ylim = ylim, log = log, axes = F, xlab='',ylab='')
                if(grepl('x',log)) {
                    xticks <- axTicks(1, axp = c(xlim, 2), log = T)
                    if(length(xticks)<3) xticks <- axTicks(1, log = T)
                }else{
                    xticks <- pretty(xlim, 5)
                }
                if(grepl('y',log)) {
                    ##        yticks <- axTicks(2, axp = c(ylim, 2), log = T)
                    yticks <- axTicks(2, log = T)
                    ##        if(length(yticks)<3) yticks <- axTicks(1)
                }else{
                    yticks <- pretty(ylim, 5)
                }
                axis(2, at = yticks, mgp=c(4,2,0))
                axis(1, at = xticks, mgp=c(4,2,0))
                ## axis(1, at = ceiling(min(xlim)):floor(max(xlim)), lab = F)
                ## axis(2, at = ceiling(min(ylim)):floor(max(ylim)), lab = F)
                outOriginal <- out
                Lout <- out
                Lxlim <- xlim
                Lxticks <- xticks
                Lylim <- ylim
                Lyticks <- yticks
                Lxr <- xr <- range(c(out[,1], xlim))
                Lyr <- yr <- range(c(out[,2], ylim))
                if(grepl('x',log)) {
                    Lout[,1] <- log(out[,1])
                    Lxlim <- log(xlim)
                    Lxticks <- log(xticks)
                    Lxr <- log(xr)
                }
                if(grepl('y',log)) {
                    Lout[,2] <- log(out[,2])
                    Lylim <- log(ylim)
                    Lyticks <- log(yticks)        
                    Lyr <- log(yr)
                }
                if(shPost) {
                    if(nrow(out)>20) {
                        z <- kde2d(Lout[,1], Lout[,2], lims = c(Lxlim, Lylim), h = c(.5,.2))
                        if(grepl('x',log)) z$x <- exp(z$x)
                        if(grepl('y',log)) z$y <- exp(z$y)
                        ## Posterior bivariate density
                        cols2Ramp <- c('black','dark green', lgreen,'white')
                        cols <- apply(colorRamp(cols2Ramp)(seq(0,1, l = nlevs)), 1, function(x) rgb(x[1],x[2],x[3], max=255))
                        .filled.contour(z$x, z$y, z$z, levels = seq(0, max(z$z), l=nlevs), col = cols)
                    }
                    if(nrow(out)<30) points(outOriginal, col = makeTransparent('light green', alpha = 100), cex = 2, pch = 16)
                }
                mtext(xparnm, 1, marLine-3, cex = 1.5)
                mtext(yparnm, 2, marLine+3, cex = 1.5)
                
                if(proptype=='block') { ## Bivariate proposer
                    bivPDF <- function(x,y) dmnorm(cbind(x,y), mean = log(lastParms),
                                                   varcov = get('covar', envir = environment(proposer$fxn)))
                    xsAt <- seq(Lxlim[1], Lxlim[2], l=70)
                    ysAt <- seq(Lylim[1], Lylim[2], l=70)
                    bivDens <- outer(xsAt, ysAt , bivPDF)
                    if(grepl('x',log)) xsAt <- exp(xsAt)
                    if(grepl('y',log)) ysAt <- exp(ysAt)
                    ## Proposal bivariate density
                    nlevs <- 40
                    colsProp <- apply(colorRamp(c('white','yellow'),
                                                bias = 1)(seq(0,1, l = nlevs+1)), 1,
                                      function(x) rgb(x[1],x[2],x[3], max=255))
                    colsProp <- diag(makeTransparent(colsProp, c(0,0, seq(50,90, l = nlevs-2))))
                    if(shHPD) {
                        logHPD <- HPDregionplot(log(out[burn:nrow(out),]), prob = c(.95),
                                                lims = c(Lxr, Lyr), n = 45, h = c(.2,.15),
                                                add = T, col = NA, lwd = 5, plot = F) #, #n = 28, h = c(1,.08),
                        with(logHPD[[1]], lines(exp(x), exp(y), col = 'purple', lwd = 5))
                    }
                    if(shPropKern) 
                        .filled.contour(xsAt,ysAt, bivDens, levels = seq(0, max(bivDens), l=nlevs), col = colsProp) #
                    ##     lines(exp(ellipse(proposer$covar, centre = log(lastParms), level = .95)),
                    ## col = makeTransparent(propDistCol,150), lwd = 4)
                    ## polygon(exp(ellipse(proposer$covar, centre = log(lastParms), level = .95)), 
                    ## col = makeTransparent(propDistCol,50), border = NA)
                }

                accepted <- sum(newParms!=proposal)==0
                pchProp <- ifelse(accepted,19,21)
                if(shCurr) 
                    points(lastParms[1],lastParms[2], pch = 19, col = curCol, cex = 3)
                points( trueParms[names(lastParms[1])], trueParms[names(lastParms[2])], pch = 19,
                       col = 'white', cex = 3)
                if(shProp) 
                    if(bb>1) points(proposal[1],proposal[2], pch = pchProp, col = propCol, cex = 3, lwd = 4)

                if(shPropKern | shHist) {
                    ## Marginal Histograms
                    xbreaks <- seq(Lxr[1], Lxr[2], l=nlevs)
                    ybreaks <- seq(Lyr[1], Lyr[2], l=nlevs)
                    xhist  <-  hist(Lout[,1], plot=FALSE, breaks = xbreaks)
                    yhist  <-  hist(Lout[,2], plot=FALSE, breaks = ybreaks)
                    XinRange <- xhist$breaks >= Lxlim[1] & xhist$breaks <= Lxlim[2]
                    YinRange <- yhist$breaks >= Lylim[1] & yhist$breaks <= Lylim[2]
                    top  <-  max(c(xhist$counts, yhist$counts))
                    scl <- .1
                    dep <- .1

                    ## Y histogram
                    par(mar=c(8,bump+2,1,1))
                    plot(0,0, type='n',ylim = Lylim, xlim = c(0,top), axes=F, xlab='',ylab='',main='')
                    attr(yhist,'class') <- 'list'
                    yhist <- within(yhist, {
                        show <- breaks <= max(Lylim) & breaks  >= min(Lylim)
                        breaks <- breaks[show]
                        counts <- counts[show]
                    })
                    if(shHist)
                        with(yhist, rect(0, breaks[1:(length(breaks) - 1)], counts, breaks[2:length(breaks)],
                                         col = lgreen, border=NA))
                    axis(2, at = Lyticks, labels = F)
                    ## Y proposal
                    if(proptype=='sequential' & shPropKern) {
                        valSeq <- seq(Lylim[1],Lylim[2], l = 1000)
                        densSeq <- dnorm(valSeq, log(lastParms[2]), sd = proposer$sdProps[2])
                        propDens <- dnorm(log(proposal[2]), log(lastParms[2]), proposer$sdProps[2])
                        valSeqP <- c(valSeq, rev(valSeq))
                        densSeqP <- c(-densSeq * scl, rep(0, length(densSeq)))-dep
                        par(xpd=NA, new=T)
                        plot(0,0, type='n',ylim = Lylim, xlim = c(0,1), axes=F, xlab='',ylab='',main='')
                        polygon(densSeqP, valSeqP, col = makeTransparent(propDistCol, alpha = 100), border=NA)
                        accepted <- newParms[2]==proposal[2]
                        if(accepted) ltyProp <- 1 else ltyProp <- 3
                        if(shCurr)
                            segments(max(densSeqP), log(lastParms[2]), min(densSeqP), log(lastParms[2]),
                                     col = curCol, lwd = 3)
                        if(onpar==1 & shProp & bb==2)
                            segments(-dep, log(proposal[2]), -propDens*scl-dep, log(proposal[2]),
                                     col = propCol, lty = ltyProp, lwd = 3)
                        par(xpd=T)
                    }
                    ## X histogram
                    par(mar=c(bump,lmar,1,1))
                    plot(0,0, type='n',xlim = Lxlim, ylim = c(0,top), axes=F, xlab='',ylab='',main='')
                    attr(xhist,'class') <- 'list'
                    xhist <- within(xhist, {
                        show <- breaks <= max(Lxlim) & breaks  >= min(Lxlim)
                        breaks <- breaks[show]
                        counts <- counts[show]
                    })
                    if(shHist)
                        with(xhist, rect(breaks[1:(length(breaks) - 1)], 0, breaks[2:length(breaks)], counts,
                                         col = lgreen, border=NA))
                    axis(1, at = Lxticks, labels = F)
                    ## X proposer
                    if(proptype=='sequential' & shPropKern) {
                        valSeq <- seq(Lxlim[1],Lxlim[2], l = 1000)
                        densSeq <- dnorm(valSeq, log(lastParms[1]), sd = proposer$sdProps[1])
                        propDens <- dnorm(log(proposal[1]), log(lastParms[1]), proposer$sdProps[1])
                        valSeqP <- c(valSeq, rev(valSeq))
                        densSeqP <- c(-densSeq * scl, rep(0, length(densSeq)))-dep
                        par(xpd=NA, new=T)
                        plot(0,0, type='n',ylim = c(0,1), xlim = Lxlim, axes=F, xlab='',ylab='',main='')
                        polygon(valSeqP, densSeqP, col = makeTransparent(propDistCol, alpha = 100), border=NA)
                        accepted <- newParms[2]==proposal[2]
                        if(accepted) ltyProp <- 1 else ltyProp <- 3
                        if(shCurr)
                            segments(log(lastParms[1]), max(densSeqP), log(lastParms[1]), min(densSeqP),
                                     col = curCol, lwd = 3)
                        if(onpar==2 & shProp & bb==2)
                            segments(log(proposal[1]), -dep, log(proposal[1]), -propDens*scl-dep,
                                     col = propCol, lty = ltyProp, lwd = 3)
                        par(xpd=T)
                    }
                    if((vv+1) > (every - 30) & (vv+1) < every) mtext('about to\nzoom in', 2, 3)
                    if((vv+1) >= (every) & (vv+1) < (every+30)) mtext('just \nzoomed in', 2, 3)
                }else{
                    plot.new()
                    plot.new()
                }
                ## Time Series
                ## True parameters
                trueDat <- as.data.frame(lsoda(init, tseq, SImod, parms=trueParms))
                trueDat$I <- rowSums(trueDat[, Is])
                ## lastParms
                lastParmsAll <- within(ref.params, { ## subs fitting parameters into reference parameter vector
                    for(nm in parnms) assign(nm, as.numeric(lastParms[nm]))
                    rm(nm)
                })
                lastParmsDat <- as.data.frame(lsoda(init, tseq, SImod, parms=lastParmsAll))
                lastParmsDat$I <- rowSums(lastParmsDat[, Is])
                ## lastParms
                propParmsAll <- within(ref.params, { ## subs fitting parameters into reference parameter vector
                    for(nm in parnms) assign(nm, as.numeric(proposal[nm]))
                    rm(nm)
                })
                propParmsDat <- as.data.frame(lsoda(init, tseq, SImod, parms=propParmsAll))
                propParmsDat$I <- rowSums(propParmsDat[, Is])
                par(mar=c(6,lmar,.7,22))
                plot(trueDat$time, trueDat$I, xlab = '', ylab = '', type = 'l', ylim = c(0,.45),
                     col='white', lwd = ts.lwd, mgp = c(4,3,0))
                if(shCurr) 
                    lines(lastParmsDat$time, lastParmsDat$I, col = curCol, lwd = ts.lwd)
                accepted <- sum(newParms!=proposal)==0
                ltyProp <- ifelse(accepted,1,3)
                if(shProp) 
                    lines(propParmsDat$time, propParmsDat$I, col = propCol, lwd = ts.lwd, lty = ltyProp)
                ## add data
                if(shDat) {
                    points(obsDat$time, obsDat$sampPrev, col = obsCol, pch = 16, cex = 4)
                    arrows(obsDat$time, obsDat$uci, obsDat$time, obsDat$lci, col = makeTransparent(obsCol),
                           len = .025, angle = 90, code = 3, lwd = 3)
                }
                mtext('HIV\nprevalence', 2, marLine+ 7, adj = .5)
                par(xpd=NA)
                show <- c(T, shDat, shCurr, shProp, shProp)
                legend(2012, .5,
                       leg = c('truth', 'observed', 'current', 'proposal (accepted)', 'proposal (rejected)')[show],
                       lwd = c(4,0,4,4,4)[show], pch = c(NA, 16, 16, 16, 21)[show],
                       lty = c(1,NA,1,1,3)[show], seg.len = 4, pt.cex = 3,
                       col = c('white', obsCol, curCol, propCol, propCol)[show],
                       cex = 1.3, ncol = 1, bty = 'n', y.intersp = 1.5)
                par(xpd=F)
                par(mar=rep(0,4))
                if(shHist) 
                    plot(0,0, type = 'n', xlim = c(0,10), ylim = c(0,10), axes = F)
                if(shAratio) 
                    if(!is.null(aratio))    text(4.5,5, paste0('acceptance ratio = ', formatC(aratio,2, format='f')),
                                                 cex = 1.5)
                if(!is.null(plotNM)) dev.off()


                ##             for(ll in 1:(numLoop-1)) {
                ##  ##browser()
                ##             dev.copy(png)
                ## #            dev.off(dev.prev())
                ##             }
                ##             graphics.off()
            } ## ll loop
        } ## bb loop
    })
}

plotterTS <- function(fit.params=NULL, ref.params=disease_params(), plotNM=NULL, obsDat) {
    parms <- within(ref.params, { ## subs fitting parameters into reference parameter vector
        for(nm in names(fit.params)) assign(nm, as.numeric(fit.params[nm]))
        rm(nm)
    })
    ## simulate model
    simDat <- as.data.frame(lsoda(init, tseq, SImod, parms=parms))
    simDat$I <- rowSums(simDat[, Is])
    plot(simDat$time, simDat$I, xlab = '', ylab = 'prevalence', type = 'l', ylim = c(0,.4), col='red')
    ## add data
    points(obsDat$time, obsDat$sampPrev, col = 'red', pch = 16, cex = 2)
    arrows(obsDat$time, obsDat$uci, obsDat$time, obsDat$lci, col = makeTransparent('red'), len = .025, angle = 90, code = 3)
}

sequential.proposer <- function(sdProps) {
    nfitted <- length(sdProps)
    on <- 0
    return(list(sdProps = sdProps, type = 'sequential',
                fxn = function(current) {
                    proposal <- current
                    proposal[on + 1] <- proposal[on + 1] + rnorm(1, mean = 0, sd = sdProps[on + 1])
                    on <<- (on+1) %% nfitted
                    list(proposal=proposal, onpar=on+1, type = 'sequential')
                }))
}

block.proposer <- function(sdProps) return(sdProps=sdProps,
                                           fxn = function(current) {
                                               nfitted <- length(sdProps)
                                               proposal <- current + rnorm(nfitted, mean = 0, sd = sdProps)
                                               list(proposal=proposal,  type = 'block')
                                           })

multiv.proposer <- function(covar, blockLS = list(rownames(covar))) {
    nblocks <- length(blockLS)
    on <- 0
    return(list(type = 'block',
                fxn = function(current) {
                    proposal <- current + rmnorm(1, mean = 0, varcov = covar)
                    propsosal <- as.vector(proposal)
                    names(proposal) <- names(current)
                    list(proposal=proposal,  type = 'block', covar=covar, onpar = NA)
                }))
       }




