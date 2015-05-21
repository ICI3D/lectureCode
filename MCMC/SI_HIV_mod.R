library(deSolve); library(ggplot2); library(MASS); library(sfsmisc); library(mnormt)
setwd('~/Documents/R Repos/lectureCode/MCMC/')
source('utilityFxns.R')

opar <- par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,2))

disease_params <- function(Beta = 0.3, alpha = 1, progRt = 1/2.5,
                           birthRt = 1/40, deathRt = birthRt, ...)
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

sampleEpidemic <- function(simDat, tseq = seq(1978, 2010, by = 2), numSamp = rep(150, length(tseq)), verbose=0) {
    if(verbose>0) browser()
    simDat$I <- rowSums(simDat[, Is])
    prev_at_sample_times <- simDat[simDat$time %in% tseq, 'I']
    numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
    lci <- mapply(function(x,n) binom.test(x,n)$conf.int[1], x = numPos, n = numSamp)
    uci <- mapply(function(x,n) binom.test(x,n)$conf.int[2], x = numPos, n = numSamp)    
    return(data.frame(time = tseq, numPos, numSamp, sampPrev =  numPos/numSamp,
                      lci = lci, uci = uci))
}


simDat <- as.data.frame(lsoda(init, tseq, SImod, parms=disease_params(Beta = .9, alpha = 8, progRt = 1/2.5)))
simDat$I <- rowSums(simDat[, Is])
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
    c(log(.06),log(2)) ## beta
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


mcmcSampler <- function(current.params, ref.params=disease_params(), obsDat, seed = 1, proposer = sequential.proposer(sdProps=sdProps),
                        plotter = plotterTS, randInit = T, niter = 100, nburn = 0, verbose=0, plotNM=NULL, tell = 100) {
    if(verbose>2) browser()
    if(randInit) current.params <- initRand(current.params)
    vv <- 2 ## mcmc iteration
    accept <- 0
    curVal <- llikePrior(current.params, ref.params = ref.params, obsDat=obsDat, verbose = verbose)
    out <- matrix(NA, nr = niter, nc=length(current.params)+1)
    out[1,] <- c(current.params, nll = -curVal)
    colnames(out) <- c(names(current.params), 'nll')
    max_index <- dim(out)[1]
    last.it <- 0
    while(vv <= max_index) {
        if ((verbose > 1) || (verbose && (vv%%tell == 0))) print(paste("on iteration",vv,"of",last.it + niter + 1))
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
        ## if(plotFitted & lmh!=-Inf) do.call(plotout, args=within(plotArgs, {objresult <- res}))
        ##            if(!is.null(plotNM)) 
        ##png(paste0(plotNM,'.png'), width = 800*resScl, height = 600*resScl)
        ## browser()
        if(!is.null(plotter)) {
            par(opar)
            plotter(out, vv, ref.params=ref.params, obsDat=obsDat, proposer = proposer, proposal = proposal, onpar=onpar,proptype=propt)
        }
        ##          do.call(plotter, args=within(plotArgs, {curState <- out[vv,]}))
        ##          if(!is.null(plotNM)) dev.off()

    }
    if(!is.null(plotNM)) graphics.off()
    aratio <- accept/((vv-nburn))
    colnames(out) <- c(names(current.params), 'nll')
    return(list(out = out[1:nrow(out)>(nburn+1),], aratio = aratio, current.params = current.params, ref.params=ref.params))
}

plotterParmDens <- function(out, vv, ref.params=disease_params(), plotNM=NULL, obsDat, verbose=0, proposer = NULL, onpar = onpar, 
                            proptype=proptype,
                            proposal = NA, propDistCol = 'yellow', propCol='brown', curCol = 'dodger blue', every = 200, burn = 100,
                            marLine = 8, lmar=23, ps = 25, xlim = c(1,50), ylim = c(.05,2), log = 'xy', bump = 5, nlevs = 50,
                            yparnm = expression(beta), xparnm=expression(alpha)) {
    out <- out[!is.na(out[,1]), colnames(out) !='nll', drop=F]
    newParms <- out[nrow(out) , colnames(out) !='nll']
    lastParms <- fit.params <- out[nrow(out)-1,]
    out <- out[1:(nrow(out)-1),,drop=F]
    parnms <- colnames(out)
    parms <- within(ref.params, { ## subs fitting parameters into reference parameter vector
        for(nm in parnms) assign(nm, as.numeric(fit.params[nm]))
        rm(nm)
    })
    if(verbose > 7) browser()
    ## simulate model
    simDat <- as.data.frame(lsoda(init, tseq, SImod, parms=parms))
    simDat$I <- rowSums(simDat[, Is])

    layout(matrix(c(3,1,4,0,2,4), 3, 2), w = c(1,.5), h = c(.6,1,1))
    par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
        cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,2), 'ps'=ps)
    par(mar=c(8,lmar,1,1))
    ## rescale axes every so many iterations
    if(vv + 1 > every) {
        axisNum <- floor(vv/every)
        ##if(vv+1==every) browser()
        tout <- out[burn:min(nrow(out), axisNum*every),]
        xlim <- quantile(tout[,1], c(.02,.98)) * c(.9, 1.1)
        ylim <- quantile(tout[,2], c(.02,.98)) * c(.9, 1.1)
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
    if(nrow(out)>20) {
        z <- kde2d(Lout[,1], Lout[,2], lims = c(Lxlim, Lylim), h = .2)
        if(grepl('x',log)) z$x <- exp(z$x)
        if(grepl('y',log)) z$y <- exp(z$y)
        ## Posterior bivariate density
        cols <- apply(colorRamp(c('black','red','orange','white'))(seq(0,1, l = nlevs)), 1, function(x) rgb(x[1],x[2],x[3], max=255))
        .filled.contour(z$x, z$y, z$z, levels = pretty(range(z$z), nlevs, xlim = Lxlim, ylim = Lylim), col = cols)
    }
    if(nrow(out)<60)    points(outOriginal, col = makeTransparent('light green', alpha = 100), cex = 2, pch = 16)
    mtext(xparnm, 1, marLine-3, cex = 1.5)
    mtext(yparnm, 2, marLine+3, cex = 1.5)
    if(proptype=='block') { ## Bivariate proposer
        bivPDF <- function(x,y) dmnorm(cbind(x,y), mean = log(lastParms), varcov = proposer$covar)
        xsAt <- seq(Lxlim[1], Lxlim[2], l=30)
        ysAt <- seq(Lylim[1], Lylim[2], l=30)
        bivDens <- outer(xsAt, ysAt , bivPDF)
        if(grepl('x',log)) xsAt <- exp(xsAt)
        if(grepl('y',log)) ysAt <- exp(ysAt)
        ## Proposal bivariate density
        colsProp <- apply(colorRamp(c('white','yellow'),
                                    bias = 1)(seq(0,1, l = nlevs)), 1, function(x) rgb(x[1],x[2],x[3], max=255))
        colsProp <- diag(makeTransparent(colsProp, c(0, seq(50,200, l = nlevs-1))))
        .filled.contour(xsAt,ysAt, bivDens, levels = pretty(range(bivDens), nlevs, xlim = Lxlim, ylim = Lylim), col = colsProp)
        accepted <- sum(newParms!=proposal)==0
        pchProp <- ifelse(accepted,19,21)
## browser()
##         points(lastParms, pch = 19, col = curCol, cex = 2.5)
##         points(proposal, pch = pchProp, col = propCol, cex = 2.5)

    }

    ## Marginal Histograms
    xbreaks <- seq(Lxr[1], Lxr[2], l=25)
    ybreaks <- seq(Lyr[1], Lyr[2], l=25)
    xhist  <-  hist(Lout[,1], plot=FALSE, breaks = xbreaks)
    yhist  <-  hist(Lout[,2], plot=FALSE, breaks = ybreaks)
    top  <-  max(c(xhist$counts, yhist$counts))

    ## Y
    scl <- .1
    dep <- .1
    par(mar=c(8,bump+2,1,1)) 
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, border = NA)
    par(new=T)
    plot(0,0, type='n',ylim = Lylim, xlim = c(0,1), axes=F, xlab='',ylab='',main='')
    axis(2, at = Lyticks, labels = F)
    if(proptype=='sequential') {
        valSeq <- seq(Lylim[1],Lylim[2], l = 1000)
        densSeq <- dnorm(valSeq, log(lastParms[2]), sd = proposer$sdProps[2])
        propDens <- dnorm(log(proposal[2]), log(lastParms[2]), proposer$sdProps[2])
        valSeqP <- c(valSeq, rev(valSeq))
        densSeqP <- c(-densSeq * scl, rep(0, length(densSeq)))-dep
        par(xpd=NA)
        polygon(densSeqP, valSeqP, col = makeTransparent(propDistCol, alpha = 100), border=NA)
        accepted <- newParms[2]==proposal[2]
        if(accepted) ltyProp <- 1 else ltyProp <- 3
        segments(max(densSeqP), log(lastParms[2]), min(densSeqP), log(lastParms[2]), col = curCol, lwd = 3)
        if(onpar==1) segments(-dep, log(proposal[2]), -propDens*scl-dep, log(proposal[2]), col = propCol, lty = ltyProp, lwd = 3)
        par(xpd=T)
    }

    ## X
    par(mar=c(bump,lmar,1,1))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0, border = NA)
    par(new=T)
    plot(0,0, type='n',xlim = Lxlim, ylim = c(0,1), axes=F, xlab='',ylab='',main='')
    axis(1, at = Lxticks, labels = F)
    if(proptype=='sequential') {
        valSeq <- seq(Lxlim[1],Lxlim[2], l = 1000)
        densSeq <- dnorm(valSeq, log(lastParms[1]), sd = proposer$sdProps[1])
        propDens <- dnorm(log(proposal[1]), log(lastParms[1]), proposer$sdProps[1])
        valSeqP <- c(valSeq, rev(valSeq))
        densSeqP <- c(-densSeq * scl, rep(0, length(densSeq)))-dep
        par(xpd=NA)
        polygon(valSeqP, densSeqP, col = makeTransparent(propDistCol, alpha = 100), border=NA)
        accepted <- newParms[2]==proposal[2]
        if(accepted) ltyProp <- 1 else ltyProp <- 3
        segments(log(lastParms[1]), max(densSeqP), log(lastParms[1]), min(densSeqP), col = curCol, lwd = 3)
        if(onpar==2) segments(log(proposal[1]), -dep, log(proposal[1]), -propDens*scl-dep, col = propCol, lty = ltyProp, lwd = 3)
        par(xpd=T)
    }
    
    ## Time Series
    par(mar=c(6,lmar,1,1))
    plot(simDat$time, simDat$I, xlab = '', ylab = '', type = 'l', ylim = c(0,.4), col='red', lwd = par()$lwd, mgp = c(4,3,0))
    ## add data
    points(obsDat$time, obsDat$sampPrev, col = 'red', pch = 16, cex = 4)
    arrows(obsDat$time, obsDat$uci, obsDat$time, obsDat$lci, col = makeTransparent('red'), len = .025, angle = 90, code = 3, lwd = 3)
    mtext('HIV\nprevalence', 2, marLine+ 7, adj = .5)
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
    return(list(covar = covar, type = 'block',
                fxn = function(current) {
                    proposal <- current + rmnorm(1, mean = 0, varcov = covar)
                    propsosal <- as.vector(proposal)
                    names(proposal) <- names(current)
                    list(proposal=proposal,  type = 'block', sdProps=covar, onpar = NA)
                }))
       }




