library(deSolve); library(ggplot2)
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
    c(log(.01),log(2))
   ,c(log(.7), log(100))
   ,c(log(1),log(1/10))))
colnames(initBounds) <- c('lower','upper')
rownames(initBounds) <- c('Beta','alpha','progRt')
class(initBounds[,2]) <- class(initBounds[,3]) <- 'numeric'
initBounds

initRand <- function(fit.params) {
    fit.params <- logParms(fit.params)
    tempnm <- names(fit.params)
    for(nm in tempnm) fit.params[nm] <- runif(1, min = initBounds[rownames(initBounds)==nm, 'lower'], 
                                              max =  initBounds[row.names(initBounds)==nm, 'upper'])
    return(unlogParms(fit.params))
}
initRand(c(alpha = 3, Beta = 1))

## MCMC
mcmcSampler <- function(current.params, ref.params=disease_params(), obsDat, seed = 1, proposer = sequential.proposer(sdProps=sdProps),
                        randInit = T, niter = 100, verbose=0) {
    if(verbose>0) browser()
    if(randInit) current.params <- initRand(current.params)
    ii <- 2
    accept <- 0
    curVal <- llikePrior(current.params, ref.params = ref.params, obsDat=obsDat, verbose = verbose)
    out <- matrix(NA, nr = niter, nc=length(current.params)+1)
    out[1,] <- c(current.params, nll = -curVal)
}

mcmcSampler(list(alpha = 1, Beta = 1), obsDat=obsDat, verbose = 6)

multiv.proposer <- function(covar, blockLS = list(rownames(covar))) {
    nblocks <- length(blockLS)
    on <- 0
    return(function(current) {
        proposal <- current
        onBlock <- blockLS[[on + 1]]
        proposal[onBlock] <- current[onBlock] + rmnorm(1, mean = 0, varcov = covar[onBlock,onBlock])
        proposal <- as.vector(proposal) #otherwise is a matrix
        names(proposal) <- names(current)
        on <<- (on+1) %% nblocks
        proposal
    })
}

block.proposer <- function(sdProps) return(function(current) {
    nfitted <- length(sdProps)
    current + rnorm(nfitted, mean = 0, sd = sdProps)
})

sequential.proposer <- function(sdProps) {
    nfitted <- length(sdProps)
    on <- 0
    return(function(current) {
  	proposal <- current
  	proposal[on + 1] <- proposal[on + 1] + rnorm(1, mean = 0, sd = sdProps[on + 1])
  	on <<- (on+1) %% nfitted
  	proposal
    })
}


