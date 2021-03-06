source('BinomMCMC.R')
source('SI_HIV_mod.R')

## Uniform proposer
sdvec <- c(.05, .2, 2)
ff <- 'logUnifPrior'
ffp <- paste0(ff, '-UnifProp')
nits <- 3000
hw <- .5

####################################################################################################
## Univariate Sampling
####################################################################################################
## Uniform Proposer

## Movie
nm <- paste0('movies/',ffp,'-',hw, '.mov')
if(file.exists(nm)) file.remove(nm)
set.seed(4)
saveVideo({
    ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
    runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), logPriorFxn = get(ff))
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 800*resScl, ani.height = 600*resScl)

## Stills
set.seed(5)
test <- runMCMC(7, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                logPriorFxn = get(ff), startvalue = logit(.25), noProps=T, plotNM = paste0(ffp,'-',sdvec[ii]))
## Final chains
set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                noProps=T, plotNM = paste0(ffp,'-',sdvec[ii]), iiShow=3000)

set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                noProps=T, plotNM = paste0(ffp,'-',sdvec[ii]), iiShow=3000, burn = 400)

## Gaussian proposer
ff <- 'logUnifPrior'
ii <- 2
nits <- 30
for(ff in c('logUnifPrior','logBetaPrior')) { ## with uniform & beta prior
    for(ii in 1:length(sdvec)) {
        set.seed(4)
        nm <- paste0('movies/',ff,'-',sdvec[ii], '.mov')
        if(file.exists(nm)) file.remove(nm)
        saveVideo({
            ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
            runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), logPriorFxn = get(ff))
        }, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 800*resScl, ani.height = 600*resScl)
        set.seed(5)
        test <- runMCMC(7, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), startvalue = logit(.25),
                        noProps=T, plotNM = paste0(ff,'-',sdvec[ii]))
    }
}
## Final chains
set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), 
                noProps=T, plotNM = paste0(ff,'-',sdvec[ii]), iiShow=3000)

set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), 
                noProps=T, plotNM = paste0(ff,'-',sdvec[ii]), iiShow=3000, burn = 400)

## Trace plots
ff <- 'logUnifPrior'
ffp <- paste0(ff, '-Gaus')
ii <- 2
nits <- 6000
nm <- paste0('movies/',ffp,'-trace',sdvec[ii],'.mov')
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.003, nmax = 300, ani.dev='png', ani.type='png')
    runMCMC(nits, plotter=mcmcHistTrace, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), startvalue = logit(.8)) 
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 800*resScl, ani.height = 600*resScl)

## Trace with multiple chains
chainLs <- list()
nchains <- 8
for(ii in 1:nchains) {
     chainLs[[ii]] <- runMCMC(5000, verbose = 0, proposer=gaussianProposal(sd=.05))
}

nm <- paste0('movies/',ffp,'-traceMultiChains',sdvec[ii],'.mov')
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.06, nmax = 300, ani.dev='png', ani.type='png')
    for(jj in seq(10, 4000, by = 20))   tracePlot(chainLs[1:4], jj, alpha = 75)#, plotNM='test')
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 800*resScl, ani.height = 600*resScl)

####################################################################################################
## Bivariate fit to SI HIV model
####################################################################################################

set.seed(7)
## Bivariate sequential sampler
nm <- paste0('movies/','HIV-BivariateSeq', '.mov')
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
    mcmcSampler(c(alpha=8, Beta=.9), ref.params=disease_params(), obsDat, seed = 1
                , proposer = sequential.proposer(sdProps=c(.15,.15)),
                , plotterLoops = list(repeats = c(20, 10, 1), breaks=c(12, 20,  10^5))
                , plotter = plotterParmDens, randInit = T, niter = 1200, nburn = 0, verbose=0, plotNM=NULL)
},
          video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 700*resScl, ani.height = 700*resScl)

## Stills for sequential sampler
for(rr in 1:6) {
    png(paste0('movies/BivariateSeqStills/BivariateSeq-',rr,'.png'), width = 700*resScl, height = 700*resScl) #
    set.seed(7)
    mcmcSampler(c(alpha=8, Beta=.9), ref.params=disease_params(), obsDat, seed = 1
              , proposer = sequential.proposer(sdProps=c(.15,.15))
              , plotter = plotterParmDens, randInit = T, niter = 2, nburn = 0, verbose=0, plotNM=NULL
               ,showing = showingFXN(shDat=rr>1, shCurr=rr>2, shPropKern =rr>3, shProp=rr>4,
                    shHist=rr>5, shPost=F, shAratio = F))
    dev.off()
}

## Block sampler
nm <- paste0('movies/','HIV-block', '.mov')
set.seed(7)
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
    mcmcSampler(c(alpha=4, Beta=.9), ref.params=disease_params(), obsDat, seed = 1, 
                proposer = multiv.proposer(covar = matrix(c(.02,.00,.00,.02),2,2)),
                plotterLoops = list(repeats = c(20, 10, 1), breaks=c(12, 20,  10^5)),
                plotter = plotterParmDens, randInit = T, niter = 1200, nburn = 0, verbose=0, plotNM=NULL)
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 700*resScl, ani.height = 700*resScl)


## Adaptive Block sampler
nm <- paste0('movies/','HIV-blockAdapt', '.mov')
set.seed(7)
if(file.exists(nm)) file.remove(nm)
saveVideo({
    ani.options(interval = 0.02, nmax = 300, ani.dev='png', ani.type='png')
    mcmcSampler(c(alpha=4, Beta=.9), ref.params=disease_params(), obsDat, seed = 1, 
                proposer = multiv.proposer(covar = matrix(c(.02,.00,.00,.02),2,2)), adaptiveMCMC = T,
                startAdapt = 300, adptBurn = 200,
                plotterLoops = list(repeats = c(20, 10, 1), breaks=c(12, 20,  10^5)),
                plotter = plotterParmDens, randInit = T, niter = 1200, nburn = 0, verbose=0, plotNM=NULL)
}, video.name = nm, other.opts = "-b 3000k -pix_fmt yuv420p", ani.width = 700*resScl, ani.height = 700*resScl)

## Stills for adaptive block sampler
for(rr in 1:6) {
    png(paste0('movies/BivariateAdaptStills/BivariateAdapt-',rr,'.png'), width = 700*resScl, height = 700*resScl) #
    set.seed(7)
    mcmcSampler(c(alpha=4, Beta=.9), ref.params=disease_params(), obsDat, seed = 1, 
                proposer = multiv.proposer(covar = matrix(c(.02,.00,.00,.02),2,2)), adaptiveMCMC = T,
                startAdapt = 300, adptBurn = 200,
                showing = showingFXN(shDat=rr>1, shCurr=rr>2, shLast = rr>2, shPropKern =rr>3, shProp=rr>4,
                    shHist=rr>5, shPost=F, shAratio = F),
                plotter = plotterParmDens, randInit = T, niter = 2, nburn = 0, verbose=0, plotNM=NULL)
    dev.off()
}
