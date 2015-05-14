setwd('~/Documents/R Repos/lectureCode/MCMC/')
source('BinomMCMC.R')

## Uniform proposer
ff <- 'logUnifPrior'
ffp <- paste0(ff, '-UnifProp')
nits <- 3000
hw <- .5
## Movie
nm <- paste0('movies/',ffp,'-',sdvec[ii], '.mov')
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
sdvec <- c(.05, .2, 2)
ff <- 'logUnifPrior'
ii <- 2
nits <- 30

for(ff in c('logUnifPrior','logBetaPrior')) {
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


runMCMC(100, verbose = 0, proposer=gaussianProposal(sd=sdvec[ii]), startvalue = logit(.8)) 

    
nchains <- 4
its <- 1000
par(mfrow=c(1,1))
plot(0, type = 'n', xlab = 'iteration', ylab = 'prevalence', col = rainbow(nchains)[1], ylim = c(0,1), xlim = c(0, its))
for(ii in 1:nchains) lines(runMCMC(1000), col = rainbow(nchains)[ii])


