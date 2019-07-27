library(deSolve); library(ggplot2); library(MASS); library(sfsmisc); library(mnormt); library(ellipse); library(emdbook)
library(boot); library(parallel); library(animation); library(coda)

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
sdmid <- median(sdvec)
set.seed(5)
test <- runMCMC(7, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                logPriorFxn = get(ff), startvalue = logit(.25), noProps=T, plotNM = paste0(ffp,'-',sdmid))
## Final chains
set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                noProps=T, plotNM = paste0(ffp,'-',sdmid), iiShow=3000)

set.seed(4)
test <- runMCMC(nits, plotter=mcmcHist, verbose = 0, proposer=unifProposal(halfwidth=hw), 
                noProps=T, plotNM = paste0(ffp,'-',sdmid), iiShow=3000, burn = 400)
