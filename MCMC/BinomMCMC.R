library(boot); library(parallel)
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

logPrior <- function(params) ## uniform prior from 0 to 1 on prevalence
    dunif(params, min = 0, max = 1, log = T)

logLikelihood <- function(params, data = list(size=size, sampPos=sampPos))
    dbinom(data$sampPos, size = data$size, prob = params, log = T)

logPosterior <- function(params, data = list(size=size, sampPos=sampPos))
    if(logPrior(params[1])==-Inf) {
        return(-Inf)
    }else{
        return(logPrior(params) + logLikelihood(params, data))
    }

posterior <- function(params, data = list(size=size, sampPos=sampPos))
    exp(logPosterior(params, data))

par(mfrow = c(4,1))
curve(logPrior, 0, 1)
curve(logLikelihood, 0, 1)
curve(logPosterior, 0, 1)
curve(posterior, 0, 1)

runMCMC <- function(iterations, startvalue = runif(1), proposer = gaussianProposal){
    chain <- array(dim = c(iterations+1, length(startvalue)))
    chain[1,] <- startvalue
    for(ii in 1:iterations){
        proposal <- proposer(chain[ii,])
        MHratio <- exp(logPosterior(proposal) - logPosterior(chain[ii,]))
        if(runif(1) < MHratio){
            chain[ii+1,] <- proposal
        }else{
            chain[ii+1,] <- chain[ii,]
        }
    }
    return(chain)
}


chain1 <- runMCMC(1000)
chain2 <- runMCMC(1000)

nchains <- 4
plot(runMCMC(1000), type = 'l', xlab = 'iteration', ylab = 'prevalence', col = rainbow(nchains)[1])
for(ii in 2:nchains) lines(runMCMC(1000), col = rainbow(nchains)[ii])
