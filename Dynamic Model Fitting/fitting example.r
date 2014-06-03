## Example of fitting prevalence data with observation error only.
## Steve Bellan, Meaningful Modeling of Epidemiologic Data 2010
## AIMS, Muizenberg

## First we need to create a system of ODEs.  The disease we are
## thinking of is something like influenza where infecteds recover to
## become immune and the birth & death rates are negligible on the
## time scale of the epidemic.  We also will consider the fact that
## not all of the cases will be reported.  So in essence, we only
## observe each case with a specified but unknown probability.

## The model will be an SEIR, meaning that a latent period is
## explicitly considered.  The C state variable is simply a way to
## keep track of cumulative incidence which is more difficult to do
## outside the model.
seir <- function(t,y,params)
  {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]                           # Cumulative Incidence
    
    beta <- params["beta"]
    N <- params["N"]
    mu <- params["mu"]
    sigma <- params["sigma"]
    gamma <- params["gamma"]
    nu <- mu*N
    
    dSdt <- nu-beta*S*I-mu*S
    dEdt <- beta*S*I-mu*E-sigma*E
    dIdt <- sigma*E-mu*I-gamma*I
    dRdt <- gamma*I-mu*R
    dCdt <- sigma*E
    
    return(list(c(dSdt,dEdt,dIdt,dRdt,dCdt)))
  }

N <- 10000                             # Population of 10,000

param.vals <- c(beta = .00004,
                N = N,
                mu = 0,                 # Assume death rate is
                                        # negligible on time scale of
                                        # outbreak.
                sigma = 1/3,
                gamma = 1/4)

times <- seq(0,300,by = 1)

S0 <- N-1
E0 <- 20
I0 <- 0
C0 <- 0

init <- c(sus = S0,exp = E0,inf = I0,rec = N-S0-E0-I0, cinc = C0)

tc <- data.frame(lsoda(init,times,seir,param.vals))
head(tc)

## Calculate incidence from cumulative incidence: We know that at time
## zero the incidence is 1 since that's when the first exposed case is
## introduced.  The rest of the incidence values are just equal to the
## cumulative incidence at time (t+1) minus those at time (t).
inc <- c(1,tc$cinc[2:nrow(tc)] - tc$cinc[1:(nrow(tc)-1)])

## Let's first plot E(t), I(t), and daily incidence(t).
par(mfrow=c(1,2))
plot(tc$time,
     tc$inf,
     type = "l",
     lwd = 3,
     col = "orange",
     bty = "n",
     xlab = "time",
     ylab = "# of people")

lines(tc$time,
      tc$exp,
      col = "purple",
      type = "l",
      lwd = 3)

legend("topright",
       c("exposed","infected"),
       bty = "n",
       lwd = 3,
       col = c("purple","orange"))
         
plot(tc$time,
     inc,
     col = "red",
     type = "h",
     lwd = 2,
     bty = "n",
     xlab = "time",
     ylab = "# of new cases")

## Let's assume that only a certain proportion of cases are reported.
## This could be due to a combination of both asymptomatic cases and
## people not going to doctors. We'll call the probability of a case
## being reported repp (reported probability).

repp <- .6
rep.case <- rbinom(length(inc), size = round(inc), prob = repp)

lines(tc$time,
      rep.case,
      col = "black",
      type = "h",
      lwd = 2)

## Note that we always observe <100% of the total # of cases.  In
## actuality, we always observe on average repp %.  But the variation
## will be much greater when the # of cases is smaller.

## So now we have this incidence plot, how can we estimate parameters
## of interest from our ODE model?  Pretend that we didn't actually
## simulate this data, but it was given to us and we were asked to
## estimate beta, given that we knew the latent and infectious periods
## exactly.  What we actually have to do is estimate both beta and the
## probability of a case being reported.

## As a start, we need to write a -log likelihood function that gives the
## probability that a given (beta,repp) would generate the observed
## data.n

nll <- function(repp,
                inc.rep,               # Incidence from data
                params,
                init,
                times)
                
  {
    ## Recall that we just need to sum the negative log likelihoods of
    ## observing the observed number of cases given the total number
    ## of cases (coming from the model & our guessed beta) and given
    ## our guessed repp. To do this we run the whole simulation for
    ## the given parameters and then extract the incidence from the
    ## model.

    ## If beta is independently specified, change its value to that in
    ## the parms vector.  This allows us to have beta as an argument
    ## of nll() which is necessary for the below optimization routine.

    dat <- data.frame(lsoda(init,times,seir,param.vals))

    inc.mod <- c(dat$cinc[2:nrow(dat)] - dat$cinc[1:(nrow(dat)-1)])

    ## Then we figure out the probability of observing the # of cases
    ## we did given the incidence from the model & the probability of
    ## observing a case for each data point.  Take the -log and sum to
    ## get the -log(likelihood)!
    nll <- - sum(dbinom(inc.rep,round(inc.mod),repp))

    return(nll)
  }

## That's it!  Writing likelihood functions is surprisingly easy for
## observation error cases like this. So let's see what the likelihood
## of the true model was given the true beta and value of repp.

nll(repp = repp,
    inc.rep = rep.case,
    params = param.vals,
    init = init,
    times = times)


## This number isn't all that helpful.  Likelihoods as numbers do not
## give any information, they only make sense when compared to other
## likelihoods.  So let's optimize the model. First let's assume that
## we know repp and we just must estimate beta.

## To do this we use R's 1-dimensional optimizer:
?optimize

optimize(f = nll,                       # function to maximize/minimize,
         interval = c(0,1),            # interval over which to search
         ## Then give all the other parameters that are fixed that
         ## feed into nll.  Beta is set to x, and optimize wil
         ## optimize over x.
         inc.rep = rep.case,
         params = param.vals,
         init = init,
         times = times,
         ## Then specify other arguments of optimize:
         maximum = FALSE,               # if FALSE, minimize
         tol = 10^-6)                   # desired accuracy

## Compare with the true value:
repp

## Notice that with daily incidence data we can figure out the
## asymptomatic proportion extremeley well.  Try and see what happens
## if you look at weekly data instead.  To do this, just change the
## ODE's time series points to be at 7-day intervals instead of 1-day
## intervals.

######################################################################
## Fitting multiple parameters
###################################################################### 
## Also note that it is rare that we know all the parameters for a
## disease except for one.  Let's pretend we know the latent &
## infectious periods very well, but we do not know the asymptomatic
## proportion (repp) or the transmission coefficient (beta).  Now we
## have to fit both at the same time.

## Because optimize() is only meant for 1-dimensional optimization, we
## know have to switch to optim() which can optimize over
## multi-dimensional parameter space.  As such we should change our
## objective function slightly so that both beta and repp can be in
## the first argument.      

nll.multi <- function(fitted.pars,
                      inc.rep,               # Incidence from data
                      params,
                      init,
                      times,
                      browse = FALSE)   # Debugging option.
  {
    if(browse)  browser()
    ## Extract parameters from the fitted parameters vector and put
    ## into the appropriate objects.
    params["beta"] <- fitted.pars["beta"]
    repp <- fitted.pars["repp"]

    ## Run ODE & extract incidence from time series.
    dat <- data.frame(lsoda(init,times,seir,params))
    inc.mod <- c(0,dat$cinc[2:nrow(dat)] - dat$cinc[1:(nrow(dat)-1)])
    inc.mod <- round(inc.mod)

    ## Only calculate NLL if there are less reported cases than cases
    ## from the model.  Any model that predicts less actual cases than
    ## reported cases must be thrown out.
    if(sum(inc.rep>inc.mod) > 0)
      {
        nll <- -130
      }else{
        
        ## Calculate negative log likelihood.  Make sure to not have
        ## negative values from the ODE model (which can occur since the
        ## optimization routine may pick weird values for the parameters
        ## to try out.

        inc.mod[inc.mod<0] <- 0
        nll <- - sum(dbinom(round(inc.rep),inc.mod,repp))

        ## We know that our parameters cannot take on any values so we add
        ## a penalty function that makes our nll big if they are outside
        ## of the correct range.  This will make the values bounce back on
        ## the next iteration.

        if(params["beta"] < 0)      nll <- nll + 20
        if(repp < 0)                nll <- nll + 20
        if(repp > 1)                nll <- nll + 20
      }

    return(nll)
  }


## Let's find out the likelihood of the true parameter values given
## our data.
nll.multi(c(param.vals["beta"],repp = repp),
          inc.rep = rep.case,               # Incidence from data
          params = param.vals,
          init = init,
          times = times,
          browse = F)

nll.multi(optim.vals$par,
          inc.rep = rep.case,               # Incidence from data
          params = param.vals,
          init = init,
          times = times)

## Check out optim()'s help file.
?optim

## Select initial values for fitted parameters from which optimization
## routine will start. If you select bad initial values the algorithm
## can get stuck on a bad set of parameters.
init.fitted <- c(beta = .00005, repp = .3)
## You can always try the true values as a starting point for this
## problem, although that's rarely possible in real problems.
## init.fitted <- c(param.vals["beta"],repp=repp)

optim.vals <- optim(par = init.fitted,
                    nll.multi,
                    inc.rep = rep.case,
                    params = param.vals,
                    init = init,
                    times = times,
                    method = "Nelder-Mead") # Try SANN as well, it much slower.
optim.vals

## Note how sensitive this optimization is to the initial conditions
## of beta and repp that we gave it. When they were way off it got
## stuck due to the penalties we imposed in the nll.multi() function.
## When it tried one guess, it got -130 then it moved slightly in
## parameter space and still got -130.  So thinking that it had
## reached a flat plane, the algorithm ended.

## Let's plot our true data vs the fitted plot. To do this we first
## must run a simulation with our optimized parameter values.
fitted.param.vals <- param.vals
fitted.param.vals["beta"] <- optim.vals$par["beta"]
fitted.sim <- data.frame(lsoda(init,times,seir,fitted.param.vals))
## Calculate incidence
inc.fitted <- c(0,fitted.sim$cinc[2:nrow(fitted.sim)] - fitted.sim$cinc[1:(nrow(fitted.sim)-1)])

## Incidence

plot(times,
     inc,
     type = "s",
     col = "red",
     bty = "n",
     ylim = c(0, max(inc.fitted, inc)),
     xlab = "time",
     lwd = 3,
     ylab = "incidence",
     main = "incidence (symptomatic & asymptomatic)")

lines(times,
      inc.fitted,
      type ="s",
      col ="black",
      lwd = 3)

legend("topright",
       c("data","fitted"),
       col = c("red","black"),
       lwd = 3)

######################################################################
## Likelihood Profiles
###################################################################### 
## Next step is to profile each of our parameters over the other.
## This means that we will maximize the likelihood for repp =
## seq(0,1,by=.01) where beta is fitted, and vice versa.  Then we can
## see how likely the values are for one parameter when the other is
## completely loose.

## This code is not yet finished! Try completing it yourself:

beta.seq <- seq(0,.01, length.out = 200)
beta.prof <- rep(NA, length(beta.seq))

for(ii in 1:length(beta.seq))
  {
    
