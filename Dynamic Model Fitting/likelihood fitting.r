### remove all items from R's workspace to start clean
rm(list=ls())

### load the ODE solving package
library(odesolve)

### note that all graphical outputs and saved data
### will be saved in the working directory which
### you can find by typing
### getwd()

### the data set:
### make sure to change your working directory to a directory that
### contains this file.  You can do this with the
### setwd() function

data <- read.csv("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/Harare.prev.with.ci.csv")

## Convert years to month & year in POSIXct
years <- data[,1]
round.year <- floor(years)
month <- years - round.year
month <- month*12+1
date.strg <- paste(round.year,"-",month,sep="","-01")
date <- strptime(date.strg, format="%Y-%m-%d")
date <- as.POSIXct(date)

data[,2:4] <- data[,2:4]/100
ANC.prev <- data[,2]

### find n for the 95% confidence intervals
### we are just using a standard calculation for
### approximate normal confidence intervals for binomial distributions
### with reasonably large n*p (e.g. sample size * prevalence)
average.95 <- (data[,2]-data[,3] + data[,4]-data[,2] )/2
norm.sd <- average.95/1.96
norm.var <- norm.sd^2

ANC.n <- (ANC.prev * (1-ANC.prev)) / norm.var
ANC.n <- round(ANC.n)
dat.len <- length(ANC.n)

### Plot the Harare data


tseq.plot <- seq(1990,2007,by=.25)
tseq.plot.ind <- rep(NA, dat.len)

for(i in 1:dat.len)	tseq.plot.ind[i] <- which(tseq.plot==years[i]) 

ANC.cnt	<-	round(ANC.prev*ANC.n)
bars.dat		<-	rbind(ANC.cnt, ANC.n-ANC.cnt)

pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/ANC N.pdf", width=10, height=6)
par(mar=c(6,7,4,4))
    
plot(date,
     ANC.prev*100,
     type = "p",
     ylab = "prevalence %",
##      main = "Harare ANC HIV Prevalence",
     bty = "n",
     pch = 19,
     col = "red",
     cex = 1.2,
     xlim=as.POSIXct(c("1985-01-01","2010-01-01")),
     ylim = c(0,40),
     cex.lab = 2.5,
     cex.main = 2.5,
     cex.axis = 2)

arrows(date, 100*data[,3],
       date, 100*data[,4],
       length = .05,
       code=3,
       col = "red",
       lwd = 1.5,
       angle=90)

dev.off()

## bars <- matrix(0,nrow=2, ncol=length(tseq.plot))
## bars[,tseq.plot.ind] <- bars.dat
## barplot(bars,
##         xlab="years",
##         ylab="sample size",
##         col=c("red","blue"),
##         names.arg=tseq.plot,
##         cex.lab=1.3,
##         cex.main=3,
##         main="Harare data by year")

legend("topright", c("HIV -", "HIV +"), col=c("blue","red"), bty="n", pch=15, bg=c("blue","red"), cex=3)

dev.off()



### the ODE model
ode.mod		<-	function(tt,	yy,	parms)
{
	bb	<-	parms["bb"]	### birth rate
	mu	<-	parms["mu"]	### death rate
	gg	<-	parms["gg"]	### Erlang rate for transition b/w HIV stages

	ff	<-	parms["ff"]	### infectivity
	aa	<-	parms["aa"]	### response to prevalence
	qq	<-	parms["qq"]	### response to mortality

	SS	<-	yy[1]		### Susceptibles
	I1	<-	yy[2]		### HIV stage 1
	I2	<-	yy[3]		### HIV stage 2
	I3	<-	yy[4]		### HIV stage 3
	I4	<-	yy[5]		### HIV stage 4

	II	<-	I1+I2+I3+I4		### total infecteds
	NN	<-	II + SS				### total population
	mort	<-	gg * I4 / NN	### HIV-related mortality
									### Force of infection
	FF	<-	ff * exp(-aa * II/NN) * exp(-qq * mort)

	deriv	<-	rep(NA,5)	### ODE system

	deriv[1]<-	bb*NN - mu*SS - FF*SS*II/NN
	deriv[2]<-	FF*SS*II/NN - gg*I1 - mu*I1
	deriv[3]<-	gg*I1 - gg*I2 - mu*I2
	deriv[4]<-	gg*I2 - gg*I3 - mu*I3
	deriv[5]<-	gg*I3 - gg*I4 - mu*I4

list(c(deriv,NULL))
	}

### Time steps to solve over
### It is important that these steps
### include the years that your data occur in
tseq		<-	seq(1970, 2040, by=.05)

### Find which time steps match the ANC surveys 
### for model fitting
dat.len		<-	dim(data)[1]
tseq.ind	<-	rep(NA, dat.len)
for(i in 1:dat.len)	tseq.ind[i]	<-	which(tseq==years[i]) 

### Create an objective function that will be minimized
### to fit the model to the data.  Here we choose parameters,
### run the model, then calculate the sum of the squared
### differences between the model predictions and the data.
### R's function optim will then search through the parameters
### we want to fit to find those that minimize this sum
### of squared differences (SoS).

obj.fun		<-	function(pars.optim, pars.fixed, binom.fit=TRUE, hessian=FALSE)
	{

	### All this code below does is determine which parameter
	### vectors contain which parameters.  This makes this function
	### flexible to handle two vectors of parameters
	### (one to be fixed & one to be fit) that contain any
	### any of the parameters you need to run the model.
	
	if(!is.na(pars.fixed["bb"]))	bb<-pars.fixed["bb"]
	if(!is.na(pars.fixed["mu"]))	mu<-pars.fixed["mu"]
	if(!is.na(pars.fixed["gg"]))	gg<-pars.fixed["gg"]
	if(!is.na(pars.fixed["SS"]))	SS<-pars.fixed["SS"]
	if(!is.na(pars.fixed["I2"]))	I2<-pars.fixed["I2"]
	if(!is.na(pars.fixed["I3"]))	I3<-pars.fixed["I3"]
	if(!is.na(pars.fixed["I4"]))	I4<-pars.fixed["I4"]
	if(!is.na(pars.fixed["ff"]))	ff<-pars.fixed["ff"]
	if(!is.na(pars.fixed["aa"]))	aa<-pars.fixed["aa"]
	if(!is.na(pars.fixed["qq"]))	qq<-pars.fixed["qq"]
	if(!is.na(pars.fixed["I1"]))	I1<-abs(pars.fixed["I1"])
	
	if(!is.na(pars.optim["bb"]))	bb<-pars.optim["bb"]
	if(!is.na(pars.optim["mu"]))	mu<-pars.optim["mu"]
	if(!is.na(pars.optim["gg"]))	gg<-pars.optim["gg"]
	if(!is.na(pars.optim["SS"]))	SS<-pars.optim["SS"]
	if(!is.na(pars.optim["I2"]))	I2<-pars.optim["I2"]
	if(!is.na(pars.optim["I3"]))	I3<-pars.optim["I3"]
	if(!is.na(pars.optim["I4"]))	I4<-pars.optim["I4"]
	if(!is.na(pars.optim["ff"]))	ff<-pars.optim["ff"]
	if(!is.na(pars.optim["aa"]))	aa<-pars.optim["aa"]
	if(!is.na(pars.optim["qq"]))	qq<-pars.optim["qq"]
	if(!is.na(pars.optim["I1"]))	I1<-abs(pars.optim["I1"])
	

	### now that you have extracted the parameters & 
	### initial conditions put
	### them all into one vector which is how the
	### odesolver lsoda() needs them to be inputted
	parms	<-	c(bb, mu, gg, ff, aa, qq)
	init	<-	c(SS=SS, I1=I1, I2=I2, I3=I3, I4=I4)
	
	### this is the line that actualy runs a dynamical simulation
	sim.dat	<-	lsoda(init, tseq, ode.mod, parms=parms)

	### now we take the prevalence values out of our ode simulation's
	### output.  Then we create a vector of prevalences from the years 
	### for which we have data.
	prev	<-	rowSums(sim.dat[,3:6]) / rowSums(sim.dat[,2:6])
	prev.yr	<-	prev[tseq.ind]

	### by default we are fitting binomial likelihoods but if 
	### you set binom.fti=FALSE then it will do a least squares fit
	### which is basically equivalent to calculating a likelihood
	### for normally distributed random variables.
	if(!binom.fit)	value	<-	sum( (ANC.prev - prev.yr)^2)
	if(binom.fit	)	value	<-	-sum( dbinom(round(ANC.prev*ANC.n), ANC.n, prev.yr, log=TRUE) )

	### Return the likelihood as output of obj.fun()
	return(value)
	}



### Next wWe create a function that can take initial parameters,
### optimize the function and then plot the optimized model.
### This function will also be run inside the likelihood profile function
### below which is why we include a plot argument.  If plot=TRUE this function,
### when run, will plot a simulation output and fit to the data.
### "prof.par" gives the parameter you will profile if input from
### prof.lik() function built below.


plot.fun	<-	function(init.pars.opt, init.pars.fixed, binom.fit=TRUE, plot=FALSE, prof.par=NA, hessian=FALSE)
	{
	
	### If there is a parameter to profile
	### (that is if prof.par is not NA)
	if(!is.na(prof.par))
	{
	
	
	### The next three lines delete the parameter to profile
	### from the optimize vector and adds it to the fixed vector.
	### The profiled parameter will be set to the value set by
	### the proflik function below.

	par.ind	<-	which(names(init.pars.opt)==prof.par)
	
	init.pars.fixed	<-	c(init.pars.fixed, init.pars.opt[par.ind])
	init.pars.opt	<-	init.pars.opt[-par.ind]
	}
	
	### This line of code optimizes the objective function over the init.par.opt
	### parameters.
	fitted.pars	<-	optim(par=init.pars.opt, obj.fun, pars.fixed=init.pars.fixed, method="Nelder-Mead",control=list(maxit=2000, parscale=init.pars.opt), hessian=hessian)
					
			
	if(plot)		### If you want to plot the simulation output.
	{	
	pars.optim	<-	fitted.pars$par
	pars.fixed	<-	init.pars.fixed

	bb	<-	pars.fixed["bb"]
	mu	<-	pars.fixed["mu"]
	gg	<-	pars.fixed["gg"]
	SS	<-	pars.fixed["SS"]
	I2	<-	pars.fixed["I2"]
	I3	<-	pars.fixed["I3"]
	I4	<-	pars.fixed["I4"]

	ff	<-	pars.optim["ff"]
	aa	<-	pars.optim["aa"]
	qq	<-	pars.optim["qq"]
	I1	<-	pars.optim["I1"]

	parms	<-	c(bb, mu, gg, ff, aa, qq)
	init	<-	c(SS=SS, I1=I1, I2=I2, I3=I3, I4=I4)

	sim.dat	<-	lsoda(init, tseq, ode.mod, parms=parms)

	SS	<-	sim.dat[,2]
	II	<-	rowSums(sim.dat[,3:6])
	I4	<-	sim.dat[,6]
	NN	<-	rowSums(sim.dat[,2:6])
	prev	<-	II / NN * 100
	mort	<-	gg*I4

	inc	<-	ff * exp(-aa * II/NN) * exp(-qq * mort) * SS*II/NN


	title				<-	"OLS fit"
	if(binom.fit)	title	<- "Maximum Likelihood Fit with Binomial Distr"

	pdf("Harare Binom 4.pdf", width=12, height=8)
	par(mfrow=c(1,2))
	plot(years, ANC.prev*100, xlab="years", ylab="HIV prevalence (%)", pch=19, ylim=c(0,100), xlim=c(1970,2040), main=title)
	lines(tseq, prev, col="red", type="l", lwd=3)
	arrows(years, 100*data[,3], years, 100*data[,4], length = .05, code=3, angle=90)

	legend("topright", bty="n", c("Uganda ANC prevalence"), pch=19)
	legend("topleft", bty="n", c("fitted prevalence"), lwd=3, col="red")

	ymax	<-	1.2*max(inc, mort)*10^5

	plot(tseq, mort*10^5, xlab="years", ylab="per 100,000", type="l", lwd=3, col="black",ylim=c(0,ymax))
	lines(tseq, inc*10^5, type="l", col="green", lwd=3)
	legend("topright", c("fitted incidence", "fitted mortality"), lwd=3, col=c("green","black"), bty="n")

	dev.off()
	}	# end plot if statement
	
	return(fitted.pars)		### This returns the output of the optimization

	}	# end plot.fun function




### This is the function that will calculate profile likelihoods
### for a specified parameter.  prof.seq is the values of the profiled
### parameter that that parameter will be constrained to.  A likelihood
### profile constrains one value of a parameter while fitting all the others
### and then shows you the values of the likelihood for different constraints
### of that parameter.

prof.lik		<-		function(init.pars.opt, init.pars.fix, prof.par, plot=FALSE, prof.seq)
{
	
### for loop for each constrained value of the profiled parameter
for(ii in 1:length(prof.seq) )
{
	### Sets the profiled parameter to its i-th value
	init.pars.opt[prof.par]	<-	prof.seq[ii]
	
	### Creates MLE fit to model for i-th value of profiled parameter.
	output	<-	plot.fun(	init.pars.opt=init.pars.opt, 
									init.pars.fixed=init.pars.fixed, 
									binom.fit=binom.fit,
									prof.par=prof.par,
									plot=plot)
### if you are doing the first iteration then set up the data frame
### to store all the output you are interested in (parameters that are fit,
### resulting likelihoods, whether the optimization converged, and how
### many iterations it took)
	if(ii==1) {
					num.pars		<-		length(output$par)
					num.rows		<-		num.pars + 3
					par.names	<-		names(output$par)
					names			<-		c(par.names, "likelihood", "iterations", "convergence")
		pars.frame	<-	data.frame(matrix(NA,nrow=num.rows,ncol=length(prof.seq)), row.names=names)
					}
					
		names(pars.frame)[ii]		<-paste("a = ", prof.seq[ii])
		
### take the output from the plot.fun() and store
### it in the data frame you have set up		
	pars.frame[1:num.pars,ii]	<-	output$par
	pars.frame[num.pars+1, ii]	<-	output$value
	pars.frame[num.pars+2, ii]	<-	round(output$counts[1])
	pars.frame[num.pars+3, ii]	<-	round(output$convergence)
	
### this line lets you know what iteration R is working on
### as it fits it because it can take some time and you
### want to know that its going.	
print(paste("iteration #",ii))
} #end for loop

return( pars.frame )
} #end proflik function
#### These are the initial values of the parameters and 
#### initial conditions that we will use for all fits by default
#### note that some or all of the first vector while be adjusted by the
#### optim() function during the fitting process.
init.pars.opt	<-	c(ff=.5, aa=8, qq=300, I1=.001)
init.pars.fix	<-	c(bb=.029, mu=.018, gg=.303, SS=1, I2=0, I3=0, I4=0)

### This line of code states the values of 
### the parameter to be profiled for which
### the other parameters will be fit
### and the likelihood calculated.
prof.seq		<-	seq(-12, 12, length.out=100)
prof.par		<-	"aa"
save(prof.seq, file=paste("prof.seq", prof.par))

### This line calculates on likelihood profile and will take a lot of time
### depending on the length of the prof.seq (e.g. # of values over which
### to profile the parameter)
profile.dat	<-	prof.lik(init.pars.opt, init.pars.fix, prof.par=prof.par, plot=FALSE, prof.seq)

### Because that took so long you want to save that data frame
### to have a look at the data frame type "profile.dat" in R's
### main window.
save(profile.dat, file=paste("prof",prof.par))


### Use these lines if you want to load a profile you've already
### calculated without redoing it.
load(paste("prof",prof.par))
load(paste("prof.seq",prof.par))

### This line optimizes with all loose parameters
#mle.fit	<-	plot.fun(init.pars.opt, init.pars.fixed, binom.fit=TRUE, plot=FALSE, prof.par=NA, hessian=F)

### What is the minimum -log(likelihood) (e.g. max(likelihood))
mle.val	<-	min(profile.dat[4,])

### which iteration gave the minimum(-log(likelihood))
mle.ind	<-	which(profile.dat[4,]==mle.val)

### what parameter value for the profiled parameter gave
### the minimum(-log(likelihood))
mle.par	<-	prof.seq[mle.ind]


### The Likelihood Ratio Test (LRT) says that
### -2*log( maximum likelihood / likelihood of constrained model)
### is approximately Chi squared distributed with degrees
### of freedom equal to the # of parameters constrained in the constrained model.
### To get 95% confidence intervals we calculate the difference
### in -log(likelihoods) between the maximum likelihood estimate
### and the -log(likelihood) of a constrained fit that would
### be at the 95% critical value.

### The 95% Chi Squared cutoff is
#qchisq(.95, df=1, lower.tail=TRUE)

### AFter some manipulation we get the following
like95	<-	mle.val + qchisq(.95,1, lower.tail=TRUE)/2

### Set scale of y axis for plotting
ymax	<-	max(profile.dat[4,])


### A bunch of plots to play with.

### Plot the likelihood profile you created.
pdf(paste("likelihood profile & collinearity for",prof.par,".pdf"), width=14, height=10)

par(mfrow=c(2,2))
plot(prof.seq, profile.dat[4,], xlab="heterogeneity parameter 'a'", ylab="negative log likelihood", type="l",col="black", lwd=3, main="profile likelihood surface", ylim=c(0,ymax), cex.main=2, cex.lab=1.3)
points(mle.par,mle.val, pch=20, cex=3,col="green")
text(mle.par,.7*mle.val,"MLE", cex=1.5)
abline(h=like95, lwd=3, col="green")
legend("topright","95% CI",lwd=3, col="green", bty="n", cex=2)

### Plot the values of the other fitted parameters vs
### the constrained values of the profiled parameter

### these 3 lines make sure that the axis labels are
### autmoatically set
par.labels	<-	c(ff="force of infection (f)", aa="heterogeneity (a)", qq="response to mortality", I1="initial # of infecteds (I1(0))")
which.pars.fitted	<-	row.names(profile.dat)[1:3]

par.labels.fitted	<-	par.labels[which.pars.fitted]

for(jj in 1:3)
{
plot(prof.seq, profile.dat[jj,], xlab=par.labels[prof.par], ylab=par.labels.fitted[jj], pch=19, cex.lab=1.3)
scale.lik	<-	profile.dat[4,]/max(profile.dat[4,])*max(profile.dat[jj,])
lines(prof.seq, scale.lik, lty=2, col="purple",lwd=2)

#if(jj==2){legend("topright", paste("profile likelihood of ",par.labels[prof.par]), col="purple", lwd=2, lty=2, bty="n",cex=1.5)}
}

dev.off()


