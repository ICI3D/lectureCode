## Make some distributions to  show on slides.



pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/BinomDistr.pdf", width = 10, height =8)

par(mar = c(5,7,6,4))
binom.seq <- 0:100
barplot(dbinom(binom.seq, size = 100, prob = .4),
        names.arg = binom.seq,
        col = "black",
        space = 0,
        xlab = "# successes in N trials",
        ylab = "probability",
        cex.lab = 2.5,
##         main = "Binomial Distribution", 
        xaxt = "n",
        cex.main = 4,
        col.main = "yellow")
dev.off()



pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/NormDistr.pdf", width = 10, height =8)
par(mar = c(5,7,6,4))
norm.seq <- seq(-3,3,length.out=3000)
barplot(dnorm(norm.seq),
        names.arg = norm.seq,
        col = "black",
        space = 0,
        xlab = "(approximately) continuous variable",
        ylab = "probability",
        cex.lab = 2.5,
##         main = "Normal Distribution",
        xaxt = "n",
        cex.main = 4,        
        col.main = "yellow")
dev.off()

pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/ExpDistr.pdf", width = 10, height =8)
par(mar = c(5,7,6,4))
exp.seq <- seq(0, 4, length.out = 2000)
barplot(dexp(exp.seq),
        names.arg = exp.seq,
        col = "black",
        space = 0,
        xlab = "time until event",
        ylab = "probability",
        cex.lab = 2.5,
##         main = "Exponential Distribution",
        xaxt = "n",
        cex.main = 4,        
        col.main = "yellow")
dev.off()

pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/PoisDistr.pdf", width = 10, height =8)
par(mar = c(5,7,6,4))
pois.seq <- 0:15
barplot(dpois(pois.seq,lambda=3),
        names.arg = pois.seq,
        col = "black",
        ylim = c(0,.3),
        space = 0,
        xlab = "# of events in time interval",
        ylab = "probability",
        cex.lab = 2.5,
##         main = "Poisson Distribution",
        xaxt = "n",
        cex.main = 4,        
        col.main = "yellow")

dev.off()

pdf("~/Documents/R files/MMED 2010/Lecture Code/Dynamic Model Fitting/fake binom data.pdf", width = 10, height = 8)
par(mar = c(5,7,6,4))

rbinom(10, mean = c(
