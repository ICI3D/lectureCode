## Makes a figure that displays difference between RR & AR
## setwd('~/Documents/R Repos/lectureCode/Study Design/')
backCol <- 'black'
mainCol <- 'white'
opar <- par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(5,6,1,2))


dat <- cbind(c(100,0),c(100,10))

pdf("AR & RR.pdf", width = 6, height = 6)
par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(4,7,4,4))
barplot(dat,
        beside = FALSE,
        col=c("grey","red"),
        main = "Relative Risk (RR) & Attributable Risk (AR)",
        ylab = "flu incidence per 1000 people",
        ylim = c(0, 150),
        space = 1,
        cex.lab = 1.5,
        cex.axis = 1.5,
        names.arg = c("> 8 hours of sleep", "< 8 hours of sleep"))
legend("topleft", c("Attributable Risk = 10 per 1,000 people","","Relative Risk = 110/100 = 1.1")
       , pch = c(19,NA,NA)
       , col = "red"
       , bty = "n")
dev.off()


## Make a fake outbreak

num.days <- 22
lambda = .5*exp(.1*1:num.days)
inc <- rpois(n = num.days, lambda)
ci <- sum(inc)
inc <- c(inc,0,ci)

total <- 60
sus <- c(rep(0,length(inc)-1),total-ci)

first.day <- as.POSIXct("2010-05-24")
days <- seq(first.day, first.day + (num.days-1)*3600*24, by = "day")
days <- format(days, "%d-%m")
days[1:length(days) %% 2 == 0] <- ''

pdf("fake incidence.pdf", width = 7, height = 5)
par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(4,7,4,4))
par(mar=c(6,6,1,2), 'ps' = 12)
barplot(rbind(inc,sus),
        col = c("red","grey"),
        names.arg = c(days, "","TOTAL"),
        beside = FALSE,
        las = 2,
        xlab = "",
        cex.lab = 1,
        cex.axis = 1, 
        ylab = "")
mtext('incidence', 2, 3, cex = 1.5, las = 3)
cum.inc.str <- substitute(paste("Cum Incidence = ",
    frac(paste(CI, " cases"), paste(TOTAL, " people")), " = ", CI.perc, "% over ", NUM.D, " days" , sep ="")
                          ,env = list(CI=ci, TOTAL=total, CI.perc = round(ci/total*100,0), NUM.D = num.days)
                          )
inc.den.str <- substitute(paste("Inc Density = ",
    frac(paste(CI, " cases", sep=""), paste("(",TOTAL, " people * ", NUM.D, " days)", sep="")),
    " = ", INC.D, frac(cases,"person*day")),
                          env = list(CI = ci,
                            TOTAL = total,
                            NUM.D = num.days,
                            INC.D =round(ci/(total*num.days),2)))
text(0,total*.8,cum.inc.str, pos = 4, cex = 1.2)
text(0,total*.5,inc.den.str, pos = 4, cex = 1.2)
dev.off()

## Prevalence Plot
pdf("fake prevalence.pdf", width = 5, height = 4)
par(bg=backCol,fg=mainCol, lwd=2, col.axis=mainCol, col.lab=mainCol, col = mainCol, col.main=mainCol, 
            cex.axis=1.5, cex.lab=1.5, 'las'=1, bty='n', 'mgp'=c(4,1,0), mar = c(4,7,4,4))
par(mar=c(6,6,4,2))
size <- 60
prev <- .1
pos <- rbinom(n=1, size = size, p = prev)
dat <- cbind(c(pos, size-pos), c(0,0))
barplot(dat,
        beside = F,
        col = c("red", "grey"),
        names.arg = c("Prevalence Study",""),
        ylab = "People",
        cex.axis = 1.5, border = NA)
legend("topright", c("Healthy", "Disease"), col=c("grey","red"), pch = 19, bty = "n", cex = 1.4)
text(2,size/3, paste("Prevalence = \n", pos, "/", size, " = ", round(pos/size*100,0), "%", sep=""), cex=1)
dev.off()
