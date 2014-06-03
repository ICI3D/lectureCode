## Introduction to Likelihood
## Meaningful Modeling of Epidemiolgic Data, 2014
## (C) Steve Bellan, 2009-2014
#################### 
## Order of lecture
#################### 
## 1) Link to binomial sampling in Jim Scott's lecture
## 
## 2) Show binomial sampling and variety of null hypotheses & probability
## of data given each of these
## 
## 3) Talk about p-values, and alpha cutoffs.
## 
## 4) Show exact binomial confidence intervals, and define them as the
## collection of null hypotheses for which a data point this extreme
## or more extreme would not occur less than 5% of the time.
## 
## 5) Define likelihood.
## 
## 6) Explain link between P(data | model) and L(model | data).
## 
## 7) Derive likelihood for binomial distribution to get x/n as best
## estimate of p.
## 
## 8) Note that -2log(L) ~ chisq(df = 1)
## 
## 9) Show how this can be used to produce confidence intervals similar
## to the exact binomial ones.

##################################################
## Plotting parameters
fgcol <- 'white'
brdcol <- 'black' ## barplot bar borders
hcol <- gray(.5) ## hypothetical distribution colors
xcol <- 'red' ## extreme value in tail color
pcol <- 'purple' ## color of probability densities or exact values
ycol <- 'yellow' ## color of emphasized text
ps <- 24
cx <- 1.2 ## ps/12 
ymax <- .15

##################################################
## Bin(100, .3), choose 28
for(ss in c(F,T)) {
  ## Population parameters for ANC clinic as a subset of a greater
  ## population example.
  size <- 100
  true.prev <- .3
  ## Sample from this distribution once (use 28 for powerpoint slides).
  samp.pos <- 28                          #rbinom(1,100,.3)
  samp.prev <- samp.pos/size
  color.vec <- rep(fgcol,size+1)
  ## more.extreme <- 0:100 <= samp.pos | 0:100 >= (abs(true.prev * size - samp.pos) + true.prev * size)
  if(ss) {
    color.vec[0:100<=samp.pos] <- xcol
    strg <- paste0(' x=', samp.pos)
  }else{
    strg <- ''
  }
### Create a histogram of the binomial distribution function
  pdf(paste0("BinPDF p=", true.prev, strg,'.pdf'), width=14, height=7)
  par(fg = fgcol, col.main = fgcol, col.axis=fgcol, col.lab=fgcol, 'ps'=ps, mgp = c(6,1.5,0), mar = c(5.3,7.8,3,4))
  mpar <- par()
  barplot.obj <- barplot(dbinom(0:size,size=size, prob=true.prev),
                         #xlab="number HIV+", #names.arg=c(0:size),
                         ylab="probability",
                         col=color.vec,
                         border = brdcol, #fgcol,
                         space =0, axis.lty=0,
                         ## main=paste("probability distribution function for  a true prevalence of", true.prev*100,"%"),
                         ylim = c(0,ymax), axes=F, xlim = c(0,100),
                         cex.axis=cx, cex.main=cx, cex.lab=cx)
  axis(1, seq(0,100, by = 10))
  axis(2, seq(0,.1, l = 3), las = 2)
  mtext("number HIV+", side = 1, line = 3.5, cex = cx)
  dev.off()
}

######################################################################
## p-values as cumulative two tails for different null hypotheses.
######################################################################
potential.prev <- seq(.15,.4, by = .05)
length.out <- length(potential.prev)
strg <- paste0(' x=', samp.pos)
file.names <- paste0("BinPDF p=", potential.prev, strg, ".pdf")
for(indiv.pnls in c(F,T)) {
  if(!indiv.pnls) {
    pdf(paste0("BinPDF p", paste0(range(potential.prev), collapse = '-'), strg, "panels.pdf"),width=12,height=8)
    suppressWarnings(par(mpar))
    par(mfrow=c(2,3), oma = c(2,2,0,0), mar = c(5.1,4.1,4,1))
  }
  for(ii in 1:length.out) {
    if(indiv.pnls)      {       pdf(file.names[ii], width=14, height=7); suppressWarnings(par(mpar))}
    ## Now calculate the total probability of finding a value in the
    ## extreme lower tail plus the extreme upper tail.
    if(samp.pos < potential.prev[ii]*size) {
      p.val.one.tail <- pbinom(samp.pos, size, potential.prev[ii], lower.tail = TRUE)
      blurb <- paste("2*pbinom(",samp.pos,", ", size, ", ", potential.prev[ii], ", lower.tail = TRUE)",
                     sep = "")
      more.extreme <- 0:100 <= samp.pos
      loc2x <- qbinom(.1, size, potential.prev[ii])*.7
    }else{
      p.val.one.tail <- pbinom(samp.pos - 1, size, potential.prev[ii], lower.tail = FALSE)
      blurb <- paste("2*pbinom(",samp.pos,"-1, ", size, ", ", potential.prev[ii], ", lower.tail = FALSE)",
                     sep = "")
      more.extreme <- 0:100 >= samp.pos
      loc2x <- qbinom(.9, size, potential.prev[ii])*1.5
    }
    main <- ifelse(indiv.pnls,
                   paste("If true prevalence were ", potential.prev[ii]*100,"%, then p(28 or more extreme) is",sep=""),
                   paste("hypothetical \nprevalence:",potential.prev[ii]*100,"%"))
    sub <- ifelse(indiv.pnls, blurb, '')
    ## Now multiply times 2, because we could have gotten reults that
    ## extreme in the other direction too. Make sure you understand what
    ## this means.
    p.val <- p.val.one.tail*2
    blurbP <- signif(p.val,3)
    ## if(potential.prev[ii]==.3) col.main <- "yellow"
    prob.dist <- dbinom(0:size, size=size, prob=potential.prev[ii])
    color.vec <- rep(hcol,size+1)
    color.vec[more.extreme] <- xcol
    barplot(prob.dist, xlab = "", ylab = ifelse(indiv.pnls, "probability",''),
            names.arg = NA, col = color.vec , border = brdcol, space = 0,
            main = main, axes = F,
            sub = blurb,
            col.main = fgcol, ylim = c(0,ymax))
    text(loc2x, .025, '2X', col = xcol, cex = 1.5)
    axis(1, seq(0,100, by = 10), ifelse(indiv.pnls, T, NA))
    axis(2, seq(0,.1, l = 3), ifelse(indiv.pnls, T, NA), las = 2)
    mtext("number HIV+", side = 1, line = ifelse(indiv.pnls, 3.5, 0), outer = !indiv.pnls, cex = cx)
    if(!indiv.pnls)     mtext("probability", side = 2, line = 0, outer = T, cex = cx)
    if(indiv.pnls)      mtext(blurb, side = 3, line = -2, col = "yellow", cex = cx)
    mtext(paste("p =",blurbP), side = 3, line = -4, col = xcol, cex = cx)
    if(indiv.pnls) dev.off()
  }
  dev.off()
}

######################################################################
## p-value vs hypothesis plot, note you can't use qbinom()! it doesn't
## give CI's unlike qnorm because of lack of symmetry
######################################################################
pdf("P value Confidence Intervals.pdf", width = 14, height = 7)
suppressWarnings(par(mpar)); par(mar = c(8,9,1,1))
potential.prev.vector <- seq(0,1, by = .001)
p.val <- rep(NA, length(potential.prev.vector))
for(jj in 1:length(potential.prev.vector))  {
    if(samp.pos < potential.prev.vector[jj]*size) {
        p.val.one.tail <- pbinom(samp.pos, size, potential.prev.vector[jj], lower.tail = TRUE)
      }else{
        p.val.one.tail <- pbinom(samp.pos - 1, size, potential.prev.vector[jj], lower.tail = FALSE)
      }
    ## Now multiply times 2, because we could have gotten reults that
    ## extreme in the other direction too. Make sure you understand what
    ## this means.
    p.val[jj] <- p.val.one.tail*2
  }
non.rejected.nulls <- potential.prev.vector[p.val > .05]
ci <- range(non.rejected.nulls)
plot(potential.prev.vector, p.val,
     xlab = "hypothetical prevalence", ylab = "p-value", las = 2,
     lwd = 4, bty = "n", type = "l", col = xcol)
segments(0,.05, 1,.05, lty = 2, lwd = 3)
cip.l <- round(min(potential.prev.vector[p.val> .05]),3)
cip.u <- round(max(potential.prev.vector[p.val> .05]),3)
arrows(cip.l, .4, cip.l, -.03, length = .2)
arrows(cip.u, .4, cip.u, -.03, length = .2)
text(cip.l,.4, cip.l,cex = 1.5,pos = 3)
text(cip.u,.4, cip.u,cex = 1.5,pos = 3)
text(.5, .05, expression(alpha==.05), pos = 3)
text(.7,.7, paste("95% CI includes HIV prevalences \nof ", cip.l*100, "% to ", cip.u*100, "%", sep=""), cex = 1.5,
     col = "yellow")
dev.off()

###################################################################### 
### LIKELIHOOD APPROACH
######################################################################
for(ss in c(F,T)) {
  color.vec <- rep(fgcol,size+1)
  if(ss) {
    color.vec[samp.pos + 1] <- pcol
    strg <- paste0(' x=', samp.pos)
  }else{
    strg <- ''
  }
### Create a histogram of the binomial distribution function
  pdf(paste0("BinLIK p=", true.prev, strg,'.pdf'), width=14, height=7)
  suppressWarnings(par(mpar))
  barplot.obj <- barplot(dbinom(0:size,size=size, prob=true.prev),
                         #xlab="number HIV+", #names.arg=c(0:size),
                         ylab="probability",
                         col=color.vec,
                         border = brdcol, #fgcol,
                         space =0, axis.lty=0,
                         ## main=paste("probability distribution function for  a true prevalence of", true.prev*100,"%"),
                         ylim = c(0,ymax), axes=F, xlim = c(0,100),
                         cex.axis=cx, cex.main=cx, cex.lab=cx)
  axis(1, seq(0,100, by = 10))
  axis(2, seq(0,.1, l = 3), las = 2)
  mtext("number HIV+", side = 1, line = 3.5, cex = cx)
  dev.off()
}

######################################################################
## Likelihood values for different null hypotheses.
######################################################################
potential.prev <- seq(.15,.4, by = .05)
length.out <- length(potential.prev)
strg <- paste0(' x=', samp.pos)
file.names <- paste0("BinLIK p=", potential.prev, strg, ".pdf")
for(indiv.pnls in c(F,T)) {
  if(!indiv.pnls) {
    pdf(paste0("BinLIK p=", paste0(range(potential.prev), collapse = '-'), strg, "panels.pdf"),width=12,height=8)
    suppressWarnings(par(mpar))
    par(mfrow=c(2,3), oma = c(2,2,0,0), mar = c(5.1,4.1,4,1))
  }
  for(ii in 1:length.out) {
    if(indiv.pnls)      {       pdf(file.names[ii], width=14, height=7); suppressWarnings(par(mpar))}
    ## Now calculate the total probability of finding a value in the
    ## extreme lower tail plus the extreme upper tail.
    dbn <- signif(dbinom(samp.pos, size, potential.prev[ii]),3)
    blurb <- paste0("dbinom(",samp.pos,", ", size, ", ", potential.prev[ii], ') = ', dbn)
    main <- paste("hypothetical", "\n"[!indiv.pnls], " prevalence:",potential.prev[ii]*100,"%")  
    sub <- ifelse(indiv.pnls, blurb, '')
    ## if(potential.prev[ii]==.3) col.main <- "yellow"
    prob.dist <- dbinom(0:size, size=size, prob=potential.prev[ii])
    color.vec <- rep(hcol,size+1)
    color.vec[samp.pos + 1] <- pcol
    barplot(prob.dist, xlab = "", ylab = ifelse(indiv.pnls, "probability",''),
            names.arg = NA, col = color.vec , border = brdcol, space = 0,
            main = main, axes = F, sub = blurb, ylim = c(0,ymax))
    axis(1, seq(0,100, by = 10), ifelse(indiv.pnls, T, NA))
    axis(2, seq(0,.1, l = 3), ifelse(indiv.pnls, T, NA), las = 2)
    mtext("number HIV+", side = 1, line = ifelse(indiv.pnls, 3.5, 0), outer = !indiv.pnls, cex = cx)
    if(!indiv.pnls)     mtext("probability", side = 2, line = 0, outer = T, cex = cx)
    if(indiv.pnls)      mtext(blurb, side = 3, line = -3.5, col = pcol, cex = cx)
    if(!indiv.pnls)     mtext(dbn, side = 3, line = -4, col = pcol, cex = cx) 
    if(indiv.pnls) dev.off()
  }
  if(!indiv.pnls)       dev.off()
}

######################################################################
## Show likelihood plot
###################################################################### 
######################################################################
## Likelihood-based CI's
## give CI's unlike qnorm because of lack of symmetry
######################################################################
pdf("Likelihood Confidence Intervals.pdf", width = 14, height = 7)
suppressWarnings(par(mpar)); par(mar = c(8,9,5,1))
potential.prev.vector <- seq(0,1, by = .001)
likelihoods <- dbinom(samp.prev*size,size = size, prob = potential.prev.vector)
plot(potential.prev.vector, likelihoods,     
     xlab = "hypothetical prevalence", ylab = "likelihood", las = 1,
     lwd = 4, bty = "n", type = "l", col = pcol)
mtext("p(our data given prevalence) = LIKELIHOOD", side = 3, line = 3, cex= cx)
dev.off()

pdf("logLikelihood Confidence Intervals.pdf", width = 14, height = 7)
suppressWarnings(par(mpar)); par(mar = c(8,9,5,1))
potential.prev.vector <- seq(0,1, by = .001)
likelihoods <- dbinom(samp.prev*size,size = size, prob = potential.prev.vector)
plot(potential.prev.vector, log(likelihoods), 
     xlab = "hypothetical prevalence", ylab = "log(likelihood)", las = 1,
     lwd = 4, bty = "n", type = "l", col = pcol)
axis(2, seq(-400,0, l = 5))
mtext("p(our data given prevalence) = LIKELIHOOD", side = 3, line = 3, cex= cx)
dev.off()

######################################################################
## L to logL to -logL
######################################################################
pdf("Likelihood Confidence Intervals (small).pdf", width = 5, height = 3)
suppressWarnings(par(mpar)); par(mar = c(4,4,1,1), mgp = c(2,1,0))
plot(potential.prev.vector, likelihoods,     
     xlab = "hypothetical prevalence", ylab = "likelihood", las = 1,
     axes = F,
     lwd = 4, bty = "n", type = "l", col = pcol)
axis(1, seq(0,1, by = .2), NA)
axis(2, seq(0,.08, by = .02), NA)
mtext("p(our data given prevalence) = LIKELIHOOD", side = 3, line = 3, cex= cx)
dev.off()

pdf("logLikelihood Confidence Intervals (small).pdf", width = 5, height = 3)
suppressWarnings(par(mpar)); par(mar = c(4,4,1,1), mgp = c(2,1,0))
plot(potential.prev.vector, log(likelihoods), 
     xlab = "hypothetical prevalence", ylab = "log(likelihood)", las = 1,
     axes = F,
     lwd = 4, bty = "n", type = "l", col = pcol)
axis(1, seq(0,1, by = .2), NA)
axis(2, seq(-400,0, l = 5), NA)
mtext("p(our data given prevalence) = LIKELIHOOD", side = 3, line = 3, cex= cx)
dev.off()


pdf("NeglogLikelihood Confidence Intervals (small).pdf", width = 5, height = 3)
suppressWarnings(par(mpar)); par(mar = c(4,4,1,1), mgp = c(2,1,0))
plot(potential.prev.vector, -log(likelihoods), 
     xlab = "hypothetical prevalence", ylab = "-log(likelihood)", las = 1,
     axes = F,
     lwd = 4, bty = "n", type = "l", col = pcol)
axis(1, seq(0,1, by = .2), NA)
axis(2, seq(0,400, l = 5), NA)
mtext("p(our data given prevalence) = LIKELIHOOD", side = 3, line = 3, cex= cx)
dev.off()


######################################################################
## Likelihood-based CI's BIG
## give CI's unlike qnorm because of lack of symmetry
######################################################################
pdf("Likelihood Confidence Intervals big.pdf", width = 14, height = 7)
suppressWarnings(par(mpar)); par(mar = c(8,9,5,1))
potential.prev.vector <- seq(0,1, by = .001)
likelihoods <- dbinom(samp.prev*size,size = size, prob = potential.prev.vector)
plot(potential.prev.vector, likelihoods,     
     xlab = "hypothetical prevalence", ylab = "likelihood", las = 1,
     axes=F,
     lwd = 4, bty = "n", type = "l", col = pcol)
axis(1, seq(0,1, by = .2))
axis(2, seq(0,.08, by = .02), las = 2)
mtext("Likelihood", side = 3, line = 3, cex= cx)
dev.off()

for(crit in c(T,F)) {
  pdf(paste0("NeglogLikelihood Confidence Intervals", 'crit'[crit], "big.pdf"), width = 14, height = 7)
  suppressWarnings(par(mpar)); par(mar = c(8,9,5,1))
  potential.prev.vector <- seq(0,1, by = .001)
  likelihoods <- dbinom(samp.prev*size,size = size, prob = potential.prev.vector)
  plot(potential.prev.vector, -log(likelihoods), 
       xlab = "hypothetical prevalence", ylab = "-log(likelihood)", las = 1, axes=F,
       lwd = 4, bty = "n", type = "l", col = pcol)
  axis(1, seq(0,1, by = .2))
  axis(2, seq(0,400, l = 5), las = 2)
  mtext("we usually minimize the -log(likelihood)", side = 3, line = 3, cex= cx)
  if(crit) {
    chisq.crit <- qchisq(.95, df = 1)
    min.l <- min(-log(likelihoods))
    abline(h = min.l + chisq.crit/2,
           lwd = 3,
           lty = 2,
           col = fgcol)
  }
  dev.off()
}

##################################################
## Chi squared cutoff
######################################################################
## Chi sq plot
###################################################################### 
pdf("chisq.pdf",width = 8, height = 6)
suppressWarnings(par(mpar)); par(mar=c(12,9,4,2), 'ps'=25)
plot.seq <- seq(0,5,by = .001)
above.crit <- plot.seq > qchisq(.95,1)
col.vec <- rep(fgcol, length(plot.seq))
col.vec[above.crit] <- xcol
plot(plot.seq, dchisq(plot.seq, df = 1),
     type = "h", col = col.vec,
     ylim=c(0,3), bty = "n",
     main = expression(paste("PDF for", chi["df=1"] ^2)),
     xlab = "", ylab = "probability density")
mtext(expression(paste(-2*log(frac(L[null],L[alternative])))), side = 1, line = 8)
arrows(qchisq(.95,1), 1, qchisq(.95,1), dchisq(qchisq(.95,1),1)+.1, length = .2, lwd = 2)
text(qchisq(.95,1),1, expression(paste(chi["df=1",alpha==.05] ^ 2)), pos = 3)
dev.off()


######################################################################
## Zoom in on last plot
for(ii in 1:3) {
  pdf(paste0("NeglogLikelihood Confidence Intervals", 'crit zoom', ii, ".pdf"), width = 14, height = 7)
  suppressWarnings(par(mpar)); par(mar = c(8,9,5,1))
  potential.prev.vector <- seq(0,1, by = .001)
  likelihoods <- dbinom(samp.prev*size,size = size, prob = potential.prev.vector)
  plot(potential.prev.vector, -log(likelihoods), ylim = c(2,7),
       xlab = "hypothetical prevalence", ylab = "-log(likelihood)", las = 1, axes=F,
       lwd = 4, bty = "n", type = "l", col = pcol)
  axis(1, seq(0,1, by = .2))
  axis(2, seq(2,6, b = 2), las = 2)
  chisq.crit <- qchisq(.95, df = 1)
  min.l <- min(-log(likelihoods))
  abline(h = min.l + chisq.crit/2, lwd = 3, lty = 2, col = fgcol)
  if(ii==2) {
    arrows(.28, min.l, .28, min.l+chisq.crit/2, len = .15, code = 3, angle = 45, lwd = 4, col = ycol)
    text(.28, min.l+chisq.crit/4, '1.92', pos = 4, col = ycol)
  }
  if(ii==3) {
    ci.likelihood <- range(potential.prev.vector[-log(likelihoods) < min.l + chisq.crit/2])
    ci.l <- signif(ci.likelihood[1],3)
    ci.u <- signif(ci.likelihood[2],3)
    mtext(paste("95% CI includes HIV prevalences of ", ci.l*100, "% to ", ci.u*100, "%", sep=""), side = 3, line =0,
          col = ycol)
    arrows(ci.l, min.l + 1.92, ci.l, min.l, length = .2, col = ycol)
    arrows(ci.u, min.l + 1.92, ci.u, min.l, length = .2, col = ycol)
    text(ci.l, min.l + 1.92, ci.l,cex = 1,pos = 3, col = ycol)
    text(ci.u, min.l + 1.92, ci.u,cex = 1,pos = 3, col = ycol)
  }
dev.off()
}
 
######################################################################
## Let's also make a 2 panel plot so we can compare this with the
## P-Value approach from above.
pdf("Both CIs.pdf", width = 12, h= 8)
suppressWarnings(par(mpar)); par(mfrow=c(2,1), mar = c(8,8,3,4), 'ps' = 20)
xlim <- c(0, 1)
## Redraw the P-Value plot from above
## Plot the p-values vs the null hypotheses.
plot(potential.prev.vector, p.val, xlab = "hypothetical prevalence (null hypothesis)", ylab = "p-value",
     xlim = xlim, lwd = 2, bty = "n", type = "l",  col = xcol, las = 1, yaxt='n')
axis(2, seq(0,1,l=3), las = 2)
## Draw a line at the alpha = .05 cutoff.
segments(0, .05, 1, .05, lty = 2, lwd = 3)
## Now plot the CI bounds
arrows(cip.l, .4, cip.l, -.03, length = .2)
arrows(cip.u, .4, cip.u, -.03, length = .2)
text(cip.l,.4, cip.l,cex = 1,pos = 3)
text(cip.u,.4, cip.u,cex = 1,pos = 3)
text(.4,.2, paste("95% CI includes ", cip.l*100, "% to ", cip.u*100, "%", sep=""),
     col = "yellow", cex = 1, pos = 4)
plot(potential.prev.vector, -log(likelihoods),
     type = "l", col = pcol, lwd = 3, xlim = xlim, las = 1,
     ylim = c(2, 7), bty = "n", yaxt='n',
     xlab = "potential prevalences (our models)", ylab = "-log(likelihood)",
     main = "")
axis(2, seq(2,6,b=2), las = 2)
abline(h = min.l+chisq.crit/2, lwd = 3, lty = 2)
text(.4, 5, paste("95% CI includes ", ci.l*100, "% to ", ci.u*100, "%", sep=""),
     col = "yellow", cex = 1, pos = 4)
arrows(ci.l, min.l + 3, ci.l, min.l, length = .2)
arrows(ci.u, min.l + 3, ci.u, min.l, length = .2)
text(ci.l, min.l + 3, ci.l,cex = 1,pos = 3)
text(ci.u, min.l + 3, ci.u,cex = 1,pos = 3)
dev.off()
