source("DOprob.R")

# load data
r <- seq(0, 0.49, by=0.01)
oRA <- array(dim=c(12,11,50))
dimnames(oRA) <- list(1:12, 0:10, myround(r, 2))
oRXf <- oRXm <- oRA
for(i in 1:12) {
  cat(i,"\n")
  attachfile(i, "diversity_outcross_results")
  oRA[i,,] <- cbind(0, t(sapply(gA,  apply, 3, function(a) 1-sum(diag(a))/sum(a))))
  oRXf[i,,] <- cbind(0, t(sapply(gXf,  apply, 3, function(a) 1-sum(diag(a))/sum(a)))) 
  oRXm[i,,] <- cbind(0, t(sapply(gXm,  apply, 3, function(a) 1-sum(diag(a))/sum(a)))) 
  detach(2)
}

RA <- RXf <- RXm <- oRA
for(k in 1:12) {
  RA[k,1,] <- rikRprob(r, k, "autosome")
  RXf[k,1,] <- rikRprob(r, k, "femaleX")
  RXm[k,1,] <- rikRprob(r, k, "maleX")
  for(s in 1:10) {
    RA[k,s+1,] <- doRprob(r, k, s, "autosome")
    RXf[k,s+1,] <- doRprob(r, k, s, "femaleX")
    RXm[k,s+1,] <- doRprob(r, k, s, "maleX")
  }
}

max((abs(oRA - RA)/sqrt(RA*(1-RA)/160000)), na.rm=TRUE)
max((abs(oRXf - RXf)/sqrt(RXf*(1-RXf)/80000)), na.rm=TRUE)
max((abs(oRXm - RXm)/sqrt(RXm*(1-RXm)/40000)), na.rm=TRUE)
median((abs(oRA - RA)/sqrt(RA*(1-RA)/160000)), na.rm=TRUE)
median((abs(oRXf - RXf)/sqrt(RXf*(1-RXf)/80000)), na.rm=TRUE)
median((abs(oRXm - RXm)/sqrt(RXm*(1-RXm)/40000)), na.rm=TRUE)

hist(((oRA - RA)/sqrt(RA*(1-RA)/160000)), breaks=150)
hist(((oRXf - RXf)/sqrt(RXf*(1-RXf)/80000)), breaks=150)
hist(((oRXm - RXm)/sqrt(RXm*(1-RXm)/40000)), breaks=150)

pdf("results.pdf", height=7.5, width=10)
par(las=1)
for(k in 1:12) {
  par(mfrow=c(2,2))
  for(s in 0:10) {
    plot(r, oRA[k, s+1,], ylim=c(0, 7/8), xlab="r", ylab="Pr(rec't haplotype)",
         main=paste("k = ", k, ", s = ", s, sep=""), cex=0.5)
    lines(r, RA[k,s+1,])
    points(r, oRXf[k, s+1,], col="red", cex=0.5)
    lines(r,  RXf[k, s+1,], col="red")
    points(r, oRXm[k, s+1,], col="blue", cex=0.5)
    lines(r, RXm[k, s+1,], col="blue")
  }
}
dev.off()

pdf("simulations.pdf", height=7.5, width=10)
par(las=1)
par(mfrow=c(2,2))
for(k in c(5,12)) {
  for(s in c(3,10)) {
    plot(r, oRA[k, s+1,], ylim=c(0, 7/8), xlab="recombination fraction (at meiosis)",
         ylab="Pr(recombinant haplotype)",
         main=paste("k = ", k, ", s = ", s, sep=""), cex=0.5, type="n")
    abline(lty=3, h=(0:7)/8, col="gray80", v=seq(0, 0.5, by=0.05))
    points(r, oRA[k, s+1,], cex=0.5)
    lines(r, RA[k,s+1,])
    points(r, oRXf[k, s+1,], col="red", cex=0.5)
    lines(r,  RXf[k, s+1,], col="red")
    points(r, oRXm[k, s+1,], col="blue", cex=0.5)
    lines(r, RXm[k, s+1,], col="blue")
  }
}
dev.off()
