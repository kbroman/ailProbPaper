source("../Calculations/DOprob.R")
source("../Calculations/AILprob.R")
library(ricalc)

dogen <- read.table("../Calculations/Fk_distribution.txt", header=TRUE, as.is=TRUE)
dogen[,1] <- as.numeric(substr(dogen[,1], 2, nchar(dogen[,1])))
dogen[,2] <- dogen[,2]

# average k
avek <- sum(dogen[,1]*dogen[,2])/sum(dogen[,2]) # about 6

r <- seq(0, 0.5, by=0.01)

args=(commandArgs(TRUE))
eps <- FALSE
if(length(args)==0) {
  eps <- FALSE
} else {
  for(i in seq(along=args))
    eval(parse(text=args[[i]]))
}

if(eps) {
  postscript("../Figs/fig1.eps", height=5.5, width=6.5, pointsize=12, onefile=TRUE, horizontal=FALSE)
} else {
  pdf("../Figs/fig1.pdf", height=5.5, width=6.5, pointsize=12)
}
par(las=1, mar=c(4.6,4.1,1.1,1.1))
plot(r, doRprob(r, k=rep(dogen[,1], dogen[,2]), s=5, "autosome"), type="l",
     xaxs="i", yaxs="i", xlab="Recombination fraction (at meiosis)",
     ylab="Pr(recombinant haplotype)", lwd=2, ylim=c(0,1))
lines(r, doRprob(r, k=rep(dogen[,1], dogen[,2]), s=5, "female"), lwd=2, col="red")
lines(r, doRprob(r, k=rep(dogen[,1], dogen[,2]), s=5, "male"), lwd=2, col="blue")

lines(r, doRprob(r, k=1, s=10, "autosome"), lty=2, lwd=2)
lines(r, doRprob(r, k=1, s=10, "female"), lwd=2, col="red", lty=2)
lines(r, doRprob(r, k=1, s=10, "male"), lwd=2, col="blue", lty=2)

lines(r, ailRprob(r, 12), lwd=2, lty=3)
lines(r, ailRprob(r, 12, "female"), lwd=2, lty=3, col="red")
lines(r, ailRprob(r, 12, "male"  ), lwd=2, lty=3, col="blue")

# ail at F12
# HS at generation s=10
# DO as in Svenson et al (ave k = 6), s=5
legend("bottomright", lwd=rep(2,7), lty=c(1,2,3,1,1,1,2),
       col=c(rep("black", 4), "blue", "red","green"),
       c("DO at s=5", "HS at s=10", "AIL at s=12", "autosome", "female", "male", "Mott et al. (2000)"))

r2 <- c(r[-51], 0.495, 0.499)
d <- imf.gam(r2, 11.3)

lines(r2, 7*(1-exp(-(12)*d/100))/8, col="green", lwd=2, lty=2)

dev.off()
