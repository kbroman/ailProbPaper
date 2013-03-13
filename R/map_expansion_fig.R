source("map_expansion.R")

maxs <- 20
plot(1:maxs, doA(maxs, k=12), type="l", lwd=2,
     xlab="No. generations of random mating",
     ylab="Map expansion", ylim=c(0, 23.24),
     xaxs="i", yaxs="i", col="green")
lines(1:maxs, doA(maxs, k=6), lwd=2, col="red")
lines(1:maxs, hsA(maxs), lwd=2, col="blue")
lines(1:maxs, doX(maxs, 12), lty=2, lwd=2, col="green")
lines(1:maxs, doX(maxs, 6), lty=2, lwd=2, col="red")
lines(1:maxs, hsX(maxs), lty=2, lwd=2, col="blue")
lines(1:maxs, ailA(maxs), lwd=2, col="black")
lines(1:maxs, ailXbal(maxs), lty=2, lwd=2, col="black")
lines(1:maxs, ailXubal(maxs), lty=3, lwd=2, col="black")

######################################################################
# actual value in Svenson et al (2011)
######################################################################
dogen <- read.table("../Calculations/Fk_distribution.txt", header=TRUE, as.is=TRUE)
dogen[,1] <- as.numeric(substr(dogen[,1], 2, nchar(dogen[,1])))
dogen[,2] <- dogen[,2]/sum(dogen[,2])
