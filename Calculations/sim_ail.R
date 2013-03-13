ailA <- simAIL(n.pairs=1000)
ailX <- simAIL(n.pairs=1000, xchr=TRUE)

r <- seq(0, 0.49, by=0.01)
loc <- imf.h(r)
ailA <- convertall2discrete(ailA, loc)
ailX <- convertall2discrete(ailX, loc)

g.ailA <- vector("list", length(ailA))
for(i in 1:length(ailA))
  g.ailA[[i]] <- gtable(lapply(ailA[[i]], function(a) a[[1]])) +
                      gtable(lapply(ailA[[i]], function(a) a[[2]]))

g.ailXf <- vector("list", length(ailX))
for(i in 1:length(ailX))
  g.ailXf[[i]] <- gtable(lapply(ailX[[i]], function(a) a[[1]]))

g.ailXm <- vector("list", length(ailX))
for(i in 1:length(ailX))
  g.ailXm[[i]] <- gtable(lapply(ailX[[i]], function(a) a[[1]]), "first")

######################################################################

system.time(ailXbal <- simAIL(n.pairs=1000, balanced=TRUE, xchr=TRUE))
system.time(ailXbal <- convertall2discrete(ailXbal, loc))

g.ailXfbal <- vector("list", length(ailXbal))
for(i in 1:length(ailXbal))
  g.ailXfbal[[i]] <- gtable(lapply(ailXbal[[i]], function(a) a[[1]]))

g.ailXmbal <- vector("list", length(ailXbal))
for(i in 1:length(ailXbal))
  g.ailXmbal[[i]] <- gtable(lapply(ailXbal[[i]], function(a) a[[1]]), "first")

######################################################################

R.ailA <- cbind(0, t(sapply(g.ailA, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
R.ailXfbal <- cbind(0, t(sapply(g.ailXfbal, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
R.ailXmbal <- cbind(0, t(sapply(g.ailXmbal, apply, 3, function(a) 1-sum(diag(a))/sum(a))))

f <- function(r, k) 
{
  if(length(r) > 1) return(sapply(r, f, k))
  A <- cbind(c(0, 1-r, r/4), c(1/2, (1-r)/2, r/8), c(0, 0, 1))
  pi0 <- cbind((1-r)/2, 1/4+(1-r)/4, 1)
  if(k==2) return(pi0[1:2])
  cur <- pi0
  for(j in 3:k) 
    cur <- cur %*% A
  cur[1:2]
}


k <- 4
plot(r, R.ailA[k-1,], lwd=2)
lines(r, 1-2*(1/4 + (1-2*r)*(1-r)^(k-2)/4), lwd=2)
points(r, R.ailXfbal[k-1,], lwd=2, col="red")
points(r, R.ailXmbal[k-1,], lwd=2, col="blue")
temp <- f(r, k)
lines(r, 1-2*temp[2,], lwd=2, col="red")
lines(r, 1-2*temp[1,], lwd=2, col="blue")

######################################################################

par(ask=TRUE)
for(k in 2:20) {
  plot(r, ailRprob(r, k=k, "autosome"), type="l", lwd=2, main=paste("k =", k),
       xlab="recombination fraction (at meiosis)",
       ylab="Pr(Recombinant haplotype)", las=1)
  lines(r, ailRprob(r, k=k, "femaleX"), lwd=2, col="red")
  lines(r, ailRprob(r, k=k, "maleX"), lwd=2, col="blue")
}

######################################################################
# comparisons to simulations
######################################################################
rm(g.ailA, g.ailXf, g.ailXfbal, g.ailXm, g.ailXmbal)

for(i in 13:20) {
  attachfile(i, "ail_results")
  if(i==13) {
    gailA <- g.ailA
    gailXf <- g.ailXf
    gailXfbal <- g.ailXfbal
    gailXm <- g.ailXm
    gailXmbal <- g.ailXmbal
  }
  else {
    for(j in seq(along=g.ailA)) {
      gailA[[j]] <- gailA[[j]] + g.ailA[[j]]
      gailXf[[j]] <- gailXf[[j]] + g.ailXf[[j]]
      gailXfbal[[j]] <- gailXfbal[[j]] + g.ailXfbal[[j]]
      gailXm[[j]] <- gailXm[[j]] + g.ailXm[[j]]
      gailXmbal[[j]] <- gailXmbal[[j]] + g.ailXmbal[[j]]
    }
  }
  detach(2)
}
oR.ailA <- cbind(0, t(sapply(gailA, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
oR.ailXfbal <- cbind(0, t(sapply(gailXfbal, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
oR.ailXmbal <- cbind(0, t(sapply(gailXmbal, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
oR.ailXf <- cbind(0, t(sapply(gailXf, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
oR.ailXm <- cbind(0, t(sapply(gailXm, apply, 3, function(a) 1-sum(diag(a))/sum(a))))

r <- seq(0, 0.49, by=0.01)

pdf("ail_results.pdf", height=7.5, width=10)
for(k in 2:12) {
  plot(r, oR.ailA[k-1,], lwd=2, main=paste("k =", k),
       xlab="recombination fraction (at meiosis)",
       ylab="Pr(recombinant haplotype)")
  lines(c(r, 0.5), ailRprob(c(r, 0.5), k), lwd=2)
  points(r, oR.ailXfbal[k-1,], lwd=2, col="red")
  lines(c(r, 0.5), ailRprob(c(r, 0.5), k, "femaleX"), lwd=2, col="red")
  points(r, oR.ailXmbal[k-1,], lwd=2, col="blue")
  lines(c(r, 0.5), ailRprob(c(r, 0.5), k, "maleX"), lwd=2, col="blue")

  points(r, oR.ailXf[k-1,], lwd=2, col="red", pch=4)
  points(r, oR.ailXm[k-1,], lwd=2, col="blue", pch=4)
}
abline(h=c(1/2,4/9), lty=2, col="gray")
dev.off()

######################################################################
# haplotype probabilities: unbalanced case
######################################################################
hapAAXf <- t(sapply(gailXf, apply, 3, function(a) a[1,1]/sum(a)))
hapBBXf <- t(sapply(gailXf, apply, 3, function(a) a[2,2]/sum(a)))
hapABXf <- t(sapply(gailXf, apply, 3, function(a) a[1,2]/sum(a)))
hapBAXf <- t(sapply(gailXf, apply, 3, function(a) a[2,1]/sum(a)))

hapAAXm <- t(sapply(gailXm, apply, 3, function(a) a[1,1]/sum(a)))
hapBBXm <- t(sapply(gailXm, apply, 3, function(a) a[2,2]/sum(a)))
hapABXm <- t(sapply(gailXm, apply, 3, function(a) a[1,2]/sum(a)))
hapBAXm <- t(sapply(gailXm, apply, 3, function(a) a[2,1]/sum(a)))

# single locus frequencies
plot(2:12, (hapAAXf + hapABXf)[,1], col="red", lwd=2, ylim=c(0.5, 0.75),
     xlab="Generation", ylab="Pr(A)", las=1)
lines(2:12, 2/3 + (1/3)*(-1/2)^(2:12), col="red", lwd=2)
points(2:12, (hapAAXm + hapABXm)[,1], col="blue", lwd=2)
lines(2:12, 2/3 + (1/3)*(-1/2)^(2:12-1), col="blue", lwd=2)

# haplotype frequencies
plot(2:12, hapAAXf[,1], lwd=2, ylim=c(0.25, 0.75), col="red")
lines(2:12, hapAAailub(0.01, 12)[-1,1], col="red", lwd=2,lty=1)
points(2:12, hapAAXm[,1],  lwd=2, pch=1, col="blue")
lines(2:12, hapAAailub(0.01, 12)[-1,2], col="blue", lwd=2, lty=1)
points(2:12, hapAAXf[,25], lwd=2, pch=4, col="red")
lines(2:12, hapAAailub(0.25, 12)[-1,1], col="red", lwd=2, lty=2)
points(2:12, hapAAXm[,25], lwd=2, pch=4, col="blue")
lines(2:12, hapAAailub(0.25, 12)[-1,2], col="blue", lwd=2, lty=2)
points(2:12, hapAAXf[,49], lwd=2, pch=0, col="red")
lines(2:12, hapAAailub(0.49, 12)[-1,1], col="red", lwd=2, lty=3)
points(2:12, hapAAXm[,49], lwd=2, pch=0, col="blue")
lines(2:12, hapAAailub(0.49, 12)[-1,2], col="blue", lwd=2, lty=3)


######################################################################
d <- imf.h(0.25)
f1 <- create.par(d, c(1,2))
n.pairs <- 1000
f2m <- f2f <- vector("list", n.pairs)
for(i in 1:n.pairs) {
  f2f[[i]] <- cross(f1, f1, xchr=TRUE, male=FALSE, obl=FALSE, m=0)
  f2m[[i]] <- cross(f1, f1, xchr=TRUE, male=TRUE, obl=FALSE, m=0)
}
f3m <- f3f <- vector("list", n.pairs)
for(i in 1:n.pairs) {
  f3f[[i]] <- cross(f2f[[i]], f2m[[i]], xchr=TRUE, male=FALSE, obl=FALSE, m=0)
  f3m[[i]] <- cross(f2f[[i]], f2m[[i]], xchr=TRUE, male=TRUE, obl=FALSE, m=0)
}
