source("cross_func.R")
source("DOprob.R")

doA <- simDO(n.pairs=500) # 22 sec
doX <- simDO(n.pairs=500, xchr=TRUE) # 10 sec

r <- seq(0, 0.49, by=0.01)
loc <- imf.h(r)
doA <- convertall2discrete(doA, loc) # 15 sec
doX <- convertall2discrete(doX, loc) # 15 sec

gA <- vector("list", length(doA))
for(i in 1:length(doA))
  g[[i]] <- gtable(lapply(doA[[i]], function(a) a[[1]])) +
               gtable(lapply(doA[[i]], function(a) a[[2]]))

gXf <- vector("list", length(doX))
for(i in 1:length(doX))
  gXf[[i]] <- gtable(lapply(doX[[i]], function(a) a[[1]]))

gXm <- vector("list", length(doX))
for(i in 1:length(doX))
  gXm[[i]] <- gtable(lapply(doX[[i]], function(a) a[[1]]), "first")


RA  <- cbind(0, t(sapply(gA,  apply, 3, function(a) 1-sum(diag(a))/sum(a))))
RXf <- cbind(0, t(sapply(gXf, apply, 3, function(a) 1-sum(diag(a))/sum(a))))
RXm <- cbind(0, t(sapply(gXm, apply, 3, function(a) 1-sum(diag(a))/sum(a))))


plot(r, RA[1,], type="l", ylim=c(0,1))
for(i in 2:length(x)) lines(r, RA[i,])

plot(r, RA[1,], type="l", ylim=c(0,1))
for(i in 2:length(x)) lines(r, RA[i,])
for(i in 1:length(x)) lines(r, RXf[i,], col="red")
for(i in 1:length(x)) lines(r, RXm[i,], col="blue")


######################################################################

n.sim <- 1000
v <- vector("list", n.sim)
for(i in 1:n.sim)
  v[[i]] <- ri8(L=200, n.gen=50, m=0, obligate.chiasma=FALSE)

r <- seq(0, 0.49, by=0.01)
loc <- imf.h(r)
vd <- convertall2discrete(v, loc)

g <- gtable(vd)
R <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))

fsnp <- simFounderSnps(list(a=loc), n.str="8")
x <- sim.cross(list(a=loc), n.ind=n.sim, type="ri8sib", foun=fsnp)
x <- x$geno[[1]]$truegeno
v <- vector("list", n.sim)
for(i in 1:n.sim) v[[i]] <- rbind(x[i,], x[i,])
g <- gtable(v)
R <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))
plot(r[-1], R)
lines(r, 7*r/(1+6*r))

fsnp4 <- simFounderSnps(list(a=loc), n.str="4")
x <- sim.cross(list(a=loc), n.ind=n.sim, type="ri4sib", foun=fsnp4)
x <- x$geno[[1]]$truegeno
v <- vector("list", n.sim)
for(i in 1:n.sim) v[[i]] <- rbind(x[i,], x[i,])
g <- gtable(v)
R <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))
points(r[-1], R, col="green")
lines(r, 6*r/(1+6*r), col="green")


n.sim <- 1000
v <- vector("list", n.sim)
for(i in 1:n.sim)
  v[[i]] <- ri8(L=200, n.gen=50, m=0, obligate.chiasma=FALSE)
vd <- convertall2discrete(v, loc)
g <- gtable(vd)
R <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))

doA <- simDO()
doA <- convertall2discrete(doA, loc)
R <- matrix(ncol=length(loc)-1, nrow=length(doA))
for(i in 1:length(doA)) {
  cat(i,"\n")
  g <- gtable(lapply(doA[[i]], function(a) a[[1]])) +
    gtable(lapply(doA[[i]], function(a) a[[2]]))
  R[i,] <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))
}

plot(r[-1], R[1,], type="l", ylim=c(0,1))
for(i in 2:length(doA)) lines(r[-1], R[i,])

source("DOprob.R")


doA <- simDO(maxs=10, k=1)
doA <- convertall2discrete(doA, loc)
R <- matrix(ncol=length(loc)-1, nrow=length(doA))
for(i in 1:length(doA)) {
  cat(i,"\n")
  g <- gtable(lapply(doA[[i]], function(a) a[[1]])) +
    gtable(lapply(doA[[i]], function(a) a[[2]]))
  R[i,] <- apply(g, 3, function(a) 1-sum(diag(a))/sum(a))
}
plot(r[-1], R[1,])
lines(r, rikRprob(r, 1))

par(mfrow=c(3,3))
for(i in 1:9) {
  plot(r[-1], R[i+1,])
  lines(r, doRprob(r, k=1, s=i-1))
}

par(mfrow=c(3,3))
for(i in 1:9) {
  plot(r[-1], R[i+1,])
  lines(r, doRprob(r, k=0, s=i+1))
}

par(mfrow=c(3,3))
for(i in 1:9) {
  plot(r[-1], R[i+1,])
  lines(r, doRprob(r, k=0, s=i))
}
