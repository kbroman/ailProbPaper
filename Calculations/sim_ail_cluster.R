source("cross_func.R")
set.seed(18439413+SUB)
n.sim <- 4000
ailA <- simAIL(n.pairs=n.sim)

r <- seq(0, 0.49, by=0.01)
loc <- imf.h(r)
ailA <- convertall2discrete(ailA, loc)

g.ailA <- vector("list", length(ailA))
for(i in 1:length(ailA))
  g.ailA[[i]] <- gtable(lapply(ailA[[i]], function(a) a[[1]])) +
                      gtable(lapply(ailA[[i]], function(a) a[[2]]))

save(g.ailA, file="ail_resultsSUB.RData")

ailX <- simAIL(n.pairs=n.sim, xchr=TRUE)
ailX <- convertall2discrete(ailX, loc)

g.ailXf <- vector("list", length(ailX))
for(i in 1:length(ailX))
  g.ailXf[[i]] <- gtable(lapply(ailX[[i]], function(a) a[[1]]))

g.ailXm <- vector("list", length(ailX))
for(i in 1:length(ailX))
  g.ailXm[[i]] <- gtable(lapply(ailX[[i]], function(a) a[[1]]), "first")

save(g.ailA, g.ailXf, g.ailXm, file="ail_resultsSUB.RData")

######################################################################

ailXbal <- simAIL(n.pairs=n.sim, balanced=TRUE, xchr=TRUE)
system.time(ailXbal <- convertall2discrete(ailXbal, loc))

g.ailXfbal <- vector("list", length(ailXbal))
for(i in 1:length(ailXbal))
  g.ailXfbal[[i]] <- gtable(lapply(ailXbal[[i]], function(a) a[[1]]))

g.ailXmbal <- vector("list", length(ailXbal))
for(i in 1:length(ailXbal))
  g.ailXmbal[[i]] <- gtable(lapply(ailXbal[[i]], function(a) a[[1]]), "first")

save(g.ailA, g.ailXf, g.ailXm, g.ailXfbal, g.ailXmbal, file="ail_resultsSUB.RData")

