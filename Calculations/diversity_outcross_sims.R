set.seed(27729779+SUB)
source("cross_func.R")

doA <- simDO(n.pairs=40000, k=SUB) # 22 sec
doX <- simDO(n.pairs=40000, xchr=TRUE, k=SUB) # 10 sec

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


save(gA, gXf, gXm, file="diversity_outcross_resultsSUB.RData")
