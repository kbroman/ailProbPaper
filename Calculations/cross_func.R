######################################################################
#
# func.R
#
# Karl W Broman, Johns Hopkins University
# 18 Nov 2003
#
# functions for simulating meiotic products and RI lines
#
######################################################################

##############################
# meiosis.sub <- 
#
# Simulate the locations of crossovers on a meiotic product
# via the chi-square model (m=0 corresponds to no interference)
# and potentially with an obligate chiasma.
##############################

meiosis.sub <-
function(L, m=10, obligate.chiasma=TRUE)
{
  if(obligate.chiasma) { # adjust mean no. chiasmata
    if(L <= 50) stop("L must be > 50 cM")
    if(m==0) f <- function(Lstar,f.L,f.m) f.L-Lstar/(1-exp(-Lstar/50)) 
    else {
      f <- function(Lstar,f.L,f.m=0)
        {
          lambdastar <- Lstar/50*(f.m+1)
          temp <- lambdastar
          for(i in 1:length(temp))
            temp[i] <- sum(exp(-lambdastar[i] + (0:f.m)*log(lambdastar[i])-
                               lgamma((0:f.m)+1)) * (f.m+1-(0:f.m))/(f.m+1))
          f.L - Lstar/(1-temp)
        }
    }

    Lstar <- uniroot(f,c(1e-5,L+1),f.m=m,f.L=L)$root
  }
  else Lstar <- L

  if(m==0) { # no interference
    if(!obligate.chiasma)  # no obligate chiasma
      n.xo <- rpois(1,Lstar/100)
    else {
      up <- qpois(1e-14,Lstar/50,lower=FALSE)
      p <- dpois(1:up,Lstar/50)/ppois(0,Lstar/50)
      n.chi <- sample(1:up,1,prob=p)
      n.xo <- rbinom(1,n.chi,0.5)
    }
    if(n.xo==0) xo <- NULL
    else xo <- sort(runif(n.xo,0,L))
  }
  else { # chi-square model
    n.chi <- 0
    while(n.chi == 0) {
      n.pts <- rpois(1,Lstar/50*(m+1))
      first <- sample(1:(m+1),1)
      if(first <= n.pts || !obligate.chiasma) n.chi <- 1
    }
    if(first > n.pts)
      xo <- NULL
    else {
      pt.loc <- sort(runif(n.pts,0,L))
      chi.loc <- pt.loc[seq(first,length(pt.loc),by=m+1)]
      n.xo <- rbinom(1,length(chi.loc),0.5)
      if(n.xo==0) xo <- NULL
      else if(length(chi.loc)==1) xo <- chi.loc
      else xo <- sort(sample(chi.loc,n.xo,repl=FALSE))
    }
  }
    
  if(length(xo) == 0) xo <- NULL
  xo
}

##############################
# create.par
#
# create a parental individual
##############################

create.par <-
function(L, allele=1)
{
  if(length(allele) == 1) allele <- rep(allele,2)
  if(length(allele) != 2)
    stop("allele should be of length 1 or 2")
  
  list(mat=rbind(c(0,L),allele[1]),
       pat=rbind(c(0,L),allele[2]))
}


##############################
# meiosis
#
# Output a random meiotic product from an
# input individual.
##############################

meiosis <-
function(parent, m=10, obligate.chiasma=TRUE)
{
  L <- parent$mat[1,ncol(parent$mat)]
  if(abs(parent$pat[1,ncol(parent$pat)] - L) > 1e-13)
    stop("There is a problem with the parent's data structure.")

  product <- meiosis.sub(L, m, obligate.chiasma)
  a <- sample(1:2,1)
#  print(product)
#  print(a)

  if(length(product)==0) return(parent[[a]])

  else {
    for(i in 1:length(product)) {
#      cat("\t",i,a,"\n")
      if(i == 1) 
        result <- parent[[a]][,parent[[a]][1,]<product[1],drop=FALSE]
      else {
        temp <- parent[[a]][1,]>=product[i-1] & parent[[a]][1,]<product[i]
        result <- cbind(result,parent[[a]][,temp])
      }
      u <- parent[[a]][2,parent[[a]][1,]>=product[i]]
      result <- cbind(result,c(product[i],u[1]))
      a <- 3-a
    }
    temp <- parent[[a]][1,]>=product[length(product)]
    result <- cbind(result,parent[[a]][,temp])
  }

  # clean out excess stuff in the result
  if(ncol(result)>2) {
    keep <- rep(TRUE,ncol(result))
    for(i in 2:(ncol(result)-1)) 
      if(result[2,i] == result[2,i+1])
        keep[i] <- FALSE
  }
#  print(result)
  result[,keep,drop=FALSE]
}


##############################
# cross
#
# cross two individuals to create a
# single progeny
##############################

cross <-
function(mom, dad, m=10, obligate.chiasma=TRUE, xchr=FALSE, male=FALSE)
{
  if(!xchr) {
    return(list(mat=meiosis(mom,m,obligate.chiasma),
                pat=meiosis(dad,m,obligate.chiasma)))
  }
  else {
    if(male)
      return(list(mat=meiosis(mom,m,obligate.chiasma),
                  pat=dad$pat))
    else
      return(list(mat=meiosis(mom,m,obligate.chiasma),
                  pat=dad$mat))
  }
}

ri2 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE)
{
  f1 <- create.par(L,c(1,2))
  par1 <- cross(f1,f1,m,obligate.chiasma)
  par2 <- cross(f1,f1,m,obligate.chiasma)
  for(i in 1:n.gen) {
    c1 <- cross(par1,par2,m,obligate.chiasma)
    c2 <- cross(par1,par2,m,obligate.chiasma)
    par1 <- c1
    par2 <- c2
  }
  par1
}

ri8 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE)
{
  f1a <- create.par(L,c(1,2))
  f1b <- create.par(L,c(3,4))
  f1c <- create.par(L,c(5,6))
  f1d <- create.par(L,c(7,8))
  par1 <- cross(f1a,f1b,m,obligate.chiasma)
  par2 <- cross(f1c,f1d,m,obligate.chiasma)
  if(length(n.gen)==1) {
    for(i in 1:(n.gen+1)) {
      c1 <- cross(par1,par2,m,obligate.chiasma)
      c2 <- cross(par1,par2,m,obligate.chiasma)
      par1 <- c1
      par2 <- c2
    }
    return(par1)
  }
  else {
    result <- vector("list",length(n.gen))
    names(result) <- n.gen
    n.gen <- c(-1,n.gen)
    for(j in 2:length(n.gen)) {
      for(i in (n.gen[j-1]+2):(n.gen[j]+1)) {
        c1 <- cross(par1,par2,m,obligate.chiasma)
        c2 <- cross(par1,par2,m,obligate.chiasma)
        par1 <- c1
        par2 <- c2
      }
      result[[j-1]] <- par1
    }
    return(result)
  }
        
}
    
    
##############################
# where.het
#
# find regions of heterozygosity
# in an individual
##############################
where.het <-
function(ind)
{
  if(ncol(ind$mat)==ncol(ind$pat) && all(ind$mat == ind$pat)) {
#    cat(" --No regions of heterozygosity\n")
    return(NULL)
  }
  u <- sort(unique(c(ind$mat[1,],ind$pat[1,])))
  het <- NULL
  for(i in 2:length(u)) {
    mat <- ind$mat[,ind$mat[1,] >= u[i],drop=FALSE]
    mat <- mat[2,1]

    pat <- ind$pat[,ind$pat[1,] >= u[i],drop=FALSE]
    pat <- pat[2,1]

    if(mat!=pat) { # heterozygous
      if(is.null(het)) het <- cbind(u[i-1],u[i])
      else het <- rbind(het,c(u[i-1],u[i]))
    }
  }

  # clean up
  if(nrow(het) > 1) {
    keep <- rep(TRUE,nrow(het))
    for(j in 2:nrow(het)) {
      if(het[j,1] == het[j-1,2]) {
        het[j,1] <- het[j-1,1]
        keep[j-1] <- FALSE
      }
    }
    het <- het[keep,,drop=FALSE]
  }
  het
}

simDO <-
function(n.pairs=500, k=6, maxs=10, L=200, xchr=FALSE, m=0, obligate.chiasma=FALSE)
{
  if(k < 1) stop("k must be >= 1")
  if(n.pairs < 2) stop("n.pairs must be >= 2")
  if(maxs < 1) stop("maxs must be >= 1")
  if(L < 0) stop("L must be > 0")

  nalle <- 8
  res <- vector("list", 1+maxs)
  for(i in 1:(maxs+1))
    res[[i]] <- vector("list", n.pairs)

  for(i in 1:n.pairs) {
    # simulate F0 of collaborate cross
    crossdir <- sample(1:nalle)
    v <- vector("list", nalle/2)
    mom <- cross(create.par(L, crossdir[1:2]),
                 create.par(L, crossdir[3:4]), xchr=xchr, male=FALSE,
                 m=m, obligate.chiasma=obligate.chiasma)
    dad <- cross(create.par(L, crossdir[5:6]),
                 create.par(L, crossdir[7:8]), xchr=xchr, male=TRUE,
                 m=m, obligate.chiasma=obligate.chiasma)

    # simulate k generations of inbreeding
    for(j in 1:k) {
      fkid <- cross(mom, dad, xchr=xchr, male=FALSE,
                 m=m, obligate.chiasma=obligate.chiasma)
      mkid <- cross(mom, dad, xchr=xchr, male=TRUE,
                 m=m, obligate.chiasma=obligate.chiasma)
      mom <- fkid
      dad <- mkid
    }
    res[[1]][[i]] <- list(fkid, mkid)
  }

  # simulate DO
  for(i in 2:(maxs+1)) {
    # randomly swap dads
    z <- sample(n.pairs)
    while(any(z==1:n.pairs))
      z <- sample(n.pairs)

    # generate two offspring per pair
    for(j in 1:n.pairs)
      res[[i]][[j]] <- list(cross(res[[i-1]][[j]][[1]], res[[i-1]][[z[j]]][[2]], xchr=xchr, male=FALSE,
                                  m=m, obligate.chiasma=obligate.chiasma),
                            cross(res[[i-1]][[j]][[1]], res[[i-1]][[z[j]]][[2]], xchr=xchr, male=TRUE,
                                  m=m, obligate.chiasma=obligate.chiasma))
  }
  res
}

convert2discrete <-
function(individual, locations) {
  f <-
    function(x,y) {
      if(!any(y[1,]<x)) return(y[2,1])
      y[2,max(which(y[1,] < x))+1]
    }
  res <- t(sapply(individual, function(a,b) sapply(b, f, a), locations))
  colnames(res) <- myround(locations, 2)
  res
}

convertall2discrete <-
function(x, locations)
{
  if(!is.matrix(x[[1]][[1]]))
    return(lapply(x, convertall2discrete, locations))
  lapply(x, convert2discrete, locations)
}

gtable <-
function(ind, what=c("both", "first"))
{
  what <- match.arg(what)
  nloc <- ncol(ind[[1]])
  maxa <- max(unlist(ind))
  g <- array(dim=c(maxa, maxa, nloc-1))
  for(i in 1:(nloc-1)) {
    g[,,i] <- table(factor(sapply(ind, function(a) a[1,1]), levels=1:maxa),
                    factor(sapply(ind, function(a) a[1,i+1]), levels=1:maxa))
    if(what=="both")
      g[,,i] <- g[,,i] + table(factor(sapply(ind, function(a) a[2,1]), levels=1:maxa),
                               factor(sapply(ind, function(a) a[2,i+1]), levels=1:maxa))
  }
  g
}
                    
    

simAIL <-
function(n.pairs=500, maxk=12, L=200, xchr=FALSE, m=0, obligate.chiasma=FALSE,
         balanced=FALSE)
{
  if(maxk < 2) stop("maxk must be >= 2")
  if(n.pairs < 2) stop("n.pairs must be >= 2")
  if(L < 0) stop("L must be > 0")

  res <- vector("list", maxk-1)
  for(i in 1:(maxk-1))
    res[[i]] <- vector("list", n.pairs)

  # create F1 
  mom <- create.par(L, 1:2)
  dad <- create.par(L, 1:2)
  dad2 <- create.par(L, 2:1)

  # generate F2
  for(i in 1:n.pairs) {
    if(!balanced || i <=n.pairs/2) {
      fkid <- cross(mom, dad, xchr=xchr, male=FALSE,
                    m=m, obligate.chiasma=obligate.chiasma)
      mkid <- cross(mom, dad, xchr=xchr, male=TRUE,
                    m=m, obligate.chiasma=obligate.chiasma)
    }
    else {
      fkid <- cross(mom, dad2, xchr=xchr, male=FALSE,
                    m=m, obligate.chiasma=obligate.chiasma)
      mkid <- cross(mom, dad2, xchr=xchr, male=TRUE,
                    m=m, obligate.chiasma=obligate.chiasma)
    }
    res[[1]][[i]] <- list(fkid, mkid)
  }
  if(maxk==2) return(res[[1]])

  for(i in 2:(maxk-1)) {
    # randomly swap dads
    z <- sample(n.pairs)
    while(any(z==1:n.pairs))
      z <- sample(n.pairs)

    # generate two offspring per pair
    for(j in 1:n.pairs)
      res[[i]][[j]] <- list(cross(res[[i-1]][[j]][[1]], res[[i-1]][[z[j]]][[2]], xchr=xchr, male=FALSE,
                                  m=m, obligate.chiasma=obligate.chiasma),
                            cross(res[[i-1]][[j]][[1]], res[[i-1]][[z[j]]][[2]], xchr=xchr, male=TRUE,
                                  m=m, obligate.chiasma=obligate.chiasma))
  }
  res
}
