
ailA <- function(maxs) c(0, (2:maxs)/2)

ailXbal <- function(maxs) ailA(maxs)*2/3

ailXubal <-
function(maxs)
{  
  if(maxs < 2) stop("maxs >= 2")
  me <- rep(NA, maxs)
  me[1] <- 0
  me[2] <- 2/3

  freqAfem <- 2/3 + 1/3 * (-1/2)^(1:maxs)
  freqAmal <- 2/3 + 1/3 * (-1/2)^(1:maxs-1)
  
  for(s in 3:maxs) 
    me[s] <- me[s-1] + 4/3*(freqAfem[s-1] - freqAfem[s-2]*freqAmal[s-2])
  me
}

hsA <- function(maxs) (7*(1:maxs)+17)/8

hsX <- function(maxs) hsA(maxs)*2/3

doA <-
function(maxs, k)
{
  s <- 1:maxs
  (7*s+49)/8 - (15+7*sqrt(5))/5*((1+sqrt(5))/4)^(k+1) - (15-7*sqrt(5))/5*((1-sqrt(5))/4)^(k+1)
}

doX <- function(maxs, k) doA(maxs, k)*2/3
  
