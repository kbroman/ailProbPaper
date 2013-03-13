A <- cbind(c(0, 1-0.05, 0.05/64),
           c(1/2, (1-0.05)/2, 0.05/128),
           c(0, 0, 1))
pi0 <- c(0.03, 0.08, 1)
b <- cbind(c(1, 0, 0), c(0, 1, 0))
e <- eigen(A)
D <- diag(e$values)
v <- e$vectors
vi <- solve(v)
pi0 %*% v %*% D %*% vi %*% b

mat <- matrix(ncol=2, nrow=51)
mat[1,] <- pi0[1:2]
for(i in 1:50)
  mat[i+1,] <- pi0 %*% v %*% D^i %*% vi %*% b

plot(0:50, mat[,2], type="l", col="red")
lines(0:50, mat[,1], col="blue")
