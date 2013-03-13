source("map_expansion.R")

y <- ailXubal(30)

y3 <- y[-(1:2)] - y[2]
x3 <- 1:28
summary(lm(y3 ~ -1 + x3))

y4 <- y[-(1:3)] - y[3]
x4 <- 1:27
summary(o <- lm(y4 ~ -1 + x4))
plot(x4, y4); abline(0, o$coef)
plot(o$fitted, o$resid/o$fitted*100)

y5 <- y[-(1:4)] - y[4]
x5 <- 1:26
summary(o <- lm(y5 ~ -1 + x5))
plot(x5, y5); abline(0, o$coef)
plot(o$fitted, o$resid/o$fitted*100)

