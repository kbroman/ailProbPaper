source("rifig_func.R")
png("broman_g3_thumb_big.png", height=1500, width=1500, res=288)

color <- c(rgb(240,240,  0,maxColorValue=255),
           rgb(128,128,128,maxColorValue=255),
           rgb(240,128,128,maxColorValue=255),
           rgb( 16, 16,240,maxColorValue=255),
           rgb(  0,160,240,maxColorValue=255),
           rgb(  0,160,  0,maxColorValue=255),
           rgb(240,  0,  0,maxColorValue=255),
           rgb(144,  0,224,maxColorValue=255))

file <- "ail_sim.RData"
if(!file.exists(file)) {
  f1 <- create.par(100,c(1,2))
  set.seed(9951035+3)
  f2 <- vector("list",200)
  for(i in 1:200) f2[[i]] <- cross(f1,f1,m=10,obl=TRUE)
  f3 <- vector("list",200)
  for(i in 1:100) {
    f3[[i*2-1]] <- cross(f2[[i*2-1]],f2[[i*2]],m=10,obl=TRUE)
    f3[[i*2]] <- cross(f2[[i*2-1]],f2[[i*2]],m=10,obl=TRUE)
  }
  f3p <- f3
  f3p[-(1:10)] <- f3p[sample(11:200)]
  mate3 <- c(1,3,2,5,4,7,6,9,8,10)
  f3p[1:10] <- f3p[mate3]
  f4 <- vector("list",200)
  for(i in 1:100) {
    f4[[i*2-1]] <- cross(f3p[[i*2-1]],f3p[[i*2]],m=10,obl=TRUE)
    f4[[i*2]] <- cross(f3p[[i*2-1]],f3p[[i*2]],m=10,obl=TRUE)
  }
  f4p <- f4
  f4p[-(1:10)] <- f4p[sample(11:200)]
  mate4 <- c(1,7,3,5,2,9,6,8,4,10)
  f4p[1:10] <- f4p[mate4]
  f5 <- vector("list",200)
  for(i in 1:100) {
    f5[[i*2-1]] <- cross(f4p[[i*2-1]],f4p[[i*2]],m=10,obl=TRUE)
    f5[[i*2]] <- cross(f4p[[i*2-1]],f4p[[i*2]],m=10,obl=TRUE)
  }

  temp <- f5
  fk <- vector("list", 200)
  for(j in 1:5) {
    temp <- temp[sample(length(temp))]
    for(i in 1:100) {
      fk[[i*2-1]] <- cross(temp[[i*2-1]],temp[[i*2]],m=10,obl=TRUE)
      fk[[i*2]] <- cross(temp[[i*2-1]],temp[[i*2]],m=10,obl=TRUE)
    }
    temp <- fk
  }
  save(f1,f2,f3,f4,f5,fk,mate3,mate4,file=file)
} else load(file)


par(mar=rep(0,4))
plot(0,0,xlim=c(-50,1097),ylim=c(25-80,485),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
u <- par("usr")
abline(v=u[1:2], h=u[3:4], lwd=6, lend=1, ljoin=1)
u[1] <- u[1]+15


rect(c(300,326)+114,c(480,480),c(312,338)+114,c(440,440),col=color[1],border=color[1], lend=1, ljoin=1)
rect(c(526,552)+114,c(480,480),c(538,564)+114,c(440,440),col=color[5],border=color[5], lend=1, ljoin=1)
rect(c(300,326)+114,c(480,480),c(312,338)+114,c(440,440),lend=1, ljoin=1)
rect(c(526,552)+114,c(480,480),c(538,564)+114,c(440,440),lend=1, ljoin=1)
points(432+114,465,pch=4,cex=0.7)
segments(432+114,450,432+114,430)
segments(319+114,430,545+114,430)
arrows(c(319,545)+114,c(430,430),c(319,545)+114,c(410,410),len=0.04)
text(300-25+114,460,expression(A),adj=c(1,0.5),cex=1.0)
text(564+25+114,460,expression(B),adj=c(0,0.5),cex=1.0)
text(u[1]/1.2,460,expression(F[0]),adj=c(0,0.5),cex=1.0)

rect(300+114,400,312+114,360,col=color[1],border=color[1], lend=1, ljoin=1)
rect(526+114,400,538+114,360,col=color[1],border=color[1], lend=1, ljoin=1)
rect(326+114,400,338+114,360,col=color[5],border=color[5], lend=1, ljoin=1)
rect(552+114,400,564+114,360,col=color[5],border=color[5], lend=1, ljoin=1)
rect(300+114,400,312+114,360,lend=1, ljoin=1)
rect(526+114,400,538+114,360,lend=1, ljoin=1)
rect(326+114,400,338+114,360,lend=1, ljoin=1)
rect(552+114,400,564+114,360,lend=1, ljoin=1)
text(u[1]/1.2,380,expression(F[1]),adj=c(0,0.5),cex=1.0)
points(432+114,385,pch=4,cex=0.7)
segments(432+114,370,432+114,350)


xloc <- 38
mult <- 40/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
xxloc <- NULL
for(i in 1:5) {
  rect(xloc-6,320,xloc+6,280,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+20,320,xloc+32,280,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+76,320,xloc+88,280,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+102,320,xloc+114,280,col=color[1],border=color[1], lend=1, ljoin=1)

  f2m <- f2[[2*i-1]]$mat
  for(j in 2:ncol(f2m)) {
    if(f2m[2,j]==2)
      rect(xloc-6,280+f2m[1,j]*mult,xloc+6,280+f2m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f2p <- f2[[2*i-1]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+20,280+f2p[1,j]*mult,xloc+32,280+f2p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  f2m <- f2[[2*i]]$mat
  for(j in 2:ncol(f2m)) {
    if(f2m[2,j]==2)
      rect(xloc+76,280+f2m[1,j]*mult,xloc+88,280+f2m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f2p <- f2[[2*i]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+102,280+f2p[1,j]*mult,xloc+114,280+f2p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  rect(xloc-6,320,xloc+6,280,    lend=1, ljoin=1)
  rect(xloc+20,320,xloc+32,280,  lend=1, ljoin=1)
  rect(xloc+76,320,xloc+88,280,  lend=1, ljoin=1)
  rect(xloc+102,320,xloc+114,280,lend=1, ljoin=1)

  xxloc <- c(xxloc,xloc+(6+20)/2,xloc+(88+102)/2)

  points(xloc+54,300,pch=4,cex=0.7)
  segments(xloc+54,290,xloc+54,270)
  segments(xxloc[2*i-1],270,xxloc[2*i],270)
  
  xloc <- xloc+78+120+30
}

segments(min(xxloc),350,max(xxloc),350)
arrows(xxloc,c(350,350),xxloc,c(330,330),len=0.04)
text(u[1]/1.2,300,expression(F[2]),adj=c(0,0.5),cex=1.0)
arrows(xxloc,270,xxloc,250,len=0.04)

xloc <- 38
z <- NULL
for(i in 1:5) {
  rect(xloc-6,240,xloc+6,200,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+20,240,xloc+32,200,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+76,240,xloc+88,200,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+102,240,xloc+114,200,col=color[1],border=color[1], lend=1, ljoin=1)

  f3m <- f3[[2*i-1]]$mat
  for(j in 2:ncol(f3m)) {
    if(f3m[2,j]==2)
      rect(xloc-6,200+f3m[1,j]*mult,xloc+6,200+f3m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f3p <- f3[[2*i-1]]$pat
  for(j in 2:ncol(f3p)) {
    if(f3p[2,j]==2)
      rect(xloc+20,200+f3p[1,j]*mult,xloc+32,200+f3p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  f3m <- f3[[2*i]]$mat
  for(j in 2:ncol(f3m)) {
    if(f3m[2,j]==2)
      rect(xloc+76,200+f3m[1,j]*mult,xloc+88,200+f3m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f3p <- f3[[2*i]]$pat
  for(j in 2:ncol(f3p)) {
    if(f3p[2,j]==2)
      rect(xloc+102,200+f3p[1,j]*mult,xloc+114,200+f3p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  rect(xloc-6,240,xloc+6,200,    lend=1, ljoin=1)
  rect(xloc+20,240,xloc+32,200,  lend=1, ljoin=1)
  rect(xloc+76,240,xloc+88,200,  lend=1, ljoin=1)
  rect(xloc+102,240,xloc+114,200,lend=1, ljoin=1)

#  points(xloc+54,220,pch=4,cex=0.7)
#  segments(xloc+54,210,xloc+54,190)
#  segments(xxloc[2*i-1],190,xxloc[2*i],190)
  
  z <- c(z, xloc+13, xloc+95)

  xloc <- xloc+78+120+30
}

#arrows(z, 195, z, 165, col="red", lwd=1, lend=1, ljoin=1, len=0.04)
#segments(z[1], 200, z[10], 200)
#segments(z[1], 160, z[10], 160)

offset <- c(0, 3, 0, 3, 0)
col <- c("green", "blue", "orange", "red", "black")

for(i in 1:5) {
  mom <- mate3[i*2-1]
  dad <- mate3[i*2]
  kids <- (i*2-1):(i*2)
  segments(z[mom], 195+5, z[mom], 185+offset[i], lend=1, ljoin=1, col=col[i])
  segments(z[dad], 195+5, z[dad], 185+offset[i], lend=1, ljoin=1, col=col[i])
  segments(z[mom], 185+offset[i], z[dad], 185+offset[i], lend=1, ljoin=1, col=col[i])
  segments(mean(z[kids]), 185+offset[i], mean(z[kids]), 180, lend=1, ljoin=1, col=col[i])
  segments(z[kids[1]], 180, z[kids[2]], 180, lend=1, ljoin=1, col=col[i])
  arrows(z[kids], 180, z[kids], 165, len=0.04, lend=1, ljoin=1, col=col[i])
}

text(u[1]/1.2,220,expression(F[3]),adj=c(0,0.5),cex=1.0)

xloc <- 38
for(i in 1:5) {
  rect(xloc-6,160,xloc+6,120,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+20,160,xloc+32,120,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+76,160,xloc+88,120,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+102,160,xloc+114,120,col=color[1],border=color[1], lend=1, ljoin=1)

  f4m <- f4[[2*i-1]]$mat
  for(j in 2:ncol(f4m)) {
    if(f4m[2,j]==2)
      rect(xloc-6,120+f4m[1,j]*mult,xloc+6,120+f4m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f4p <- f4[[2*i-1]]$pat
  for(j in 2:ncol(f4p)) {
    if(f4p[2,j]==2)
      rect(xloc+20,120+f4p[1,j]*mult,xloc+32,120+f4p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  f4m <- f4[[2*i]]$mat
  for(j in 2:ncol(f4m)) {
    if(f4m[2,j]==2)
      rect(xloc+76,120+f4m[1,j]*mult,xloc+88,120+f4m[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f4p <- f4[[2*i]]$pat
  for(j in 2:ncol(f4p)) {
    if(f4p[2,j]==2)
      rect(xloc+102,120+f4p[1,j]*mult,xloc+114,120+f4p[1,j-1]*mult,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  rect(xloc-6,160,xloc+6,120,    lend=1, ljoin=1)
  rect(xloc+20,160,xloc+32,120,  lend=1, ljoin=1)
  rect(xloc+76,160,xloc+88,120,  lend=1, ljoin=1)
  rect(xloc+102,160,xloc+114,120,lend=1, ljoin=1)

#  points(xloc+54,140,pch=4,cex=0.7)
#  arrows(xloc+54,87,xloc+54,80,len=0.04)
#  arrows(xloc+54,130,xloc+54,80,len=0.04,lty=3)
  
  xloc <- xloc+78+120+30
}

offset <- c(0, 2, 4, 2, 6)
col <- c("green", "blue", "orange", "red", "black")

for(i in 1:5) {
  mom <- mate4[i*2-1]
  dad <- mate4[i*2]
  kids <- (i*2-1):(i*2)
  segments(z[mom], 195-80+5, z[mom], 185-80+offset[i], lend=1, ljoin=1, col=col[i])
  segments(z[dad], 195-80+5, z[dad], 185-80+offset[i], lend=1, ljoin=1, col=col[i])
  segments(z[mom], 185+offset[i]-80, z[dad], 185-80+offset[i], lend=1, ljoin=1, col=col[i])
  segments(mean(z[kids]), 185+offset[i]-80, mean(z[kids]), 180-80, lend=1, ljoin=1, col=col[i])
  segments(z[kids[1]], 180-80, z[kids[2]], 180-80, lend=1, ljoin=1, col=col[i])
  arrows(z[kids], 180-80, z[kids], 165-80, len=0.04, lend=1, ljoin=1, col=col[i])
}

text(u[1]/1.2,140,expression(F[4]),adj=c(0,0.5),cex=1.0)

xloc <- 38
for(i in 1:5) {
  rect(xloc-6,160-80,xloc+6,120-80,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+20,160-80,xloc+32,120-80,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+76,160-80,xloc+88,120-80,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+102,160-80,xloc+114,120-80,col=color[1],border=color[1], lend=1, ljoin=1)

  f5m <- f5[[2*i-1]]$mat
  for(j in 2:ncol(f5m)) {
    if(f5m[2,j]==2)
      rect(xloc-6,120+f5m[1,j]*mult-80,xloc+6,120+f5m[1,j-1]*mult-80,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f5p <- f5[[2*i-1]]$pat
  for(j in 2:ncol(f5p)) {
    if(f5p[2,j]==2)
      rect(xloc+20,120+f5p[1,j]*mult-80,xloc+32,120+f5p[1,j-1]*mult-80,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  f5m <- f5[[2*i]]$mat
  for(j in 2:ncol(f5m)) {
    if(f5m[2,j]==2)
      rect(xloc+76,120+f5m[1,j]*mult-80,xloc+88,120+f5m[1,j-1]*mult-80,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  f5p <- f5[[2*i]]$pat
  for(j in 2:ncol(f5p)) {
    if(f5p[2,j]==2)
      rect(xloc+102,120+f5p[1,j]*mult-80,xloc+114,120+f5p[1,j-1]*mult-80,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  rect(xloc-6,160-80,xloc+6,120-80,    lend=1, ljoin=1)
  rect(xloc+20,160-80,xloc+32,120-80,  lend=1, ljoin=1)
  rect(xloc+76,160-80,xloc+88,120-80,  lend=1, ljoin=1)
  rect(xloc+102,160-80,xloc+114,120-80,lend=1, ljoin=1)

  if(i==3) {
    arrows(xloc+54,87-90,xloc+54,80-90,len=0.04)
    arrows(xloc+54,130-100,xloc+54,80-90,len=0.04,lty=3)
  }  
  
  xloc <- xloc+78+120+30
}

text(u[1]/1.2,140-80,expression(F[5]),adj=c(0,0.5),cex=1.0)

xloc <- 38
for(i in 1:5) {
  rect(xloc-6,160-180,xloc+6,120-180,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+20,160-180,xloc+32,120-180,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+76,160-180,xloc+88,120-180,col=color[1],border=color[1], lend=1, ljoin=1)
  rect(xloc+102,160-180,xloc+114,120-180,col=color[1],border=color[1], lend=1, ljoin=1)

  fkm <- fk[[2*i-1]]$mat
  for(j in 2:ncol(fkm)) {
    if(fkm[2,j]==2)
      rect(xloc-6,120+fkm[1,j]*mult-180,xloc+6,120+fkm[1,j-1]*mult-180,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  fkp <- fk[[2*i-1]]$pat
  for(j in 2:ncol(fkp)) {
    if(fkp[2,j]==2)
      rect(xloc+20,120+fkp[1,j]*mult-180,xloc+32,120+fkp[1,j-1]*mult-180,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  fkm <- fk[[2*i]]$mat
  for(j in 2:ncol(fkm)) {
    if(fkm[2,j]==2)
      rect(xloc+76,120+fkm[1,j]*mult-180,xloc+88,120+fkm[1,j-1]*mult-180,col=color[5],border=color[5], lend=1, ljoin=1)
  }
  fkp <- fk[[2*i]]$pat
  for(j in 2:ncol(fkp)) {
    if(fkp[2,j]==2)
      rect(xloc+102,120+fkp[1,j]*mult-180,xloc+114,120+fkp[1,j-1]*mult-180,col=color[5],border=color[5], lend=1, ljoin=1)
  }

  rect(xloc-6,160-180,xloc+6,120-180,    lend=1, ljoin=1)
  rect(xloc+20,160-180,xloc+32,120-180,  lend=1, ljoin=1)
  rect(xloc+76,160-180,xloc+88,120-180,  lend=1, ljoin=1)
  rect(xloc+102,160-180,xloc+114,120-180,lend=1, ljoin=1)

  xloc <- xloc+78+120+30
}

text(u[1]/1.2,140-180,expression(F[k]),adj=c(0,0.5),cex=1.0)


dev.off()
