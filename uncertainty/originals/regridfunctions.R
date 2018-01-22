
interp.surface.gridfix = function (obj, grid.list) 
{
  x <- grid.list$x
  y <- grid.list$y
  M <- length(x)
  N <- length(y)
  out <- matrix(NA, nrow = M, ncol = N)
  for (i in 1:M) {
  
    if(i!=21){
      out[i, ] <- interp.surface2(obj, cbind(rep(x[i], N), y))
    }  else {
      out[i, ] <- interp.surfacefix(obj, cbind(rep(x[i], N), y))
    }
    
  }
  list(x = x, y = y, z = out)
}

interp.surface2 = function (obj, loc) {
  x <- obj$x
  y <- obj$y
  z <- obj$z
  nx <- length(x)
  ny <- length(y)
  lx <- approx(x, 1:nx, loc[, 1])$y
  ly <- approx(y, 1:ny, loc[, 2])$y
  lx1 <- floor(lx)
  ly1 <- floor(ly)
  ex <- lx - lx1
  ey <- ly - ly1
  ex[lx1 == nx] <- 1
  ey[ly1 == ny] <- 1
  lx1[lx1 == nx] <- nx - 1
  ly1[ly1 == ny] <- ny - 1
  
  part1 = z[cbind(lx1, ly1)]
  part2 = z[cbind(lx1 + 1, ly1)]
  part3 = z[cbind(lx1, ly1 + 1)]
  part4 = z[cbind(lx1 + 1, ly1 + 1)]
  
  partmat = matrix(NA,nrow=length(part1),ncol=4)
  partmat[,1]=part1
  partmat[,2]=part2
  partmat[,3]=part3
  partmat[,4]=part4
  
  countmat = ifelse(is.na(partmat)==TRUE,0,1)
  countersum = apply(countmat,1,sum)
  
  out1 = rep(NA,nrow(countmat))
  
  for(j in 1:length(countersum)){
    if(countersum[j]>1 & countersum[j]<4) out1[j] = mean(partmat[j,],na.rm=TRUE) 
    if(countersum[j]==4) out1[j] = partmat[j,1] * (1 - ex[j]) * (1 - ey[j])+partmat[j,2] * ex[j] * (1 - ey[j])+partmat[j,3] * (1 - ex[j]) * ey[j]+partmat[j,4] * ex[j] * ey[j] 
    #print(out1[j])
  }
  
  return(out1)
  
}



interp.surfacefix = function (obj, loc) {
  x <- obj$x
  y <- obj$y
  z <- obj$z
  nx <- length(x)
  ny <- length(y)
  lx <- approx(x, 1:nx, loc[, 1])$y
  ly <- approx(y, 1:ny, loc[, 2])$y
  lx1 <- floor(lx)
  ly1 <- floor(ly)
  ex <- lx - lx1
  ey <- ly - ly1
  ex[lx1 == nx] <- 1
  ey[ly1 == ny] <- 1
  lx1[lx1 == nx] <- nx - 1
  ly1[ly1 == ny] <- ny - 1
  
  out1 = z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey)+
  z[cbind(lx1 + 1, ly1)] * ex * (1 - ey)+  
  z[cbind(lx1, ly1 + 1)] * (1 - ex) * ey+
  z[cbind(lx1 + 1, ly1 + 1)] * ex * ey
  
  EX = lx[10] - (lx1[10]-1)
  EY = ly[10] - (ly1[10]-1)
  
  out1[10] = z[cbind(lx1[10]-1, ly1[10]-1)] * (1 - ex[10]) * (1 - ey[10])+
    z[cbind(lx1[10] + 2, ly1[10]-1)] * ex[10] * (1 - ey[10])+  
    z[cbind(lx1[10]-1, ly1[10] + 2)] * (1 - ex[10]) * ey[10]+
    z[cbind(lx1[10] + 2, ly1[10] + 2)] * ex[10] * ey[10]
  
  return(out1)
                                                                                                                  
}
