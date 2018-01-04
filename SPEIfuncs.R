
require("SPEI")

ETcalc  = function(Tmax,Tmin,Pre,lat,na.rm=TRUE){
  
  # Tmax, Tmin, and Pre should be arrays in X,Y,T with the appropriate inputs to the
  # functions in the SPEI package for the hargreaves function
  # lat should be a vector with all the latitudes equal to the numbmer of columns
  
  ETout = array(NA,dim=dim(Tmax))
  for(r in 1:dim(Tmax)[1]){
    for(c in 1:dim(Tmax)[2]){
      if(all(is.na(Tmax[r,c,])==TRUE)==TRUE | all(is.na(Tmin[r,c,])==TRUE)==TRUE | all(is.na(Pre[r,c,])==TRUE)==TRUE){
        ETout[r,c,]=NA
      } else {
        ETout[r,c,]=hargreaves(Tmin[r,c,],Tmax[r,c,],Pre[r,c,],lat=lat[c],na.rm=na.rm)
      }
    }
  }
  ETout
}

##########

SPIcalc  = function(Pre,scale=1,na.rm=TRUE){
  
  # Tmax, Tmin, and Pre should be arrays in X,Y,T with the appropriate inputs to the
  # functions in the SPEI package for the hargreaves function
  # lat should be a vector with all the latitudes equal to the numbmer of columns
  
  SPIout = array(NA,dim=dim(Pre))
  for(r in 1:dim(Pre)[1]){
    for(c in 1:dim(Pre)[2]){
      if(all(is.na(Pre[r,c,])==TRUE)==TRUE){
        SPIout[r,c,]=NA
      } else {
        SPIout[r,c,]=spi(Pre[r,c,],scale=scale,na.rm=na.rm)
      }
    }
  }
  SPIout
}
##########
SPEIcalc  = function(Diff,scale=1,na.rm=TRUE){
  
  # Tmax, Tmin, and Pre should be arrays in X,Y,T with the appropriate inputs to the
  # functions in the SPEI package for the hargreaves function
  # lat should be a vector with all the latitudes equal to the numbmer of columns
  
  SPEIout = array(NA,dim=dim(Diff))
  for(r in 1:dim(Diff)[1]){
    for(c in 1:dim(Diff)[2]){
      if(all(is.na(Diff[r,c,])==TRUE)==TRUE){
        SPEIout[r,c,]=NA
      } else {
        SPEIout[r,c,]=spei(Diff[r,c,],scale=scale,na.rm=na.rm)
      }
    }
  }
  SPEIout
}
##########
