##########################
#
# Bias Correction Functions

#####################
# Delta method bias correction

deltaBC = function(modeldata,modclimo,obsdailyclimo,timesmod){

# modeldata - the model data to be corrected, expecting an array dimensions lon,lat,time
# modclimo - monthly climatology, expecting an array dimensions lon,lat,month
# obsdailyclimo - daily observed climatology, expecting an array dimensions lon,lat,dayofyear
# ABOVE THREE MUST HAVE MATCHING DIMENSIONS IN LON AND LAT
# timesmod - times in the model, expecting a vector character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modeldata

if(dim(modclimo)[c(1,2)] != dim(obsdailyclimo)[c(1,2)]){
stop("Model and Observed data don't have matching space dimensions",call.=TRUE)
}

if(dim(modeldata)[3] != length(timesmod)){
stop("Time vector doesn't have the same length as model data",call.=TRUE)
}

outputbc = array(NA,dim=dim(modeldata))
for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7)) # month idx for the monthly obs climatology
dayidx = which(days==substr(timesmod[i],6,10)) # day idx for the daily obs climatology
outputbc[,,i] = obsdayclimo[,,dayidx]-modclimo[,,monidx]+modeldata[,,i] # BC calc
}

outputbc

}

#####################
# Unbiasing method bias correction

unbiasingBC = function(modeldata,modclimo,obsclimo,timesmod){

# modeldata - the model data to be corrected, expecting an array dimensions lon,lat,time
# modclimo - monthly climatology, expecting an array dimensions lon,lat,month
# obsdailyclimo - monthly observed climatology, expecting an array dimensions lon,lat,month
# ABOVE THREE MUST HAVE MATCHING DIMENSIONS IN LON AND LAT
# timesmod - times in the model, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modeldata

if(dim(modclimo)[c(1,2)] != dim(obsclimo)[c(1,2)]){
stop("Model and Observed data don't have matching space dimensions",call.=TRUE)
}

if(dim(modeldata)[3] != length(timesmod)){
stop("Time vector doesn't have the same length as model data",call.=TRUE)
}

outputbc = array(NA,dim=dim(modeldata))

for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7)) # gets the monthly idx 
outputbc[,,i] = modeldata[,,i]-modclimo[,,monidx]+obsclimo[,,monidx] 
}

outputbc

}

##########################
# Scaling method bias correction

scalingBC = function(modeldata,modclimo,obsclimo,timesmod){

# modeldata - the model data to be corrected, expecting an array dimensions lon,lat,time
# modclimo - monthly climatology, expecting an array dimensions lon,lat,month
# obsdailyclimo - monthly observed climatology, expecting an array dimensions lon,lat,month
# ABOVE THREE MUST HAVE MATCHING DIMENSIONS IN LON AND LAT
# timesmod - times in the model, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modeldata

#modeldata = outputdata[[3]]
#modclimo = outputclimo[[3]]
#obsclimo = obsclimo
#timesmodold = timesmod
#timesmod = timesmod[[3]]

if(dim(modclimo)[c(1,2)] != dim(obsclimo)[c(1,2)]){
stop("Model and Observed data don't have matching space dimensions",call.=TRUE)
}

if(dim(modeldata)[3] != length(timesmod)){
stop("Time vector doesn't have the same length as model data",call.=TRUE)
}

outputbc = array(NA,dim=dim(modeldata))

for(i in 1:dim(outputbc)[3]){
monidx = as.numeric(substr(timesmod[i],6,7)) # gets the monthly idx 
outputbc[,,i] = (modeldata[,,i]/modclimo[,,monidx])*obsclimo[,,monidx] # scaling calc
}


#modeldatatest = modeldata[,,i]+273.15
#modclimotest = modclimo[,,as.numeric(substr(timesmod[i],6,7))]+273.15
#obsclimotest = obsclimo[,,as.numeric(substr(timesmod[i],6,7))]+273.15
#outputcheck = (((modeldatatest/modclimotest)*obsclimotest)-273.15

#testsfc1 = list(x=lon,y=lat,z=modeldata[,,i])
#testsfc2 = list(x=lon,y=lat,z=modclimo[,,as.numeric(substr(timesmod[i],6,7))])
#testsfc3 = list(x=lon,y=lat,z=outputbc[,,i])
#testsfc4 = list(x=lon,y=lat,z=outputcheck)

#par(mfrow=c(1,4))
#surface(testsfc1,zlim=c(-2,16))
#map("state",add=TRUE)
#surface(testsfc2,zlim=c(-2,16))
#map("state",add=TRUE)
#surface(testsfc3)
#map("state",add=TRUE)
#surface(testsfc4)
#map("state",add=TRUE)


outputbc

}

#############################
# QQmap method bias correction

####
# QQmap - unbounded variables only

qqmapBC = function(modeldata,modelhist,timesmodhist,timesmod,obsdata){

# modeldata - the model data to be corrected, expecting an array dimensions lon,lat,time
# modelhist - historical baseline period daily, expecting an array dimensions lon,lat,time
# obsdata - observed daily data, expecting an array dimensions lon,lat,time
# ABOVE THREE MUST HAVE MATCHING DIMENSIONS IN LON AND LAT
# ALSO modelhist AND obsdata SHOULD MATCH IN ALL THREE DIMENSIONS

# timesmod - times in the model data, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modeldata

# timesmodHITS - times in the model baseline period, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modelhist AND obsdata

if(dim(modelhist) != dim(obsdata)){
stop("Model Baseline and Observed data don't have matching dimensions",call.=TRUE)
}

if(dim(modeldata)[c(1,2)] != dim(obsdata)[c(1,2)]){
stop("Model data and Observed data don't have matching space dimensions",call.=TRUE)
}

if(dim(modeldata)[3] != length(timesmod)){
stop("Time vector doesn't have the same length as model data",call.=TRUE)
}

if(dim(modelhist)[3] != length(timesmodhist)){
stop("historical Time vector doesn't have the same length as model historical data",call.=TRUE)
}

outputbc = array(NA,dim=dim(modeldata))

for(m in 1:12){

monidx = which(as.numeric(substr(timesmod,6,7))==m)
monidxhist = which(as.numeric(substr(timesmodhist,6,7))==m)

for(r in 1:dim(outputbc)[1]){
for(c in 1:dim(outputbc)[2]){
if(any(!is.na(modelhist[r,c,monidxhist])) & any(!is.na(obsdata[r,c,monidxhist]))){
ePrd = ecdf(modelhist[r,c,monidxhist])
sim = quantile(obsdata[r,c,monidxhist],probs=ePrd(modeldata[r,c,monidx]),na.rm=TRUE,type=4)
outputbc[r,c,monidx] = sim
}
message("Finished calc for row ",r," and col ",c," and month ",m)
}}}

outputbc

}


####
# QQmap - bounded variables - still testing this.

qqmapBC_pr = function(modeldata,modelhist,timesmodhist,timesmod,obsdata,threshold=2.5){

# modeldata - the model data to be corrected, expecting an array dimensions lon,lat,time
# modelhist - historical baseline period daily, expecting an array dimensions lon,lat,time
# obsdata - observed daily data, expecting an array dimensions lon,lat,time
# ABOVE THREE MUST HAVE MATCHING DIMENSIONS IN LON AND LAT
# ALSO modelhist AND obsdata SHOULD MATCH IN ALL THREE DIMENSIONS

# timesmod - times in the model data, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modeldata

# timesmodHITS - times in the model baseline period, expecting a character string format "YYYY-MM-DD"
# ABOVE MUST HAVE THE SAME LENGTH AS TIME DIMENSION OF modelhist AND obsdata

modeldata = outputdata[[1]]
modelhist=outputdata[[1]][,,which(as.numeric(substr(timesmod[[1]],1,4))>=1981 & as.numeric(substr(timesmod[[1]],1,4))<=1999)]
timesmodhist = timesmod[[1]][which(as.numeric(substr(timesmod[[1]],1,4))>=1981 & as.numeric(substr(timesmod[[1]],1,4))<=1999)]
timesmod = timesmod[[1]]
obsdata=vardata
threshold=2.5


if(dim(modelhist) != dim(obsdata)){
stop("Model Baseline and Observed data don't have matching dimensions",call.=TRUE)
}

if(dim(modeldata)[c(1,2)] != dim(obsdata)[c(1,2)]){
stop("Model data and Observed data don't have matching space dimensions",call.=TRUE)
}

if(dim(modeldata)[3] != length(timesmod)){
stop("Time vector doesn't have the same length as model data",call.=TRUE)
}

if(dim(modelhist)[3] != length(timesmodhist)){
stop("historical Time vector doesn't have the same length as model historical data",call.=TRUE)
}

outputbc = array(NA,dim=dim(modeldata))

for(m in 1:12){
for(r in 1:dim(outputbc)[1]){
for(c in 1:dim(outputbc)[2]){

if(any(!is.na(modelhist[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)])) & any(!is.na(obsdata[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)]))){

nP<-sum(as.double(obsdata[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)]<=threshold & !is.na(obsdata[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)])), na.rm = TRUE)

        if (nP>0 & nP<dim(obsdata[,,which(as.numeric(substr(timesmodhist,6,7))==m)])[3]){
          ix<-sort(modelhist[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
          Ps<-sort(modelhist[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)], decreasing = FALSE, na.last = NA)
          Pth<-Ps[nP+1]
          if (Ps[nP+1]<=threshold){
            Os<-sort(obsdata[r,c,which(as.numeric(substr(timesmodhist,6,7))==m)], decreasing = FALSE, na.last = NA)
            ind<-which(Ps > threshold & !is.na(Ps))
            if (length(ind)==0){
              ind <- max(which(!is.na(Ps)))
              ind <- min(c(length(Os),ind))
            }else{
              ind<-min(which(Ps > threshold & !is.na(Ps)))
            }

            # [Shape parameter Scale parameter]
            if (length(unique(Os[(nP+1):ind]))<6){
		 Ps[(nP+1):ind] <- mean(Os[(nP+1):ind], na.rm = TRUE)
	    }else{
	      auxGamma<-fitdistr(Os[(nP+1):ind],"gamma")
              Ps[(nP+1):ind]<-rgamma(ind-nP, auxGamma$estimate[1], rate = auxGamma$estimate[2])
	    }
 	      Ps<-sort(Ps, decreasing = FALSE, na.last = NA)
          }
          ind<-min(nP,dim(modelhist[,,which(as.numeric(substr(timesmodhist,6,7))==m)])[3])
          Ps[1:ind]<-0
          modelhist[r,c,ix]<-Ps
	}else{
          if (nP==dim(obs)[3]){
            ix<-sort(modelhist[r,c,], decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
            Ps<-sort(modelhist[r,c,], decreasing = FALSE, na.last = NA)
            Pth<-Ps[nP]
            ind<-min(nP,dim(modelhist)[3])
            Ps[1:ind]<-0
            pred[r,c,ix]<-Ps
          }
        }

monidx = which(as.numeric(substr(timesmod,6,7))==m)
monidxhist = which(as.numeric(substr(timesmodhist,6,7))==m)

sim = c()
ePrd<-ecdf(modelhist[r,c,which(modelhist[r,c,]>Pth & as.numeric(substr(timesmodhist,6,7))==m)])

noRain<-which(modeldata[r,c,]<=Pth & !is.na(modeldata[r,c,]) & as.numeric(substr(timesmod,6,7))==m)
rain<-which(modeldata[r,c,]>Pth & !is.na(modeldata[r,c,]) & as.numeric(substr(timesmod,6,7))==m)
drizzle<-which(modeldata[r,c,]>Pth & modeldata[r,c,] <= min(modeldata[r,c,which(modeldata[r,c,]>Pth & as.numeric(substr(timesmod,6,7))==m)], na.rm = TRUE) & !is.na(modeldata[r,c,] & as.numeric(substr(timesmod,6,7))==m))

eFrc<-ecdf(modeldata[r,c,rain])

if (length(drizzle)>0){
   sim[drizzle]<-quantile(modeldata[r,c,which(modeldata[r,c,]>min(modelhist[r,c,which(modelhist[r,c,]>Pth & as.numeric(substr(timesmodhist,6,7))==m)], na.rm = TRUE) & !is.na(modeldata[r,c,]) & as.numeric(substr(timesmod,6,7))==m)], probs = eFrc(modeldata[r,c,drizzle]), na.rm = TRUE, type = 4)
if(length(rain)>0) sim[rain]<-quantile(obsdata[r,c,which(obsdata[r,c,] > threshold & !is.na(obsdata[r,c,]) & as.numeric(substr(timesmodhist,6,7))==m)], probs = eFrc(modeldata[r,c,rain]), na.rm = TRUE, type = 4)
} else {

if(length(rain)>0) sim[rain]<-quantile(obsdata[r,c,which(obsdata[r,c,]>threshold & !is.na(obsdata[r,c,] & as.numeric(substr(timesmodhist,6,7))==m))], probs = ePrd(modeldata[r,c,rain]), na.rm = TRUE, type = 4)


}

sim[noRain]<-0

outputbc[r,c,monidx] = sim




}
message("Finished calc for row ",r," and col ",c," and month ",m)
}}}

outputbc

}




















