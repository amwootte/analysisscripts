# modified and updated by A.Wootten, University of Oklahoma, October 2022.
# previously, Duncan created this requiring that LH and CH have the same length,
# which presents problems with the different model calendars.
# This script has been adjusted to allow for differing length of LH and CH
# change implemented 10/4/22 by AMW.
#
# EDQM v3 incorporates the crucial components from the MBC package

callDS<-function(target=NA,df.hist=NA,df.fut=NA){ 
  #'Performs an equidistant correction adjustment
  # LH: Local Historical (a.k.a. observations)
  # CH: Coarse Historical (a.k.a. GCM historical)
  # CF: Coarse Future (a.k.a GCM future)
  #'Cites Li et. al. 2010
  #' Calls latest version of the EDQM function
  #' as of 12-29
# using MBC package code

  lengthCF<-length(df.fut[[1]])
  lengthCH<-length(df.hist[[1]])
  #lengthLH<-length(target)
  
  CF.dim <- lengthCF
  CH.dim <- lengthCH
  
  
  # initialize data.frame
  temp<-data.frame(index=seq(1,CF.dim),CF=rep(NA,CF.dim))
  df<-data.frame(index=seq(1,CH.dim),CH=rep(NA,CH.dim))
  temp$CF<-df.fut[[1]]
  df$CH<-df.hist[[1]]

# jitter factor (MBC package)
jitter.factor=0
	if(jitter.factor==0 && 
      (length(unique(target))==1 ||
       length(unique(df$CH))==1 ||
       length(unique(temp$CF))==1)){
        jitter.factor <- sqrt(.Machine$double.eps)
    }
    if(jitter.factor > 0){
        target <- jitter(target, jitter.factor)
        df$CH <- jitter(df$CH, jitter.factor)
        temp$CF <- jitter(temp$CF, jitter.factor)
    }

  #temp<-temp[order(temp$CF),]
  # drop NA's due to windowing and/or kfold masking
  temp <- na.omit(temp)
  n.tau <- length(temp$CF)
  tau <- seq(0,1,length=n.tau)
  pp.type <- 7
  quant.target <- quantile(target,tau,type=pp.type,na.rm=TRUE)
  quant.CH <- quantile(df$CH,tau,type=pp.type,na.rm=TRUE)
  quant.CF <- quantile(temp$CF,tau,type=pp.type,na.rm=TRUE)
  tau.CF <- approx(quant.CF,tau,temp$CF,rule=2)$y
  delta.m <- temp$CF - approx(tau,quant.CH,tau.CF,rule=2)$y
  # EQUIDISTANT CDF (Cannon et al 2015, MBC R package)

  temp$EquiDistant<- approx(tau,quant.target,tau.CF,rule=2)$y + delta.m

  # merge back into full data frame to preserve time order
  #temp<-temp[order(temp$index),]
  full.tmp <- data.frame(index=seq(1,CF.dim),CF=rep(NA,CF.dim))
  full <- merge(full.tmp, temp, by.x="index", by.y="index", all=T)
  full$CF <- full$CF.y
  full$CF.x <- NULL
  full$CF.y <- NULL
  ds.predict <- full$EquiDistant # full$EquiDistant is the downscaled result.

  # add dimensions and name these
  # this step is necessary to 'abind' these output into the fuller DS output array
  dim(ds.predict) <- c(1,1,length(ds.predict),1,1)
  t.name <- seq(1,dims.fut[3])
  dimnames(ds.predict) <- list(i.index, j.index, t.name, w.name[[window]], k.name[[kfold]])
  # These are EDQM fit parameters
  # grabbing the summary() of CH, LH, CF and downscaled CF
  # first 4 moments would be good also, but simple use require the "moments" package
  # put these into a data frame 'df.fit' and index this with [i,j,window,kfold]
  # Add other parameters and fit statistics here as needed
  acf <- quantile(full$CF,na.rm=T)
  aeq <- quantile(full$EquiDistant,na.rm=T)
  alh <- quantile(target,na.rm=T)
  ach <- quantile(df$CH,na.rm=T)
  # All statistics on a single row right now, change format as needed

  df.fit <- data.frame(i.index,j.index,w.name[[window]],k.name[[kfold]],
                       acf[1],acf[2],acf[3],acf[4],acf[5],median(full$CF,na.rm=T),
                       aeq[1],aeq[2],aeq[3],aeq[4],aeq[5],median(full$EquiDistant,na.rm=T),
                       alh[1],alh[2],alh[3],alh[4],alh[5],median(target,na.rm=T),
                       ach[1],ach[2],ach[3],ach[4],ach[5],median(df$CH,na.rm=T))
  colnames(df.fit)[5:10]<-c("CF.0%","CF.25%","CF.50%","CF.75%","CF.100%","CF.median")
  colnames(df.fit)[11:16]<-c("CFEQ.0%","CFEQ.25%","CFEQ.50%","CFEQ.75%","CFEQ.100%","CF.median")
  colnames(df.fit)[17:22]<-c("LH.0%","LH.25%","LH.50%","LH.75%","LH.100%","LH.median")
  colnames(df.fit)[23:28]<-c("CH.0%","CH.25%","CH.50%","CH.75%","CH.100%","CH.median")
  rownames(df.fit)<-NULL
rm(full)
rm(df)
rm(temp)
gc()
gc()
  return(list("ds.predict"=ds.predict, "df.fit"=df.fit))

}
