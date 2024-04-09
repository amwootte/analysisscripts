# modified and updated by A.Wootten, University of Oklahoma, October 2022.
# previously, Duncan created this requiring that LH and CH have the same length,
# which presents problems with the different model calendars.
# This script has been adjusted to allow for differing length of LH and CH
# change implemented 10/4/22 by AMW.
#
# ERQM v3 incorporates the crucial components from the MBC package
# this version is the ERQM simulation with a cube-root transformation following the trace adjustment.

callDS<-function(target=NA,df.hist=NA,df.fut=NA){ 
  # 'Performs an equidistant correction adjustment
  # LH: Local Historical (a.k.a. observations)
  # CH: Coarse Historical (a.k.a. GCM historical)
  # CF: Coarse Future (a.k.a GCM future)
  # 'Cites Li et. al. 2010
  # ' Calls latest version of the ERQM function
  # ' as of 12-29
# changes to use MBC package code

#target=target,df.hist=df.hist,df.fut=df.fut

  lengthCF<-length(df.fut[[1]])
  lengthCH<-length(df.hist[[1]])
  #lengthLH<-length(target)
  
  CF.dim <- lengthCF
  CH.dim <- lengthCH
  
  
  # initialize data.frame
  temp<-data.frame(index=seq(1,CF.dim),CF=rep(NA,CF.dim))
  #df<-data.frame(index=seq(1,CH.dim),CH=rep(NA,CH.dim),LH=rep(NA,CH.dim)) # commented out per change described above.
  df<-data.frame(index=seq(1,CH.dim),CH=rep(NA,CH.dim))
  temp$CF<-df.fut[[1]]
  df$CH<-df.hist[[1]]
  #df$LH<-target # commented out per change described above.
  LH<-target

# adding in jitter factor in line with MBC package
jitter.factor=0 # seeing as the jitter factor is generally not used (always zero and ifs are not met) - 9/26/2023 AMW.
	#if(jitter.factor==0 && 
      #(length(unique(LH))==1 ||
       #length(unique(df$CH))==1 ||
       #length(unique(temp$CF))==1)){
        #jitter.factor <- sqrt(.Machine$double.eps)
    #}
    if(jitter.factor > 0){
        LH <- jitter(LH, jitter.factor)
        df$CH <- jitter(df$CH, jitter.factor)
        temp$CF <- jitter(temp$CF, jitter.factor)
    }
  
  if(any(target>1,na.rm=TRUE) & any(df.hist>1,na.rm=TRUE) & any(df.fut>1)){
  	trace=0.05
  } else {
	  trace=0.05/86400
  }
  trace.calc=0.5*trace                                 
  ratio.max=2 
  ratio.max.trace=10*trace

  epsilon <- .Machine$double.eps
  # Adjustment to account for precipitation values less than trace in simulations
        LH[which(LH < trace.calc)] <- runif(sum(LH < trace.calc,na.rm=TRUE), min=epsilon,max=trace.calc)
        df$CH[which(df$CH < trace.calc)] <- runif(sum(df$CH < trace.calc,na.rm=TRUE), min=epsilon, max=trace.calc)                            
        temp$CF[which(temp$CF < trace.calc)] <- runif(sum(temp$CF < trace.calc,na.rm=TRUE), min=epsilon,max=trace.calc)
 
  # drop NA's due to windowing and/or kfold masking
  temp <- na.omit(temp)

# Cubic root transformation
LH_trans = LH^(1/3)
  CH_trans = df$CH^(1/3)
  CF_trans = temp$CF^(1/3)
ratio.max.trace_trans = ratio.max.trace^(1/3)
  
  n.tau <- length(temp$CF)
  tau <- seq(0,1,length=n.tau)
  pp.type <- 7

## ERQM calculations
if(!all(is.na(LH)) & !all(is.na(CH_trans)) & !all(is.na(CF_trans)) & length(LH[which(is.na(LH)==FALSE)])>1 & length(df$CH[which(is.na(df$CH)==FALSE)])>1 & length(temp$CF[which(is.na(temp$CF)==FALSE)])>1 ){ # proceed with final calc if the LH has nonzero values post masking. Created 9/26/2023 AMW
# moved the start of this if up to address other problems with approx function when there's limited data available. 9/27/2023
# added to the if to include the GCM, to make sure it has nonzero values, also making sure they have more than 1 point to work with 9/28/2023.
  quant.LH <- quantile(LH_trans,tau,type=pp.type,na.rm=TRUE)
  quant.CH <- quantile(CH_trans,tau,type=pp.type,na.rm=TRUE)
  quant.CF <- quantile(CF_trans,tau,type=pp.type,na.rm=TRUE)
  tau.CF <- approx(quant.CF,tau,CF_trans,rule=2)$y

approx.t.qmc.tmp <- approx(tau, quant.CH, tau.CF, rule=2)$y
        delta.m <- CF_trans/approx.t.qmc.tmp
        delta.m[(delta.m > ratio.max) &
                (approx.t.qmc.tmp < ratio.max.trace_trans)] <- ratio.max

        ERQM_trans <-approx(tau, quant.LH, tau.CF, rule=2)$y*delta.m  
  # EQUIRatio CDF (Cannon et al 2015, MBC R package)
} else { #correction to account for any issue with masking that makes target to all zero. Created 9/26/2023 AMW
ERQM_trans <- rep(0,n.tau)
}

temp$EquiDistant<- ERQM_trans^3 # transform back out of cube root.
temp$EquiDistant[which(temp$EquiDistant<trace)] <-0 # values less than trace are set to zero

  # merge back into full data frame to preserve time order
  full.tmp <- data.frame(index=seq(1,CF.dim),CF=rep(NA,CF.dim))
  full <- merge(full.tmp, temp, by.x="index", by.y="index", all=T)
  full$CF <- full$CF.y
  full$CF.x <- NULL
  full$CF.y <- NULL
  ds.predict <- full$EquiDistant #full$EquiDistant is the ERQM result, temp$EquiDistant from above. 
  
  # add dimensions and name these
  # this step is necessary to 'abind' these output into the fuller DS output array
  dim(ds.predict) <- c(1,1,length(ds.predict),1,1)
  t.name <- seq(1,dims.fut[3])
  dimnames(ds.predict) <- list(i.index, j.index, t.name, w.name[[window]], k.name[[kfold]])
  # These are ERQM fit parameters
  # grabbing the summary() of CH, LH, CF and downscaled CF
  # first 4 moments would be good also, but simple use require the "moments" package
  # put these into a data frame 'df.fit' and index this with [i,j,window,kfold]
  # Add other parameters and fit statistics here as needed
  acf <- quantile(full$CF,na.rm=T)
  aeq <- quantile(full$EquiDistant,na.rm=T)
  alh <- quantile(LH,na.rm=T)
  ach <- quantile(df$CH,na.rm=T)
  # All statistics on a single row right now, change format as needed

  df.fit <- data.frame(i.index,j.index,w.name[[window]],k.name[[kfold]],
                       acf[1],acf[2],acf[3],acf[4],acf[5],median(full$CF,na.rm=T),
                       aeq[1],aeq[2],aeq[3],aeq[4],aeq[5],median(full$EquiDistant,na.rm=T),
                       alh[1],alh[2],alh[3],alh[4],alh[5],median(LH,na.rm=T),
                       ach[1],ach[2],ach[3],ach[4],ach[5],median(df$CH,na.rm=T))
  colnames(df.fit)[5:10]<-c("CF.0%","CF.25%","CF.50%","CF.75%","CF.100%","CF.median")
  colnames(df.fit)[11:16]<-c("CFEQ.0%","CFEQ.25%","CFEQ.50%","CFEQ.75%","CFEQ.100%","CF.median")
  colnames(df.fit)[17:22]<-c("LH.0%","LH.25%","LH.50%","LH.75%","LH.100%","LH.median")
  colnames(df.fit)[23:28]<-c("CH.0%","CH.25%","CH.50%","CH.75%","CH.100%","CH.median")
  rownames(df.fit)<-NULL
  return(list("ds.predict"=ds.predict, "df.fit"=df.fit))

}
