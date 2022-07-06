################
#
# Gridded weight calculator

####
# library and function load

# load libraries

library(ncdf4)
library(fields)
library(maps)

# load weight calc function

weight_calc = function(xvals,filelist=projfiles,subs = GCMs){
  
  xmdf = xvals[-length(xvals)]
  xobs = xvals[length(xvals)]
  w = c()
  
  for(i in 1:length(subs)){
    
    idx = grep(subs[i],filelist)    
    temp = xmdf[idx]
    
    w[i] = abs(1/(round(xobs,5)+abs(round(mean(temp,na.rm=TRUE),5)-round(xobs,5)))) # for inequal weights
    #w[i] = 1/length(subs) # for equal weights	
    
  }
  
  # Add in correction of the occasional Inf. value in w
  #if(length(which(is.na(w)==TRUE))>0){
  #w[which(is.na(w)==TRUE)]=min(w[-which(is.na(w)==TRUE)],na.rm=TRUE)
  #}
  
  if(length(which(w==Inf))>0){
    w[which(w==Inf)]=min(w[-which(w==Inf)],na.rm=TRUE)
  }
  
  if(all(is.na(xvals)==FALSE)==TRUE & all(w==Inf)==TRUE){
    message("Applying equal weight in situation where all values match obs")
    w = rep(1/length(subs),length(subs))
  }
  
  W = w/sum(w,na.rm=TRUE)
  
  if(any(is.na(w)==TRUE)) W=rep(NA,length(w))
  
  W  
}

###
# set parameters

varname = "TAVG_12M"
weighttype = "GCM"
obsdatasetname = "08144500"
regfilepath = "regfits"
anfilepath = "analysis"

refyear = 2005 # reference year for weighting calculation
#refperiod = 1981:2000 #reference climo period for testing with the weights
###
# Gather Projection Data - fitted anomalies for reference year only

# determine files
load(paste("/home/woot0002/uncertainty/",regfilepath,"/CPREP_",varname,"_",obsdatasetname,"_regfit.Rdata",sep=""))
load(paste("/home/woot0002/uncertainty/",regfilepath,"/",varname,"_",obsdatasetname,"_regfit.Rdata",sep=""))

memnames = names(fittedvals)[2:ncol(fittedvals)]

years = 1981:2099
yearidx = which(years==refyear) # reference period
#yearidx = which(years %in% refperiod) # climo period idx

# Gather projection data for reference year

projrefyeararray = c()

for(f in 1:length(memnames)){

projrefyeararray[f] = fittedvals[yearidx,(f+1)] # single year test

}

###
# Gather observed data - fitted anomaly for reference year only

obsrefyear = obstable[which(obstable[,1]==refyear),2]

###
# Run Weighting calculation function for GCM Weights

splitnames = do.call(rbind,strsplit(memnames,"_",fixed=TRUE))

GCMs = unique(splitnames[,1])
DSs = unique(paste(splitnames[,2],splitnames[,3],sep="_"))

#xmdf = projrefyeararray[18,12,]
#xobs = obsdatain[18,12]

xvals = c(projrefyeararray,obsrefyear)

GCMweights = weight_calc(xvals,filelist=memnames,subs=GCMs)
DSweights = weight_calc(xvals,filelist=memnames,subs=DSs)


##
# write weights to netcdf file

save(list=c("GCMweights","DSweights","GCMs","DSs"),file=paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",obsdatasetname,"_weights.Rdata",sep=""))

