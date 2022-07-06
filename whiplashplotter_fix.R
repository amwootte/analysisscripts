
library(iterators)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(zoo)

whiplashcalc = function(PrecipData,lowperc=0.25,highperc=0.75,spread=1){
  
  # PrecipData: array of precipitation data arranged as [longitude,latitude,months]
  # lowperc: low percentile, defaults to 0.25 to 25th percentile
  # highperc: high percentile, defaults to 0.75 to 75th percentile
  
  #PrecipData = precipmon
  #lowperc = 0.25
  #highperc = 0.75
  
  percentile_bottom = apply(PrecipData,c(1,2),quantile,probs=lowperc,na.rm=TRUE)
  percentile_top = apply(PrecipData,c(1,2),quantile,probs=highperc,na.rm=TRUE)
  
  #testsfc = list(x=lon,y=lat,z=percentile_top)
  #surface(testsfc,type="I")
  
  precip_hardcode = PrecipData
  
  for(r in 1:dim(PrecipData)[1]){
    for(c in 1:dim(PrecipData)[2]){
      if(all(is.na(PrecipData[r,c,])==FALSE)==TRUE){
        precip_hardcode[r,c,]=ifelse(PrecipData[r,c,]<=percentile_bottom[r,c] | PrecipData[r,c,]>=percentile_top[r,c],ifelse(PrecipData[r,c,]<=percentile_bottom[r,c],0,2),1)
      }
    }
  }
  
  PDcount = DPcount = matrix(NA,nrow=dim(PrecipData)[1],ncol=dim(PrecipData)[2])
  precip_PD = precip_DP = PrecipData
  
  for(r in 1:dim(PrecipData)[1]){
    for(c in 1:dim(PrecipData)[2]){
      if(all(is.na(PrecipData[r,c,])==FALSE)==TRUE){
        tmp = precip_hardcode[r,c,1:(dim(PrecipData)[3]-spread)]-precip_hardcode[r,c,(spread+1):dim(PrecipData)[3]]
        tmpdp = c(ifelse(tmp==-2,1,0),rep(NA,spread))
        tmppd = c(ifelse(tmp==2,1,0),rep(NA,spread))
        precip_PD[r,c,] = tmppd
        precip_DP[r,c,] = tmpdp
        PDcount[r,c] = sum(tmppd,na.rm=TRUE)
        DPcount[r,c] = sum(tmpdp,na.rm=TRUE)
      }
    }
  }
  
  #testsfc = list(x=lon,y=lat,z=PDcount)
  #surface(testsfc,type="I")
  
  #testsfc = list(x=lon,y=lat,z=DPcount)
  #surface(testsfc,type="I")
  results = list(percentile_bottom=percentile_bottom,percentile_top=percentile_top,precip_hardcode=precip_hardcode,precip_PD=precip_PD,precip_DP=precip_DP,PDcount=PDcount,DPcount=DPcount)
  results
}
calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollapply(inputdata,width=size,FUN=sum,na.rm=TRUE)
}

source("/data2/3to5/I35/scripts/analysisfunctions.R")

##########
# Set input information

message("Setting inputs")

varname = "pr" # these are all required arguments for step 1
difftype = "absolute"
futureperiod = c(2070,2099)

# set basic arguments
varplot = "PDc30sp10" # these are all required arguments for step 1
seasonin = "ann"
titlepart ="Pluvial-Drought"

step1_filename=paste("/data2/3to5/I35/all_mems/",varplot,"_allmem_",difftype,"_",futureperiod[1],"-",futureperiod[2],"_",seasonin,".nc",sep="")

# set map plot arguments
BINLIMIT=30
histcolorchoice = "whitetogreen"
diffcolorchoice = "browntoblue"
diffbartype = "difference"

plotmaps=TRUE

shapefile = "/home/woot0002/shapefiles/WBD_HU8_SanSabaLlano"

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

############

varin = varname
if(varname=="tmax85" | varname=="tmax90" | varname=="tmax95" | varname=="tmax100" | varname=="gsl") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
if(varname=="prcptot" | varname=="prmean" | varname=="prmed" | varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="r1mm" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd" | varname=="R99" | varname=="R95" | varname=="R90"  | varname=="R99freq" | varname=="R95freq" | varname=="R90freq") varin="pr"

########
# Find all file names

###########
# 1. Data Gather and conversion
histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

projnotes = paste(projnotes,collapse=",")
histnotes = paste(histnotes,collapse=",")

histlist = paste(histfilelist,collapse=",")
projlist = paste(projfilelist,collapse=",")

###

test = nc_open(step1_filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
#diffs = ncvar_get(test,"projmeandiff")
projmean = ncvar_get(test,"projmean")/30
histmean = ncvar_get(test,"histmean")/25

varunits = test$var[[1]]$units
changeunits = test$var[[3]]$units

nc_close(test)

projfilebreakdown = do.call(rbind,strsplit(projnotes,",",fixed=TRUE))
projfilebreakdown = do.call(rbind,strsplit(projfilebreakdown,"_",fixed=TRUE))
projfilebreakdown = data.frame(projfilebreakdown)
names(projfilebreakdown) = c("GCM","DS","obs","scen")

histfilebreakdown = do.call(rbind,strsplit(histnotes,",",fixed=TRUE))
histfilebreakdown = do.call(rbind,strsplit(histfilebreakdown,"_",fixed=TRUE))
histfilebreakdown = data.frame(histfilebreakdown)
names(histfilebreakdown) = c("GCM","DS","obs")

###

diffs = projmean

for(i in 1:nrow(projfilebreakdown)){
  histidx = which(histfilebreakdown$GCM==projfilebreakdown$GCM[i] & histfilebreakdown$DS==projfilebreakdown$DS[i] & histfilebreakdown$obs==projfilebreakdown$obs[i])
  diffs[,,i]=projmean[,,i]-histmean[,,histidx]
}


###

if(plotmaps==TRUE){

diffcolorbar = colorramp(diffs,colorchoice=diffcolorchoice,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = FALSE)

diffs_sort = diffs[,,order(projfilebreakdown$scen)]
projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

plotfilename = paste("/data2/3to5/I35/plots/all_mems/IndividualMembers_",varplot,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"_",seasonin,"_change.pdf",sep="")

pdf(plotfilename,onefile=TRUE,width=10,height=10)
par(mfrow=c(3,3))
for(i in 1:nrow(projfilebreakdown)){
  GCM = projfilebreakdown$GCM[i]
  scen = projfilebreakdown$scen[i]
  obs = projfilebreakdown$obs[i]
  DS = projfilebreakdown$DS[i]
  
  testsfc1 = list(x=lon,y=lat,z=diffs_sort[,,i])
  surface(testsfc1,type="I",main=paste("Projected Difference from Historical Climate\nScen: ",scen," GCM: ",GCM," DS: ", DS," Obs: ",obs,sep=""),zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

####
projfilebreakdown$group = paste(projfilebreakdown$scen,projfilebreakdown$DS,sep="_")
groupunique = unique(projfilebreakdown$group)
scens = unique(projfilebreakdown$scen)

diffsgroupmean = array(NA,dim=c(length(lon),length(lat),length(groupunique)))
diffsscenmean = array(NA,dim=c(length(lon),length(lat),length(scens)))

for(g in 1:length(groupunique)){
  idx = which(projfilebreakdown$group==groupunique[g])
  diffsgroupmean[,,g]=apply(diffs_sort[,,idx],c(1,2),mean,na.rm=TRUE)
}

for(s in 1:length(scens)){
  idx = which(projfilebreakdown$scen==scens[s])
  diffsscenmean[,,s]=apply(diffs_sort[,,idx],c(1,2),mean,na.rm=TRUE)
}

diffcolorbar_group = colorramp(diffsgroupmean,colorchoice=diffcolorchoice,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = FALSE)
diffcolorbar_scen = colorramp(diffsscenmean,colorchoice=diffcolorchoice,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = FALSE)

plotfilename = paste("/data2/3to5/I35/plots/ens_means/EnsMean_",varplot,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"_",seasonin,"_change.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=15,height=5)
par(mfrow=c(1,3))
for(i in 1:length(scens)){
  testsfc1 = list(x=lon,y=lat,z=diffsscenmean[,,i])
  surface(testsfc1,type="I",main=paste("Mean Projected Difference ",scens[i],sep=""),zlim=diffcolorbar_scen[[1]],col=diffcolorbar_scen[[3]],breaks=diffcolorbar_scen[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

plotfilename = paste("/data2/3to5/I35/plots/ens_means/EnsMean_",varplot,"_",futureperiod[1],"-",futureperiod[2],"_",difftype,"_",seasonin,"_change_DS.pdf",sep="") 

pdf(plotfilename,onefile=TRUE,width=10,height=10)
par(mfrow=c(3,3))
for(i in 1:length(groupunique)){
  testsfc1 = list(x=lon,y=lat,z=diffsgroupmean[,,i])
  surface(testsfc1,type="I",main=paste("Mean Projected Difference ",groupunique[i],sep=""),zlim=diffcolorbar_group[[1]],col=diffcolorbar_group[[3]],breaks=diffcolorbar_group[[2]],xlab="Longitude",ylab="Latitude")
  map("state",add=TRUE)
}

dev.off()

}

###########
# boxplots for San Saba / Llano region

message("opening URG shapefile")
test = readShapePoly(shapefile)
#projection(test) <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
#test = spTransform(test, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
#projection(test) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
message("got shapefile loaded")

#diffs_sort = diffs[,,order(projfilebreakdown$scen)]
#projmean = projmean[,,order(projfilebreakdown$scen)]
#projfilebreakdown = projfilebreakdown[order(projfilebreakdown$scen),]

###
# get hist values

histvalmat = matrix(NA,nrow=nrow(histfilebreakdown),ncol=5)

for(h in 1:nrow(histfilebreakdown)){
  
    tmpV = histmean[,,h]
    tmpdat = histmean[,,h]
  
  modrasV = raster(t(tmpV)[length(lat):1,])
  #rm(tmpV)
  if(all(lon>0)){
    extent(modrasV) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasV) = c(min(lon),max(lon),min(lat),max(lat))
  }
  areas = eval(parse(text=paste("test@data$","HUC_8",sep="")))
  
  tmp=c()
  
  for(a in 1:length(eval(parse(text=paste("test@data$","HUC_8",sep=""))))){ # total is 5
    test.sub <- test[as.character(eval(parse(text=paste("test@data$","HUC_8",sep="")))) %in% areas[a], ] # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
    valueout=extract(modrasV,test.sub,weights=TRUE,fun=mean,na.rm=TRUE)  
    histvalmat[h,a] = valueout
  }
}

###
# get projected and difference values

projvalmat = diffvalmat = matrix(NA,nrow=nrow(projfilebreakdown),ncol=5)

for(h in 1:nrow(projfilebreakdown)){
  
  tmpP = projmean[,,h]
  tmpD = diffs_sort[,,h]
  
  modrasP = raster(t(tmpP)[length(lat):1,])
  modrasD = raster(t(tmpD)[length(lat):1,])
  #rm(tmpV)
  if(all(lon>0)){
    extent(modrasP) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
    extent(modrasD) = c(min(lon)-360,max(lon)-360,min(lat),max(lat))
  } else {
    extent(modrasP) = c(min(lon),max(lon),min(lat),max(lat))
    extent(modrasD) = c(min(lon),max(lon),min(lat),max(lat))
  }
  areas = eval(parse(text=paste("test@data$","HUC_8",sep="")))
  
  tmp=c()
  
  for(a in 1:length(eval(parse(text=paste("test@data$","HUC_8",sep=""))))){ # total is 5
    test.sub <- test[as.character(eval(parse(text=paste("test@data$","HUC_8",sep="")))) %in% areas[a], ] # in here are two important details. 1) column name in the data array, 2) item in the data array to subset by
    valueout=extract(modrasP,test.sub,weights=TRUE,fun=mean,na.rm=TRUE)  
    projvalmat[h,a] = valueout
    valueout=extract(modrasD,test.sub,weights=TRUE,fun=mean,na.rm=TRUE)  
    diffvalmat[h,a] = valueout
  }
}

#####
# plot values

histvalmean = apply(histvalmat,1,mean,na.rm=TRUE)
projvalmean = apply(projvalmat,1,mean,na.rm=TRUE)
diffvalmean = apply(diffvalmat,1,mean,na.rm=TRUE)

histfilebreakdown$histvalmean = histvalmean
projfilebreakdown$projvalmean = projvalmean
projfilebreakdown$diffvalmean = diffvalmean

library(ggplot2)

ggplot(projfilebreakdown, aes(x=scen, y=diffvalmean)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  ggtitle(paste("2070-2099 Projected Change ",titlepart,sep=""))+xlab("RCP")+ylab("Projected Change")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")

ggplot(projfilebreakdown, aes(x=scen, y=diffvalmean,fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  ggtitle(paste("2070-2099 Projected Change ",titlepart,sep=""))+xlab("RCP")+ylab("Projected Change")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed")

ggplot(projfilebreakdown, aes(x=DS, y=diffvalmean,fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  ggtitle(paste("2070-2099 Projected Change ",titlepart,sep=""))+xlab("Downscaling Technique")+ylab("Projected Change")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_grid(cols=vars(scen)) +geom_hline(yintercept=0,linetype="dashed")
  
aggregate(projfilebreakdown[,6:7],by=list(scen=projfilebreakdown$scen),mean,na.rm=TRUE)

aggregate(projfilebreakdown[,6:7],by=list(scenGCM=paste(projfilebreakdown$scen,projfilebreakdown$GCM,sep="_")),mean,na.rm=TRUE)

ggplot(histfilebreakdown, aes(x=DS, y=histvalmean,fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) +
  theme_bw() + 
  ggtitle(paste("1981-2005 ",titlepart,sep=""))+xlab("Downscaling Technique")+ylab("Historical Percentage")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed")
