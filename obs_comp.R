########
# Step 3 analysis 3^5
# Plot historical comparisons

#plot_indiv = function(step1_filename,projnotes,colorchoicediff,BINLIMIT,diffbartype,use_fixed_scale,fixed_scale){

###################################################################################################
## must pass in the following: step1_filename ,projnotes, colorchoicediff, BINLIMIT, diffbartype ##
###################################################################################################
filename = "/data2/3to5/I35/all_mems/tasmax_allmem_absolute_2006-2035_ann.nc"
colorchoicehist = "yellowtored"
colorchoicediff = "bluetored"
BINLIMIT = 30
diffbartype = "difference"
use_fixed_scale = FALSE
varin="tasmax"

library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}


histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
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

histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

#projnotes = do.call("c",strsplit(projnotes,",",fixed=TRUE))

split1 = strsplit(filename,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
difftype = split2[3]
seasonin = substr(split2[5],1,3)

###

test = nc_open(filename) #data for period 1
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
histvals = ncvar_get(test,"histmean")
varunits = test$var[[1]]$units
nc_close(test)

########
# Plotting Differences

if(all(lon>0)) lon =lon-360

if(use_fixed_scale==TRUE){
  print(fixed_scale)
  fixed_scale = do.call("c",strsplit(fixed_scale,",",fixed=TRUE))
  print(fixed_scale)
  if(1 %in% grep("n",fixed_scale)) fixed_scale[1] =  paste("-",substr(fixed_scale[1],2,nchar(fixed_scale[1])),sep="")
  if(2 %in% grep("n",fixed_scale)) fixed_scale[2] =  paste("-",substr(fixed_scale[2],2,nchar(fixed_scale[2])),sep="")
  print(fixed_scale)
  fixed_scale=as.numeric(fixed_scale)
  print(fixed_scale)
}

histvals = histvals[,,1:3]
obsnames = histfilebreakdown$obs[1:3]
histmean = apply(histvals,c(1,2),mean,na.rm=TRUE)
histcolorbar = colorramp(c(histvals,histmean),colorchoice=colorchoicehist,Blimit=30,type="raw")

diffs = histvals
for(i in 1:3) diffs[,,i] = histvals[,,i] - histmean

diffcolorbar = colorramp(diffs,colorchoice=colorchoicediff,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = use_fixed_scale,fixed_scale=fixed_scale)

#######################
# plot raw values and differences
plotfilename = paste("/data2/3to5/I35/plots/obs_plots/obsplots_",varname,"_",seasonin,".pdf",sep="")

pdf(plotfilename,onefile=TRUE,width=20,height=20)
par(mfrow=c(3,3),mar=c(2,1,2,3))

for(h in 1:3){
  testsfc1 = list(x=lon,y=lat,z=histvals[,,h])
  surface(testsfc1,type="I",zlim=histcolorbar[[1]],col=histcolorbar[[3]],breaks=histcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=FALSE)
  map("state",add=TRUE)
  text(-108.95,27,labels=obsnames[h],cex=2,pos = 4)
}

#plot(1,1,type="n",main="",ylab="",xlab="",xaxt='n', yaxt='n')

frame()

testsfc1 = list(x=lon,y=lat,z=histmean)
surface(testsfc1,type="I",zlim=histcolorbar[[1]],col=histcolorbar[[3]],breaks=histcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=FALSE)
map("state",add=TRUE) #main=paste("GCM: ",GCM," DS: ", DS," Obs: ",obs,sep="")
text(-108.95,27,labels="Mean",cex=2,pos = 4)

frame()

for(h in 1:3){
  testsfc1 = list(x=lon,y=lat,z=diffs[,,h])
  surface(testsfc1,type="I",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=FALSE)
  map("state",add=TRUE)
  text(-108.95,27,labels=paste(obsnames[h],"-Mean",sep=""),cex=2,pos = 4)
}

dev.off()
