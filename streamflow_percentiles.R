#####
#
# streamflow testing - regression with C-PrEP

library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)
library(zoo)
library(lars)
library(mailR)

gaugenum = "08146000"
futureperiod = c(2040,2069)

# load DS files

prfile_hist = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_*historical*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_*historical*.nc",intern=TRUE))

tasmaxfile_hist = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_*historical*.nc",intern=TRUE))

tasminfile_hist = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_*historical*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_*historical*.nc",intern=TRUE))

prfile_proj = c(system("ls /data2/3to5/I35/pr/DeltaSD/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/EDQM/pr_day*rcp*.nc",intern=TRUE),
                system("ls /data2/3to5/I35/pr/PARM/pr_day*rcp*.nc",intern=TRUE))

tasmaxfile_proj = c(system("ls /data2/3to5/I35/tasmax/DeltaSD/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/EDQM/tasmax_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmax/PARM/tasmax_day*rcp*.nc",intern=TRUE))

tasminfile_proj = c(system("ls /data2/3to5/I35/tasmin/DeltaSD/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/EDQM/tasmin_day*rcp*.nc",intern=TRUE),
                    system("ls /data2/3to5/I35/tasmin/PARM/tasmin_day*rcp*.nc",intern=TRUE))

#####
# Create file breakdown tables

filebreakdown = do.call(rbind,strsplit(prfile_proj,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(prfile_hist,"_",fixed=TRUE))
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

################
# load data

load(file=paste("/home/woot0002/CPREP_",gaugenum,"_",futureperiod[1],"-",futureperiod[2],"_fixed.Rdata",sep=""))

################
# Get historical quantiles

histfilebreakdown$histlow1 = NA
histfilebreakdown$histlow5 = NA
histfilebreakdown$histhigh1 = NA
histfilebreakdown$histhigh5 = NA

for(i in 1:length(histresults)){
  histfilebreakdown$histlow1[i] = quantile(histresults[[i]][[1]],probs=0.01,na.rm=TRUE)
  histfilebreakdown$histlow5[i] = quantile(histresults[[i]][[1]],probs=0.05,na.rm=TRUE)
  histfilebreakdown$histhigh5[i] = quantile(histresults[[i]][[1]],probs=0.95,na.rm=TRUE)
  histfilebreakdown$histhigh1[i] = quantile(histresults[[i]][[1]],probs=0.99,na.rm=TRUE)
}

################
# projected changes and comparison with obs

streamflow = read.table(paste("/home/woot0002/streamflow_",gaugenum,sep=""),header=TRUE,sep="\t",fill=TRUE)

if(gaugenum!="08151500" & gaugenum!="08144500" & gaugenum!="08148500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}

streamflow$DATE = as.character(streamflow$DATE)
streamflow$streamflow_mean=as.numeric(streamflow$streamflow_mean)
tmp = c(NA,rollapply(streamflow$streamflow_mean,3,mean,na.rm=TRUE),NA)
streamflow$streamflow_mean=tmp

streamflow$YEAR = as.numeric(substr(streamflow$DATE,1,4))
streamflow = subset(streamflow,YEAR>=1981 & YEAR<=2005)

OBSlow1 = quantile(streamflow$streamflow_mean,probs=0.01,na.rm=TRUE)
OBSlow5 = quantile(streamflow$streamflow_mean,probs=0.05,na.rm=TRUE)
OBShigh1 = quantile(streamflow$streamflow_mean,probs=0.99,na.rm=TRUE)
OBShigh5 = quantile(streamflow$streamflow_mean,probs=0.95,na.rm=TRUE)

boxplot(histfilebreakdown$histlow5,ylab="5th Percentile of Daily Streamflow (ft^3/s)",main="DS 5th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$histlow5,OBSlow5),na.rm=TRUE))
points(OBSlow5,pch=19,col="blue")

boxplot(histfilebreakdown$histhigh5,ylab="95th Percentile of Daily Streamflow (ft^3/s)",main="DS 95th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$histhigh5,OBShigh5),na.rm=TRUE))
points(OBShigh5,pch=19,col="blue")

boxplot(histfilebreakdown$histhigh1,ylab="99th Percentile of Daily Streamflow (ft^3/s)",main="DS 99th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$histhigh1,OBShigh1),na.rm=TRUE))
points(OBShigh1,pch=19,col="blue")

###########
# find number of historical and projected percentile events

# combined historical daily obs and modeled

num_OBSlow1 = length(which(streamflow$streamflow_mean<=OBSlow1))/25
num_OBSlow5 = length(which(streamflow$streamflow_mean<=OBSlow5))/25
num_OBShigh1 = length(which(streamflow$streamflow_mean>=OBShigh1))/25
num_OBShigh5 = length(which(streamflow$streamflow_mean>=OBShigh5))/25

histfilebreakdown$num_histlow1 = NA
histfilebreakdown$num_histlow5 = NA
histfilebreakdown$num_histhigh1 = NA
histfilebreakdown$num_histhigh5 = NA

for(i in 1:nrow(histfilebreakdown)){
  histfilebreakdown$num_histlow1[i] = length(which(histresults[[i]][[1]]<=histfilebreakdown$histlow1[i]))/25
  histfilebreakdown$num_histlow5[i] = length(which(histresults[[i]][[1]]<=histfilebreakdown$histlow5[i]))/25
  histfilebreakdown$num_histhigh1[i] = length(which(histresults[[i]][[1]]>=histfilebreakdown$histhigh1[i]))/25
  histfilebreakdown$num_histhigh5[i] = length(which(histresults[[i]][[1]]>=histfilebreakdown$histhigh5[i]))/25
}

boxplot(histfilebreakdown$num_histlow5,ylab="Number of Events",main="Number of Events \n DS 5th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$num_histlow5,num_OBSlow5),na.rm=TRUE))
points(num_OBSlow5,pch=19,col="blue")

boxplot(histfilebreakdown$num_histhigh5,ylab="Number of Events",main="Number of Events \n DS 95th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$num_histhigh5,num_OBShigh5),na.rm=TRUE))
points(num_OBShigh5,pch=19,col="blue")

boxplot(histfilebreakdown$num_histhigh1,ylab="Number of Events",main="Number of Events \n DS 99th percentile of streamflow vs. Obs",ylim=range(c(histfilebreakdown$num_histhigh1,num_OBShigh1),na.rm=TRUE))
points(num_OBShigh1,pch=19,col="blue")


projfilebreakdown$num_histlow1 = NA
projfilebreakdown$num_histlow5 = NA
projfilebreakdown$num_histhigh1 = NA
projfilebreakdown$num_histhigh5 = NA

for(i in 1:nrow(projfilebreakdown)){
  
  OBSin = projfilebreakdown$obs[i]
  DSin = projfilebreakdown$DS[i]
  GCMin = projfilebreakdown$GCM[i]
  
  histidx = which(histfilebreakdown$obs==OBSin & histfilebreakdown$DS==DSin & histfilebreakdown$GCM==GCMin)
  
  tmp = as.vector(projresults[[i]])[[1]]
  
  num_hlow1 = length(which(tmp<=histfilebreakdown$histlow1[histidx]))/length(futureperiod[1]:futureperiod[2])
  num_hlow5 = length(which(tmp<=histfilebreakdown$histlow5[histidx]))/length(futureperiod[1]:futureperiod[2])
  num_hhigh1 = length(which(tmp>=histfilebreakdown$histhigh1[histidx]))/length(futureperiod[1]:futureperiod[2])
  num_hhigh5 = length(which(tmp>=histfilebreakdown$histhigh5[histidx]))/length(futureperiod[1]:futureperiod[2])
  
  projfilebreakdown$num_histlow1[i] = num_hlow1-histfilebreakdown$num_histlow1[histidx]
  projfilebreakdown$num_histlow5[i] = num_hlow5-histfilebreakdown$num_histlow5[histidx]
  projfilebreakdown$num_histhigh1[i] = num_hhigh1-histfilebreakdown$num_histhigh1[histidx]
  projfilebreakdown$num_histhigh5[i] = num_hhigh5-histfilebreakdown$num_histhigh5[histidx]
  
}


####
# boxplots

library(ggplot2)

ggplot(projfilebreakdown, aes(x=scen, y=num_histlow5,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change (",futureperiod[1],"-",futureperiod[2],")\nNumber of Events/year below 5th percentile ",sep=""))+xlab("")+ylab("Number of Events") 

ggplot(projfilebreakdown, aes(x=scen, y=num_histhigh5,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change (",futureperiod[1],"-",futureperiod[2],")\nNumber of Events/year above 95th percentile ",sep=""))+xlab("")+ylab("Number of Events") 

ggplot(projfilebreakdown, aes(x=scen, y=num_histhigh1,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change (",futureperiod[1],"-",futureperiod[2],")\nNumber of Events/year above 99th percentile ",sep=""))+xlab("")+ylab("Number of Events") 





