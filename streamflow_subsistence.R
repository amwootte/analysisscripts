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

gaugenum = "08144500"
futureperiod = c(2070,2099)

streamflow = read.table(paste("/home/woot0002/streamflow_",gaugenum,sep=""),header=TRUE,sep="\t",fill=TRUE) # watch it with reordered tables based on gauge ID

if(gaugenum!="08151500" & gaugenum!="08144500" & gaugenum!="08148500"){
  names(streamflow) = c("agency","site","DATE","streamflow_mean","streamflow_mean_QC","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC")
  streamidx = 4
} else {
  names(streamflow) = c("agency","site","DATE","streamflow_max","streamflow_max_QC","streamflow_min","streamflow_min_QC","streamflow_mean","streamflow_mean_QC")
  streamidx = 8
}

tmp = c(NA,rollapply(streamflow$streamflow_mean,3,mean,na.rm=TRUE),NA)
streamflow$streamflow_mean=tmp

if(gaugenum=="08144500"){
  subnums = c(2.1,1.5,1,1)
}

if(gaugenum=="08144600"){
  subnums = c(20,6.8,1,1)
}

if(gaugenum=="08146000"){
  subnums = c(29,22,3,13)
}

if(gaugenum=="08151500"){
  subnums = c(44,35,3,20)
}

####
# historical values

streamflow2 = subset(streamflow,as.numeric(substr(streamflow$DATE,1,4))>=1981 & as.numeric(substr(streamflow$DATE,1,4))<=2005)

yearlycalcs = aggregate(streamflow2$streamflow_mean,by=list(year=substr(streamflow2$DATE,1,4)),mean,na.rm=TRUE)
mean(yearlycalcs[,2],na.rm=TRUE)

histvals = c()

for(s in 1:4){
  
  if(s==1){
    streamsub = subset(streamflow2,as.numeric(substr(streamflow2$DATE,6,7))>=12 | as.numeric(substr(streamflow2$DATE,6,7))<=2)  
  }
  if(s==2){
    streamsub = subset(streamflow2,as.numeric(substr(streamflow2$DATE,6,7))>=3 & as.numeric(substr(streamflow2$DATE,6,7))<=5)  
  }
  if(s==3){
    streamsub = subset(streamflow2,as.numeric(substr(streamflow2$DATE,6,7))>=6 & as.numeric(substr(streamflow2$DATE,6,7))<=8)  
  }
  if(s==4){
    streamsub = subset(streamflow2,as.numeric(substr(streamflow2$DATE,6,7))>=9 & as.numeric(substr(streamflow2$DATE,6,7))<=11)  
  }
  
  streamsub$subday = ifelse(streamsub$streamflow_mean<=subnums[s],1,0)
  
  yearlycalcs = aggregate(streamsub$subday,by=list(year=substr(streamsub$DATE,1,4)),sum,na.rm=TRUE)
  histvals[s]=mean(yearlycalcs[,2],na.rm=TRUE)
  
}




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

########
# Figure out obs percentile for subsistence flows

#streamflow = subset(streamflow,as.numeric(substr(DATE,1,4))>=1981 & as.numeric(substr(DATE,1,4))<=2005)

subpercs = c()

for(s in 1:4){
  
  if(s==1){
    idx = which((as.numeric(substr(streamflow$DATE,6,7))>=12 | as.numeric(substr(streamflow$DATE,6,7))<=2) & streamflow$streamflow_mean>0)
  }
  if(s==2){
    idx = which(as.numeric(substr(streamflow$DATE,6,7))>=3 & as.numeric(substr(streamflow$DATE,6,7))<=5 & streamflow$streamflow_mean>0)
  }
  if(s==3){
    idx = which(as.numeric(substr(streamflow$DATE,6,7))>=6 & as.numeric(substr(streamflow$DATE,6,7))<=8 & streamflow$streamflow_mean>0)
  }
  if(s==4){
    idx = which(as.numeric(substr(streamflow$DATE,6,7))>=9 & as.numeric(substr(streamflow$DATE,6,7))<=11 & streamflow$streamflow_mean>0)
  }
  
  streamtmp = streamflow$streamflow_mean[idx]
  
  quants = seq(0.000001,0.999999,by=0.000001)
  streamquants = quantile(streamtmp,probs=quants,na.rm=TRUE)
  
  diffstreamquants = streamquants-subnums[s]
  quantidx = which(diffstreamquants==0)
  
  if(length(quantidx)>=1){
    subpercs[s] = quants[quantidx[length(quantidx)]]
  } else {
    minval = min(abs(diffstreamquants))
    foundit = 0
    q=1
    while(foundit==0){
      dq = diffstreamquants[q]
      
      if(abs(dq)==minval){
        foundit=1
      }
      q=q+1
    }
    subpercs[s]=quants[q]
  }
  
}

########
# Figure out matching percentiles numbers and number of days below for each.

subvals = matrix(NA,nrow=length(histresults),ncol=4)
numdayslowsub = matrix(NA,nrow=length(histresults),ncol=4)

histdates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:27){
  
  histtmp = histresults[[i]][[1]]
  
  if(length(histtmp)<length(histdates)){
    datesin = histdates[-which(substr(histdates,6,10)=="02-29")]
  } else {
    datesin = histdates
  }
  
  for(s in 1:4){
    if(s==1){
      idx = which((as.numeric(substr(datesin,6,7))>=12 | as.numeric(substr(datesin,6,7))<=2) & histtmp>0)
    }
    if(s==2){
      idx = which(as.numeric(substr(datesin,6,7))>=3 & as.numeric(substr(datesin,6,7))<=5 & histtmp>0)
    }
    if(s==3){
      idx = which(as.numeric(substr(datesin,6,7))>=6 & as.numeric(substr(datesin,6,7))<=8 & histtmp>0)
    }
    if(s==4){
      idx = which(as.numeric(substr(datesin,6,7))>=9 & as.numeric(substr(datesin,6,7))<=11 & histtmp>0)
    }
    
    streamtmp = histtmp[idx]
    subvals[i,s] = quantile(streamtmp,probs=subpercs[s],na.rm=TRUE)
    
    numdayslowsub[i,s] = length(which(streamtmp<subvals[i,s]))/25
    
  }
}

###########
# find number of historical and projected percentile events

########
# Figure out matching percentiles numbers and number of days below for each.

numdayslowsub_proj = matrix(NA,nrow=length(projresults),ncol=4)
numdayslowsub_anom = matrix(NA,nrow=length(projresults),ncol=4)

projdates = seq(as.Date(paste(futureperiod[1],"-01-01",sep="")),as.Date(paste(futureperiod[2],"-12-31",sep="")),by="day")

for(i in 1:81){
  
  projtmp = projresults[[i]][[1]]
  
  if(length(projtmp)<length(projdates)){
    datesin = projdates[-which(substr(projdates,6,10)=="02-29")]
  } else {
    datesin = projdates
  }
  
  histidx = which(histfilebreakdown$GCM == projfilebreakdown$GCM[i] & histfilebreakdown$DS==projfilebreakdown$DS[i] & histfilebreakdown$obs==projfilebreakdown$obs[i])
  
  for(s in 1:4){
    if(s==1){
      idx = which(as.numeric(substr(datesin,6,7))>=12 | as.numeric(substr(datesin,6,7))<=2)
    }
    if(s==2){
      idx = which(as.numeric(substr(datesin,6,7))>=3 & as.numeric(substr(datesin,6,7))<=5)
    }
    if(s==3){
      idx = which(as.numeric(substr(datesin,6,7))>=6 & as.numeric(substr(datesin,6,7))<=8)
    }
    if(s==4){
      idx = which(as.numeric(substr(datesin,6,7))>=9 & as.numeric(substr(datesin,6,7))<=11)
    }
    
    streamtmp = projtmp[idx]
    
    numdayslowsub_proj[i,s] = length(which(streamtmp<subvals[histidx,s]))/30
    numdayslowsub_anom[i,s] = numdayslowsub_proj[i,s]-numdayslowsub[histidx,s]
    
  }
}

########
# Create tables

pfilebreakdown = NULL

for(s in 1:4){
  
  numlowsubflow = numdayslowsub_anom[,s]
  
  tmp = projfilebreakdown
  
  if(s==1){
    tmp$season = "Winter"
  } 
  if(s==2){
    tmp$season = "Spring"
  } 
  if(s==3){
    tmp$season = "Summer"
  } 
  if(s==4){
    tmp$season = "Fall"
  } 
  tmp$numlowsubflow = numlowsubflow
  pfilebreakdown = rbind(pfilebreakdown,tmp)
  
}

pfilebreakdown$season = factor(pfilebreakdown$season,levels=c("Winter","Spring","Summer","Fall"))


aggregate(pfilebreakdown$numlowsubflow,by=list(scen_sea = paste(pfilebreakdown$scen,pfilebreakdown$season,sep="_")),mean,na.rm=TRUE)

aggregate(pfilebreakdown$numlowsubflow,by=list(scen_sea = paste(pfilebreakdown$scen,pfilebreakdown$season,sep="_")),sd,na.rm=TRUE)


####
# boxplots

library(ggplot2)

ggplot(pfilebreakdown, aes(x=season, y=numlowsubflow,fill=scen)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=0.5,width=0.5,position=position_dodge2(width=0.5)) + scale_fill_manual(values=c( "green", "blue","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + geom_hline(yintercept=0,linetype="dashed") +
  ggtitle(paste("Projected Change (",futureperiod[1],"-",futureperiod[2],")\nNumber of Events/year below Subsistence ",sep=""))+xlab("")+ylab("Number of Events") + ylim(0,80)
ggsave(paste("/home/woot0002/streamflow_projchangesubsist_",gaugenum,".pdf",sep=""),device = "pdf",width=5,height=5)





