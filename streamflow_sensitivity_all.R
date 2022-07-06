
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

gauges = c("08144500","08144600","08146000")

#####
# Read in climo files

climoprojall = NULL

for(i in 1:length(gauges)){
  
  filename=paste("/home/woot0002/climotable_",gauges[i],"_2070-2099.csv",sep="")
  climoproj = read.table(filename,header=TRUE,sep=",",colClasses = c(rep("character",8),rep("numeric",67)))
  
  if(gauges[i]=="08144500"){
    gaugename = "San Saba at Menard"  
  }
  
  if(gauges[i]=="08144600"){
    gaugename = "San Saba near Brady"  
  }
  
  if(gauges[i]=="08146000"){
    gaugename = "San Saba near San Saba"  
  }
  
  climoproj$gauge = gauges[i]
  climoproj$gaugename = gaugename
  
  climoprojall = rbind(climoprojall,climoproj)
  
}

#####
# pull historical

climohistall = NULL

for(i in 1:length(gauges)){
  
  filename=paste("/home/woot0002/climotable_",gauges[i],"_1981-2005.csv",sep="")
  climohist = read.table(filename,header=TRUE,sep=",",colClasses = c(rep("character",8),rep("numeric",44))) #,colClasses = c(rep("character",8),rep("numeric",44))
  
  if(gauges[i]=="08144500"){
    gaugename = "San Saba at Menard"  
  }
  
  if(gauges[i]=="08144600"){
    gaugename = "San Saba near Brady"  
  }
  
  if(gauges[i]=="08146000"){
    gaugename = "San Saba near San Saba"  
  }
  
  climohist$gauge = gauges[i]
  climohist$gaugename = gaugename
  
  climohistall = rbind(climohistall,climohist)
  
}

#####

climoprojall$DS[which(climoprojall$DS=="QDM")] = "EDQM"

climoprojall$namein = paste(climoprojall$GCM,climoprojall$DS,climoprojall$obs,sep="_")

climoprojall$namein <- factor(climoprojall$namein,
                                 levels = c('CCSM4_DeltaSD_Daymet','MIROC5_DeltaSD_Daymet','MPI-ESM-LR_DeltaSD_Daymet',
                                            'CCSM4_EDQM_Daymet','MIROC5_EDQM_Daymet','MPI-ESM-LR_EDQM_Daymet',
                                            'CCSM4_PARM_Daymet','MIROC5_PARM_Daymet','MPI-ESM-LR_PARM_Daymet',
                                            'CCSM4_DeltaSD_Livneh','MIROC5_DeltaSD_Livneh','MPI-ESM-LR_DeltaSD_Livneh',
                                            'CCSM4_EDQM_Livneh','MIROC5_EDQM_Livneh','MPI-ESM-LR_EDQM_Livneh',
                                            'CCSM4_PARM_Livneh','MIROC5_PARM_Livneh','MPI-ESM-LR_PARM_Livneh',
                                            'CCSM4_DeltaSD_PRISM','MIROC5_DeltaSD_PRISM','MPI-ESM-LR_DeltaSD_PRISM',
                                            'CCSM4_EDQM_PRISM','MIROC5_EDQM_PRISM','MPI-ESM-LR_EDQM_PRISM',
                                            'CCSM4_PARM_PRISM','MIROC5_PARM_PRISM','MPI-ESM-LR_PARM_PRISM'),ordered = TRUE)
                                            

climoprojin = subset(climoprojall,scen=="rcp85")


climohistall$DS[which(climohistall$DS=="QDM")] = "EDQM"

climohistall$namein = paste(climohistall$GCM,climohistall$DS,climohistall$obs,sep="_")

climohistall$namein <- factor(climohistall$namein,
                              levels = c('CCSM4_DeltaSD_Daymet','MIROC5_DeltaSD_Daymet','MPI-ESM-LR_DeltaSD_Daymet',
                                         'CCSM4_EDQM_Daymet','MIROC5_EDQM_Daymet','MPI-ESM-LR_EDQM_Daymet',
                                         'CCSM4_PARM_Daymet','MIROC5_PARM_Daymet','MPI-ESM-LR_PARM_Daymet',
                                         'CCSM4_DeltaSD_Livneh','MIROC5_DeltaSD_Livneh','MPI-ESM-LR_DeltaSD_Livneh',
                                         'CCSM4_EDQM_Livneh','MIROC5_EDQM_Livneh','MPI-ESM-LR_EDQM_Livneh',
                                         'CCSM4_PARM_Livneh','MIROC5_PARM_Livneh','MPI-ESM-LR_PARM_Livneh',
                                         'CCSM4_DeltaSD_PRISM','MIROC5_DeltaSD_PRISM','MPI-ESM-LR_DeltaSD_PRISM',
                                         'CCSM4_EDQM_PRISM','MIROC5_EDQM_PRISM','MPI-ESM-LR_EDQM_PRISM',
                                         'CCSM4_PARM_PRISM','MIROC5_PARM_PRISM','MPI-ESM-LR_PARM_PRISM'),ordered = TRUE)



library(ggplot2)

#####
# historical plots

varnames = c("HOT1S_1W","HOT1S_3M","HOT90_12M","HOT90_24M","HOT90_3M",
             "MAXPRCP1W_12M","MAXPRCP1W_2W","MAXPRCP1W_6M",
             "MAXTMAX1W_2W","MAXTMAX1W_3W",
             "MDRN_3M","MDRN_4W","PRCP","PRCPAVG_12M","PRCPAVG_3D","PRCPAVG_6M",
             "SF","TAVG_12M","TAVG_3M","TMAX","TMAXAVG_1W","TMIN")
 units = c(rep("days",5),rep("in.",3),rep("deg F",2),rep("days",2),rep("in.",4),"cfs",rep("deg F",5))


for(i in 1:length(varnames)){
  
  histidx = which(names(climohistall)==varnames[i])
  projidx = which(names(climoprojin)==varnames[i])
  projanomidx = which(names(climoprojin)==paste(varnames[i],"anom",sep="_"))
  # hist plot
  
  ggplot(climohistall, aes(x=namein, y=climohistall[,histidx], fill=GCM)) + 
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
    ggtitle(paste(varnames[i]," hist - San Saba River",sep=""))+xlab("")+ylab(units[i])+ylim(range(climohistall[,histidx],na.rm=TRUE)) 
  ggsave(paste("/home/woot0002/SanSaba_hist_",varnames[i],".pdf",sep=""),device = "pdf",width=8,height=7)
  
  # future plot
  ggplot(climoprojin, aes(x=namein, y=climoprojin[,projidx], fill=GCM)) + 
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
    ggtitle(paste(varnames[i]," future - San Saba River",sep=""))+xlab("")+ylab(units[i])+ylim(range(climoprojin[,projidx],na.rm=TRUE)) 
  ggsave(paste("/home/woot0002/SanSaba_proj_",varnames[i],".pdf",sep=""),device = "pdf",width=8,height=7)
  
  # future anomaly plot
  ggplot(climoprojin, aes(x=namein, y=climoprojin[,projanomidx], fill=GCM)) + 
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
    ggtitle(paste(varnames[i]," future anomaly - San Saba River",sep=""))+xlab("")+ylab(units[i])+ylim(range(climoprojin[,projanomidx],na.rm=TRUE)) 
  ggsave(paste("/home/woot0002/SanSaba_projanom_",varnames[i],".pdf",sep=""),device = "pdf",width=8,height=7)
  
}










ggplot(climoprojin, aes(x=namein, y=TMAX, fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
  ggtitle("tmax Future - San Saba River")+xlab("")+ylab("TMAX (F)")+ylim(range(climoprojin$TMAX,na.rm=TRUE)) 

ggplot(climoprojin, aes(x=namein, y=TMIN, fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
  ggtitle("tmin Future - San Saba River")+xlab("")+ylab("TMIN (F)")+ylim(range(climoprojin$TMIN,na.rm=TRUE)) 



ggplot(climoprojin, aes(x=namein, y=R1MM, fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
  ggtitle("r1mm Future - San Saba River")+xlab("")+ylab("R1mm (days)")+ylim(range(climoprojin$R1MM,na.rm=TRUE)) 

ggplot(climoprojin, aes(x=namein, y=PRCP, fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
  ggtitle("prcp Future - San Saba River")+xlab("")+ylab("PRCP (in)")+ylim(range(climoprojin$PRCP,na.rm=TRUE)) 

ggplot(climoprojin, aes(x=namein, y=SF, fill=GCM)) + 
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
  ggtitle("Streamflow Future - San Saba River")+xlab("")+ylab("streamflow (cfs)")+ylim(range(climoprojin$SF,na.rm=TRUE)) 


