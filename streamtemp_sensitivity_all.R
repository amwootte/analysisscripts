
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

gauges = c("boisdearc","charlies","cr340")

#####
# Read in climo files

climoprojall = NULL

for(i in 1:length(gauges)){
  
  filename=paste("/home/woot0002/CPREP_streamtemp_",gauges[i],"_2006-2099_fixed.Rdata",sep="")
  
  load(filename)
  
  climoprojall = rbind(climoprojall,projclimoanom)
  
}

#####

split1 = do.call("rbind",strsplit(as.character(climoprojall[,1]),"_",fixed=TRUE))
names(split1) = c("scen","GCM","DS","obs")

climoprojall = cbind(climoprojall,split1[,2:4])

names(climoprojall)[21:23] = c("GCM","DS","obs")

#####

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
                                            

climoprojin = subset(climoprojall,scen=="rcp85" & as.character(period)=="end-of-century")


library(ggplot2)

#####
# historical plots

varnames = c("TMAX","TMIN","tmaxfit5d","tminfit5d","tavgfit5d",
             "LT05gclimo","LT05jclimo")
 units = c(rep("C",5),"days/year","events/year")


for(i in 1:length(varnames)){
  
  projanomidx = which(names(climoprojin)==varnames[i])
  # hist plot
  
  if(i<=5){
    rangein = range(climoprojin[,c(2,3,6:8)],na.rm=TRUE)
  } else {
    rangein = range(climoprojin[,c(10,12)],na.rm=TRUE)
  }
  
  
  
  # future anomaly plot
  ggplot(climoprojin, aes(x=namein, y=climoprojin[,projanomidx], fill=GCM)) + 
    stat_boxplot(geom ='errorbar', width = 0.5) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=0.1,width=0.5,position=position_dodge2(width=0.5)) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=1) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_manual(values=c("brown", "cyan", "yellow"))+
    ggtitle(paste(varnames[i]," future anomaly - San Saba River",sep=""))+xlab("")+ylab(units[i])+ylim(rangein) 
  ggsave(paste("/home/woot0002/SanSaba_projanom_streamtmp_",varnames[i],".pdf",sep=""),device = "pdf",width=8,height=7)
  
}












