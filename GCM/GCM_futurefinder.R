source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "tasmin"
scen = "rcp85"
period=c(2070,2099)

startpoint = as.Date(paste(period[1],"-01-01",sep=""))
endpoint = as.Date(paste(period[2],"-12-31",sep=""))
datesout = seq(startpoint,endpoint,by="day")

GCMtable = read.table("/home/woot0002/GCMfilesforanalysis_LOCA.csv",sep=",",header=TRUE,colClasses = "character")
ALLGCMdat = read.table("https://climatedata.oscer.ou.edu/climatedata_updated-oscer_path.csv",sep=",",header=TRUE,colClasses = "character")

domain = unique(GCMtable$domain)
GCMexp = unique(paste(GCMtable$model,GCMtable$ensemble,sep="_"))

GCMsub= subset(ALLGCMdat,domain==domain & experiment==scen & variable==var)

filesout = NULL

for(i in 1:length(GCMexp)){
  split1=strsplit(GCMexp[i],"_",fixed=TRUE)[[1]]
  
  if(split1[1]!="GISS-E2-H" & split1[1]!="GISS-E2-R" & split1[1]!="EC-EARTH"){
    tmpsub = subset(GCMsub,model==split1[1] & ensemble==split1[2] & domain=="day"  & latest=="1.0")
  } else {
    
    if(split1[1]!="EC-EARTH"){
      tmpsub = subset(GCMsub,model==split1[1] & ensemble=="r2i1p1" & domain=="day"  & latest=="1.0")  
    } else {
      tmpsub = subset(GCMsub,model==split1[1] & ensemble=="r2i1p1" & domain=="day")  
    }
    
  }
  
  #tmpsub = subset(GCMsub,model==split1[1])
  if(length(unique(tmpsub$time))<nrow(tmpsub)){
    if(split1[1]=="HadGEM2-CC" | split1[1]=="CCSM4" | split1[1]=="GFDL-ESM2M"){
      if(split1[1]=="HadGEM2-CC") tmpsub=tmpsub[grep("v20120531",tmpsub$SCCASC_climatedata_path),]  
      if(split1[1]=="CCSM4") tmpsub=tmpsub[grep("v20130703",tmpsub$SCCASC_climatedata_path),] 
      if(split1[1]=="GFDL-ESM2M") tmpsub=tmpsub[grep("v20120228",tmpsub$SCCASC_climatedata_path),] 
    } else {
      tmpsub=tmpsub[grep("v4",tmpsub$SCCASC_climatedata_path),]   
    }
  }
  startdate = as.Date(substr(tmpsub$time,1,8),"%Y%m%d")
  enddate = as.Date(substr(tmpsub$time,10,18),"%Y%m%d")
  
  idxout=c()
  for(r in 1:nrow(tmpsub)){
    datestmp = seq(startdate[r],enddate[r],by="day")
    if( any(datestmp %in% datesout)==TRUE){
      idxout=c(idxout,r)
    }
  }
  
  filesout = rbind(filesout,tmpsub[idxout,])
}

write.table(filesout,"/home/woot0002/tasmin_GCMfilelist_LOCA.csv",sep=",",row.names=FALSE)
