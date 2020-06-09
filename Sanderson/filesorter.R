
source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

load("/home/woot0002/DS_ind/GCMlist.Rdata")

for(i in 1:length(GCMlist)){
  
  GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/*",GCMlist[i],"_*histclimo*.nc",sep=""),intern=TRUE)
  LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/*",GCMlist[i],"_*histclimo*.nc",sep=""),intern=TRUE)
  
  for(f in 1:length(GCMfiles)){
    system(paste("cp ",GCMfiles[f]," /data/static_web/RCMES/CMIP5/historical/.",sep=""))
  }
  
  for(f in 1:length(LOCAfiles)){
    system(paste("cp ",LOCAfiles[f]," /data/static_web/RCMES/LOCA/historical/.",sep=""))
  }
  
}

for(i in 1:length(GCMlist)){
  
  GCMfiles = system(paste("ls /home/woot0002/GCMs/regrid/*",GCMlist[i],"_*projclimo*.nc",sep=""),intern=TRUE)
  LOCAfiles = system(paste("ls /home/woot0002/LOCA/regrid/*",GCMlist[i],"_*projclimo*.nc",sep=""),intern=TRUE)
  
  for(f in 1:length(GCMfiles)){
    system(paste("cp ",GCMfiles[f]," /data/static_web/RCMES/CMIP5/future/.",sep=""))
  }
  
  for(f in 1:length(LOCAfiles)){
    system(paste("cp ",LOCAfiles[f]," /data/static_web/RCMES/LOCA/future/.",sep=""))
  }
  
}



LIVNEHfiles = system("ls /home/woot0002/monthlyclimo/*_day*livneh*.nc",intern=TRUE)
for(f in 1:length(LIVNEHfiles)){
  system(paste("cp ",LIVNEHfiles[f]," /data/static_web/RCMES/OBS/.",sep=""))
}




