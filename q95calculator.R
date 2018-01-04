
library(ncdf4) # loading necessary libraries and extra functions
library(maps)
library(fields)
library(sp)
source("/home/woot0002/analysisfunctions.R")

########################################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!!!

###########
# 1. Data Gather and conversion

  histfilelist = system("ls /data2/3to5/I35/tasmax/EDQM/*historical*.nc",intern=T)
  tmaxhistfilelist = histfilelist
  tminhistfilelist = system("ls /data2/3to5/I35/tasmin/EDQM/*historical*.nc",intern=T)
 
filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),3)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  
  test1 = nc_open(tmaxhistfilelist[i])
  if(i==1){
    lon = ncvar_get(test1,"lon")
    lat = ncvar_get(test1,"lat")
    q95max = q95min= matrix(NA,nrow=length(lon),ncol=length(lat))
  }
  test2 = nc_open(tminhistfilelist[i])
  
  message("starting calcs for q95")
  
  for(R in 1:length(lon)){
    for(C in 1:length(lat)){
      tmp1 = ncvar_get(test1,"tasmax",start = c(R,C,1),count=c(1,1,-1))
      tmp2 = ncvar_get(test2,"tasmin",start = c(R,C,1),count=c(1,1,-1))
      if(all(is.na(tmp1)==TRUE)==FALSE){
        q95max[R,C] = quantile(tmp1,probs=0.95,na.rm=TRUE)
        q95min[R,C] = quantile(tmp2,probs=0.95,na.rm=TRUE)
      } 
    message("Finished Calcs for R ",R," and C ",C)
      }
  }
  
  nc_close(test1)
  nc_close(test2)
  
  filesplit = strsplit(tmaxhistfilelist[i],"/")[[1]]
  newq95maxfile = paste(substr(filesplit[length(filesplit)],1,(nchar(filesplit[length(filesplit)])-3)),"_q95.nc",sep="")
  filesplit = strsplit(tminhistfilelist[i],"/")[[1]]
  newq95minfile = paste(substr(filesplit[length(filesplit)],1,(nchar(filesplit[length(filesplit)])-3)),"_q95.nc",sep="")
  
  londef = ncdim_def("lon",units="degrees E",vals=lon,longname="Longitude")
  latdef = ncdim_def("lat",units="degrees N",vals=lat,longname="Latitude")
  
  q95maxdef = ncvar_def("tmaxq95",units="degrees K",dim=list(londef,latdef),missval=1E20,longname="95th Percentile of Daily High Temperatures")
  q95mindef = ncvar_def("tminq95",units="degrees K",dim=list(londef,latdef),missval=1E20,longname="95th Percentile of Daily Low Temperatures")
  
  nc1 = nc_create(newq95maxfile,q95maxdef)
  nc2 = nc_create(newq95minfile,q95mindef)
  
  ncvar_put(nc1,q95maxdef,q95max)
  ncvar_put(nc2,q95mindef,q95min)
  
  nc_close(nc1)
  nc_close(nc2)
  
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
