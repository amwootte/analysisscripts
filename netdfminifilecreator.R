############################
#
# NetCDF mini file creator
#
#############################

library(ncdf4)
library(sp)
library(fields)

varname = "pr" # short name for the variable of interest options include tasmax, tasmin, pr, tmax95, tmax100, tmin32, tmin28, pr25, and pr50
varin = varname # don't change this

unitsin="mm"

projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/EDQM/*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

years = 2006:2025
prdat = list()

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")
leapdays = which(substr(dates,6,10)=="02-29")
dates = dates[-leapdays]

timeidx = which(as.numeric(substr(dates,1,4))>=years[1] & as.numeric(substr(dates,1,4))<=years[length(years)])

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  test = nc_open(projfilelist[i])
  tempdata = ncvar_get(test,"pr",start = c(1,1,timeidx[1]),count=c(-1,-1,length(timeidx)))
  
  if(i==1){
    lat = ncvar_get(test,"lat")
    lon = ncvar_get(test,"lon")-360
    times = ncvar_get(test,"time")
    times = times[timeidx]
  }
  nc_close(test)
  
  tempdata=tempdata*86400
  
  newfilename = paste(paste(varname,projfilebreakdown$DS[i],projfilebreakdown$obs[i],projfilebreakdown$scen[i],projfilebreakdown$GCM[i],years[1],years[2],sep="_"),".nc",sep="")
  
  if(i==1){
    dimX <- ncdim_def("lon","degrees_east",lon)
    dimY <- ncdim_def("lat","degrees_north",lat)
    dimT <- ncdim_def("times","days since 1961-01-01 00:00",times)
  }
  var1 <- ncvar_def(varname,unitsin,dim=list(dimX,dimY,dimT),missval=999999)
  nc <- nc_create(paste("/home/woot0002/miniNetCDFS/",newfilename,sep=""),var1)
  ncvar_put(nc,var1,tempdata)
  nc_close(nc)
  
  rm(tempdata)
  rm(var1)
  rm(nc)
  gc()
  gc()
  gc()
  gc()
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to trim file: ",ptmend[3]-ptm[3]," secs")
}



