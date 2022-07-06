source("/data2/3to5/I35/scripts/analysisfunctions.R")
source("/data2/3to5/I35/scripts/springpheno_v0.5.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(multiApply)

filesin = system("ls /data4/data/DS_proj/LOCA/pr/future/*.nc",intern=TRUE)

filesplit = do.call("rbind",strsplit(filesin,split="/",fixed=TRUE))
filesplit2 = do.call("rbind",strsplit(filesplit[,ncol(filesplit)],"_",fixed=TRUE))
filesplit3 = do.call("rbind",strsplit(filesplit2[,ncol(filesplit2)],"-",fixed=TRUE))
filesplit3[,ncol(filesplit3)] = substr(filesplit3[,ncol(filesplit3)],1,4)

filetab = cbind(filesplit,filesplit2)
filetab = cbind(filetab,filesplit3)

filetab = data.frame(filetab)
names(filetab) = c("blank","dir1","dir2","dir3","dir4","dir5","dir6","filename","var","GCM","exp","group","scen","startyear","endyear")
GCMs = unique(as.character(filetab$GCM))
scens = unique(as.character(filetab$scen))

filetab$varin = paste(filetab$var,filetab$GCM,filetab$exp,filetab$scen,sep="_")

datesall = seq(as.Date("2006-01-01"),as.Date("2100-12-31"),by="day")
yearmonall = unique(substr(datesall,1,7))

for(i in 3:length(GCMs)){
  for(j in 1:length(scens)){
    
    filesused = filesin[which(filetab$GCM==GCMs[i] & filetab$scen==scens[j])]
    filetabused = filetab[which(filetab$GCM==GCMs[i] & filetab$scen==scens[j]),]
    
    timeoutmon = c()
    
    for(f in 1:length(filesused)){
      
      nctest = nc_open(filesused[f])
      
      if(i==1 & j==1 & f==1){
        dataunits = nctest$var[[1]]$units
        lon = ncvar_get(nctest,"lon")
        lat = ncvar_get(nctest,"lat")
        timeunits = nctest$var[[1]]$dim[[3]]$units
        latunits = nctest$var[[1]]$dim[[2]]$units
        lonunits = nctest$var[[1]]$dim[[1]]$units
        
        timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
        indexdate = as.Date(substr(timeunits,12,21))
      }  
        
      time=ncvar_get(nctest,"time")
        dates=seq(as.Date(paste(filetabused$startyear[f],"-01-01",sep="")),as.Date(paste(filetabused$endyear[f],"-12-31",sep="")),by="day")
        timeout= time
        yearmon = unique(substr(dates,1,7))
        timeoutmon = c(timeoutmon,timeout[which(substr(dates,9,10)=="01")])
      
      if(f==1){
        tmp = array(NA,dim=c(length(lon),length(lat),length(yearmonall)))
      }
        
        if(length(yearmon)>120){
          yearmon=yearmon[1:120]
        }
        
        for(y in 1:length(yearmon)){
          idxmon = which(substr(dates,1,7)==yearmon[y])
          idxout = which(yearmonall == yearmon[y])
          message("idxout is ", idxout)
          startpoint = idxmon[1]
          endpoint = length(idxmon)
          if(idxmon[length(idxmon)]>length(time)){
            endpoint = length(idxmon[1:which(idxmon==length(time))])-1
          }
          testdat = ncvar_get(nctest,filetabused$varin[f],start=c(1,1,startpoint),count=c(-1,-1,endpoint))
          #if(mean(testdat,na.rm=TRUE)<1){
            testdat = testdat*86400
          #}
          
          tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)  
          tmp[,,idxout]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
          rm(testdat)
          gc()
          message("Finished with yearmon ",y," / ",length(yearmon))
        }
        nc_close(nctest)
    }
    
    filenamepart = paste(filetabused$var[1],"_","mon","_",filetabused$GCM[1],"_",filetabused$exp[1],"_",filetabused$scen[1],".nc",sep="")
    filenameout = paste("/data4/data/DS_proj/LOCA/prmon/future/",filenamepart,sep="")
    
    dimX <- ncdim_def( "lon", lonunits, lon)
    dimY <- ncdim_def( "lat", latunits, lat)
    dimT <- ncdim_def("time",timeunits,timeoutmon)
    
    # Make varables of various dimensionality, for illustration purposes
    mv <- 1E20 # missing value to use
    
    var1d <- ncvar_def(filetabused$varin[1],dataunits,longname=filetabused$varin[1], list(dimX,dimY,dimT), mv ,compression=9)
    
    #######
    # Create netcdf file
    
    nc <- nc_create(filenameout ,  var1d )
    # Write some data to the file
    ncvar_put(nc, var1d, tmp) # no start or count: write all values\
    # close ncdf
    nc_close(nc)
    rm(tmp)
    gc()
    message("Finished calcs for GCM ",i," / ",length(GCMs)," and for scen ",j," / ",length(scens))
    
  }
}
  
 
 

