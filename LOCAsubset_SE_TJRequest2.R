source("/data2/3to5/I35/scripts/analysisfunctions.R")
source("/data2/3to5/I35/scripts/springpheno_v0.5.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

var = "tasmin"
tempres = "annual"
varout = "frd"

filesin = system(paste("ls /data4/data/DS_proj/LOCA/",var,"/historical/*.nc",sep=""),intern=TRUE)

#if(var=="tasmax" & varout!="heatwaves"){
#  filesin = filesin[-grep("MRI-CGCM3",filesin)]
#}

if(var=="tasmin" & varout!="heatwaves"){
  filesin = filesin[-grep("HadGEM2-ES",filesin)]
}

if(varout=="heatwaves" | varout=="FLI" | varout=="FBI"){
  #filesin = filesin[-grep("MRI-CGCM3",filesin)]
  filesin = filesin[-grep("HadGEM2-ES",filesin)]
}

latbnds = c(25,37)
lonbnds = c(-107.90,-88.770)

for(i in 1:length(filesin)){
  
  ptm = proc.time()
  
  filesplit = strsplit(filesin[i],split="/",fixed=TRUE)
  filesplit2 = strsplit(filesplit[[length(filesplit)]],"_",fixed=TRUE)
  filesplit2 = filesplit2[[length(filesplit2)]]
  
  varin = paste(filesplit2[1],filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  varout2 = paste(varout,filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
  
  nctest = nc_open(filesin[i])
  
  if(i==1){
    dataunits = nctest$var[[1]]$units
    lon = ncvar_get(nctest,"lon")
    lat = ncvar_get(nctest,"lat")
    time=ncvar_get(nctest,"time")
    timeunits = nctest$var[[1]]$dim[[3]]$units
    latunits = nctest$var[[1]]$dim[[2]]$units
    lonunits = nctest$var[[1]]$dim[[1]]$units
    
    timeunits = paste(substr(timeunits,1,22),"12:00:00",sep="")
    
    indexdate = as.Date(substr(timeunits,12,21))
    dates=indexdate+time
    timesin = which(as.numeric(substr(dates,1,4))>=1976)
    dates = dates[timesin]
    timestart = timesin[1]
    timecount = length(timesin)
    timeout= time[timesin]
    
    if(tempres=="monthly"){
      yearmon = unique(substr(dates,1,7))
      timeout = timeout[which(substr(dates,9,10)=="01")]
    } 
    if(tempres=="annual"){
      years = unique(as.numeric(substr(dates,1,4)))
    } 
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    idxuse = which(modelgrid$lat>=latbnds[1] & modelgrid$lat<=latbnds[2] & modelgrid$lon>=lonbnds[1] & modelgrid$lon<=lonbnds[2])
    modelgriduse = modelgrid[idxuse,]
    
    lonstart = min(modelgriduse$R)
    latstart = min(modelgriduse$C)
    loncount = length(unique(modelgriduse$R))
    latcount = length(unique(modelgriduse$C))
    
    lonout = unique(modelgriduse$lon)
    latout = unique(modelgriduse$lat)
    
  }
  
  testdat = ncvar_get(nctest,varin,start=c(lonstart,latstart,timestart),count=c(loncount,latcount,timecount))
  if(var=="tasmax" | var=="tasmin"){
    if(testdat[140,140,1]>200){
      testdat = testdat-273.15
    }
  }
  
  if(mean(testdat,na.rm=TRUE)<1 & var=="pr"){
    testdat = testdat*86400
  }
  
  #if(var=="pr"){
  #  testdat=testdat*86400
  #}
  nc_close(nctest)
  
  if(varout=="FLI" | varout=="FBI"){
    file2 = paste("/data4/data/DS_proj/LOCA/tasmin/historical/tasmin",filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
    varin2 = paste("tasmin",filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
    nctest = nc_open(file2)
    testdat2 = ncvar_get(nctest,varin2,start=c(lonstart,latstart,timestart),count=c(loncount,latcount,timecount))
    nc_close(nctest)
    if(testdat2[140,140,1]>200){
      testdat2 = testdat2-273.15
    }
    testdat = testdat*9/5+32 # C to F conversion for springpheno
    testdat2 = testdat2*9/5+32
  }
  
  if(varout=="heatwaves"){
    file2 = paste("/data4/data/DS_proj/LOCA/tasmin/historical/tasmin",filesplit2[2],filesplit2[3],filesplit2[4],sep="_")
    varin2 = paste("tasmin",filesplit2[2],filesplit2[3],substr(filesplit2[4],1,nchar(filesplit2[4])-3),sep="_")
    nctest = nc_open(file2)
    testdat2 = ncvar_get(nctest,varin2,start=c(lonstart,latstart,timestart),count=c(loncount,latcount,timecount))
    nc_close(nctest)
    if(testdat2[140,140,1]>200){
      testdat2 = testdat2-273.15
    }
    
    tmaxq95file = paste("/home/woot0002/TJZenzal/tmaxq95",filesplit2[2],filesplit2[3],"historical_climo_SE.nc",sep="_")
    tminq95file = paste("/home/woot0002/TJZenzal/tminq95",filesplit2[2],filesplit2[3],"historical_climo_SE.nc",sep="_")
    
    nctest = nc_open(tmaxq95file)
    tmaxq95 = ncvar_get(nctest,"tmaxq95")
    nc_close(nctest)
    
    nctest = nc_open(tminq95file)
    tminq95 = ncvar_get(nctest,"tminq95")
    nc_close(nctest)
    
  }
  
  #testsfc = list(x=lon[unique(modelgriduse$R)]-360,y=lat[unique(modelgriduse$C)],z=testdat[,,1])
  #surface(testsfc,type="I")
  #map("state",add=TRUE)
  
  if(tempres=="climo"){
    if(varout=="pr"){
      tmp1 = apply(testdat,c(1,2),sum,na.rm=TRUE)  
    } 
    if(varout=="rx1day"){
      tmp1 = apply(testdat,c(1,2),max,na.rm=TRUE)
    }
    if(varout=="rx5day"){
      tmp1 = apply(testdat,c(1,2),calcrollsum,size=5)
    }
    if(varout=="tasmax" | varout=="tasmin"){
      tmp1 = apply(testdat,c(1,2),mean,na.rm=TRUE)
    }
    if(varout=="tmax100"){
      td = ifelse(testdat>=37.7778,1,0)
      tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
    }
    if(varout=="tmaxq95" | varout=="tminq95"){
      tmp1 = apply(testdat,c(1,2),quantile,probs=0.95,na.rm=TRUE)
    }
    
    tmp=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
    testdat=tmp
  }
  
  if(tempres=="annual"){
    tmp = array(NA,dim=c(loncount,latcount,length(years)))
    for(y in 1:length(years)){
      idxmon = which(substr(dates,1,4)==years[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="frd"){
        tmp1 = array(NA,dim=c(loncount,latcount))
        for(r in 1:lonout){
          for(c in 1:latout){
            message("Working on r ",r," and c ",c)
            tmp1[r,c] = lastfreeze(testdat[r,c,idxmon],startdate=startdate,inputtimes=timeout[idxmon],threshold=0)
          }
        }
      }
      if(varout=="FLI" | varout=="FBI"){
        tmp1 = array(NA,dim=c(loncount,latcount))
        for(r in 1:length(lonout)){
          for(c in 1:length(latout)){
            
            TMAX = testdat[r,c,idxmon]
            TMIN = testdat2[r,c,idxmon]
            if(length(idxmon)<366){
              TMAX = matrix(c(TMAX[1:59],NA,TMAX[60:365]),nrow=366,ncol=1)
              TMIN = matrix(c(TMIN[1:59],NA,TMIN[60:365]),nrow=366,ncol=1)
            } else {
              TMAX = matrix(TMAX,nrow=366,ncol=1)
              TMIN = matrix(TMIN,nrow=366,ncol=1)
            }
            
            if(all(is.na(TMAX)==TRUE)==FALSE & all(is.na(TMIN)==TRUE)==FALSE){
              
              calcsires = calc_si(TMAX=TMAX,TMIN=TMIN,lat=lat[c])
              if(varout=="FLI"){
                tmp1[r,c] = calcsires$FLImat[1,1]  
              }
              if(varout=="FBI"){
                tmp1[r,c] = calcsires$FBImat[1,1]   
              }
            }
            #spell_length_calc(inputdata,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="count")
            #message("Calcs complete for r ",r," and c ",c)
          }
        }
      }
      
      if(varout=="heatwaves"){
        tmp1 = array(NA,dim=c(loncount,latcount))
        for(r in 1:length(lonout)){
          for(c in 1:length(latout)){
            
            #tmaxmask = ifelse(testdat[r,c,idxmon]>=tmaxq95[r,c],1,0)
            #tminmask = ifelse(testdat2[r,c,idxmon]>=tminq95[r,c],1,0)
            #inputdata = ifelse(tmaxmask==1 & tminmask==1,1,0)
            
            tmp1[r,c] = heatwave.calc(tmaxdata=testdat[r,c,idxmon],tmindata=testdat2[r,c,idxmon],tmaxq95=tmaxq95[r,c],tminq95=tminq95[r,c])
            #spell_length_calc(inputdata,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="count")
            #message("Calcs complete for r ",r," and c ",c)
          }
        }
      }
      if(varout=="tmax100"){
        td = ifelse(testdat[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(years))
    }
    testdat=tmp
  }
  
  if(tempres=="monthly"){
    tmp = array(NA,dim=c(loncount,latcount,length(yearmon)))
    for(y in 1:length(yearmon)){
      idxmon = which(substr(dates,1,7)==yearmon[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="tmax100"){
        td = ifelse(testdat[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(yearmon))
    }
    testdat=tmp
  }
  
  if(tempres=="seasonal"){
    years = as.numeric(substr(dates,1,4))
    seas = rep(NA,length(dates))
    seas[which(as.numeric(substr(dates,6,7))>=12 | as.numeric(substr(dates,6,7))<=2)] = 1
    seas[which(as.numeric(substr(dates,6,7))>=3 & as.numeric(substr(dates,6,7))<=5)] = 2
    seas[which(as.numeric(substr(dates,6,7))>=6 & as.numeric(substr(dates,6,7))<=8)] = 3
    seas[which(as.numeric(substr(dates,6,7))>=9 & as.numeric(substr(dates,6,7))<=11)] = 4
    
    yearseasall = paste(years,seas,sep="-")
    yearseas = unique(yearseasall)
    
    tmp = array(NA,dim=c(loncount,latcount,length(yearseas)))
    for(y in 1:length(yearseas)){
      idxmon = which(yearseasall==yearseas[y])
      
      if(varout=="pr"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),sum,na.rm=TRUE)  
      } 
      if(varout=="rx1day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),max,na.rm=TRUE)
      }
      if(varout=="mdrn"){
        tp = ifelse(testdat[,,idxmon]>=0.254,1,0)
        tmp1 = apply(tp,c(1,2),sum,na.rm=TRUE)
      }
      if(varout=="rx5day"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),calcrollsum,size=5)
      }
      if(varout=="tasmax" | varout=="tasmin"){
        tmp1 = apply(testdat[,,idxmon],c(1,2),mean,na.rm=TRUE)
      }
      if(varout=="tmax100"){
        td = ifelse(testdat[,,idxmon]>=37.7778,1,0)
        tmp1 = apply(td,c(1,2),sum,na.rm=TRUE)
      }
      
      tmp[,,y]=ifelse(is.na(testdat[,,1])==TRUE,NA,tmp1)
      message("Finished with yearmon ",y," / ",length(yearseas))
    }
    testdat=tmp
  }
  
  #####
  if(tempres=="daily"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,".nc",sep="")
  }
  if(tempres=="monthly"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_mon_SE.nc",sep="")
  }
  if(tempres=="annual"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_ann_SE.nc",sep="")
    timeout = timeout[which(substr(dates,6,10)=="01-01")]
  }
  if(tempres=="climo"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_climo_SE.nc",sep="")
    timeout = timeout[1]
  }
  if(tempres=="seasonal"){
    filenameout = paste("/home/woot0002/TJZenzal/",varout2,"_seas_SE.nc",sep="")
    seasdiff = diff(seas)
    
    if(i==1){
      whichlist = c()
      yearsu = unique(years)
      monthsin = c(1,3,6,9)
      I=1
      for(j in 1:length(yearsu)){
        for(m in 1:length(monthsin)){
          whichlist[I] = which(as.numeric(substr(dates,1,4))==yearsu[j] & as.numeric(substr(dates,6,7))==monthsin[m] & as.numeric(substr(dates,9,10))==1)
          I=I+1
        }
      }
    }
    timeout = timeout[whichlist]
  }
  
  dimX <- ncdim_def( "lon", lonunits, lonout)
  dimY <- ncdim_def( "lat", latunits, latout)
  dimT <- ncdim_def("time",timeunits,timeout)
  
  # Make varables of various dimensionality, for illustration purposes
  mv <- 1E20 # missing value to use
  
  var1d <- ncvar_def(varout,dataunits,longname=varout, list(dimX,dimY,dimT), mv ,compression=9)
  
  #######
  # Create netcdf file
  
  nc <- nc_create(filenameout ,  var1d )
  # Write some data to the file
  ncvar_put(nc, var1d, testdat) # no start or count: write all values\
  # close ncdf
  nc_close(nc)
  rm(testdat)
  gc()
  ptmend = proc.time()-ptm
  
}

