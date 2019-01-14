########
# Step 8 analysis 3^5
# pull time series for location

ts_pull = function(varname,histfilelist,projfilelist,appfunc,loc_name,loc_lon,loc_lat,TC,TH,cond,seasonin,useobs=FALSE){

# must pass in the following: histfilelist,histfilebreakdown, projfilelist,projfilebreakdown,TC,TH,cond, appfunc, varname, futureperiod,tempperiod,varunits,changeunits,usecompound

# requirements
# varname - short varname, like tasmax, tasmin, pr, tmax95,tmin32, pr25,mdrn
# histlist - comma separated list of model historical files
# projlist - comma separated list of model future files
# appfunc - function to apply, like mean, sum, or any of those in the analysisfunctions.R file
# difftype - type of differencing function, "absolute" or "percent"
# varunits - units of the variables themselves, used for netcdf writing
# changeunits - units of the projected change, used for netcdf writing
# futureperiod - comma delimited start and end year, 1981,2005
# TC - is this a threshold calculation? TRUE or FALSE (default)
# TH - threshold value, only used if TC is TRUE. Default is NA
# cond - condition connected with the threshold. Used if TC is TRUE. gt,gte,lt,lte are options
  
  histfilelist = do.call("c",strsplit(histfilelist,",",fixed=TRUE))
  projfilelist = do.call("c",strsplit(projfilelist,",",fixed=TRUE))
 
  varin = varname
  if(varname=="tmax95" | varname=="tmax100") varin="tasmax" # for the tmax95 and other threshold variables, you need to use the base variable and calculate the number of days matching the threshold.
  if(varname=="tmin32" | varname=="tmin28" | varname=="frd") varin="tasmin"
  if(varname=="pr25" | varname=="pr50" | varname=="mdrn" | varname=="rx1day" | varname=="rx5day" | varname=="cdd" | varname=="cwd") varin="pr"
  
  if(varname=="heatwaves" | varname=="gsl"){
    if(length(grep("tasmax",projfilelist))==length(projfilelist)) varin="tasmax"
    if(length(grep("tasmin",projfilelist))==length(projfilelist)) varin="tasmin"
  }
  
  if(seasonin=="daily" & (varname=="tmax95" | varname=="tmax100" | varname=="tmin28" | varname=="tmin32" | varname=="frd" | varname=="pr50" | varname=="pr25" | varname=="mdrn" | varname=="cwd" | varname=="cdd" | varname=="rx1day" | varname=="rx5day" | varname=="gsl" | varname=="heatwaves")){
    stop("Daily timeseries can only be pulled for tasmax, tasmin, or pr",.call=TRUE)
  }
  
   
  # Creating file breakdown tables
filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9)
filebreakdown3$obs = rep(c("Daymet","Livneh","PRISM"),9)
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

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

#############
# 1a- historical calcs

dates = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")

for(i in 1:length(histfilelist)){
  ptm = proc.time()
  message("Starting work on file ",histfilelist[i])
  if(histfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  } else {
    noleap=TRUE
  }
  #### dates to use
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(appfunc!="heatwaves" & appfunc!="growing_season_length"){
  if(seasonin=="ann"){
    yearlyoutput = netcdftoyearlycalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)
  } else {
    if(seasonin=="DJF" | seasonin =="MAM" | seasonin =="JJA" | seasonin=="SON"){
      yearlyoutput = netcdftoseasonalcalcs(histfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(1981,2005),datesdataperiod=datesin,appliedfunction=appfunc)  
    }
    if(seasonin=="daily" & (varname=="tasmax" | varname=="tasmin" | varname=="pr")){
      yearlyoutput = netcdftodailypointseries(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),datesdataperiod=datesin,loc_lat=loc_lat,loc_lon=loc_lon)
    }
  }
  }
  
  if(appfunc=="growing_season_length" | appfunc=="heatwaves"){
    message("Running compound calculation, please be patient")
    yearlyoutput = netcdfcombocalcs(histfilelist[i],varname=varin,dimnames=c("lon","lat","time"),yearlydataperiod=c(1981,2005),datesdataperiod=datesin,combofunction=appfunc)
  }
  
  #print(range(yearlyoutput[[3]],na.rm=TRUE))
  
  if(i==1){
    lon = yearlyoutput[[1]]
    lat = yearlyoutput[[2]]
    
    ###
    # create model grid
    
    LON = rep(lon,each=length(lat))
    LAT = rep(lat,length(lon))
    R = rep(1:length(lon),each=length(lat))
    C = rep(1:length(lat),length(lon))
    modelgrid = data.frame(R,C,LON,LAT)
    names(modelgrid) = c("R","C","lon","lat")
    if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
    
    ###
    # get cells to use
    loc_lon = as.numeric(loc_lon)
    loc_lat = as.numeric(loc_lat)
    if(loc_lon>0) loc_lon=loc_lon-360
    
    pointarea = distfunc(loc_lon,loc_lat,modelgrid)
    locstart = c(pointarea$R-1,pointarea$C-1)
    locend = c(pointarea$R+1,pointarea$C+1)
    
    ###
    # get point values
    if(seasonin!="daily"){
    histlist = matrix(NA,nrow=length(unique(substr(datesin,1,4))),ncol=length(histfilelist))
    } else {
      histlist = matrix(NA,nrow=length(datesin),ncol=length(histfilelist))
    }
  }
  
  if(seasonin!="daily"){
    message("creating histlist")
    histlist[,i] = apply(yearlyoutput[[3]][locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
  } else {
    
    if(length(yearlyoutput[[3]])>nrow(histlist)){
      yearlyoutput[[3]] = yearlyoutput[[3]][-which(substr(dates,6,10)=="02-29")]
    }
    
    histlist[,i] = yearlyoutput[[3]]
  }
  
  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(histfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
#print(range(histlist,na.rm=TRUE))

#############
# 1b- Future calcs

dates = seq(as.Date("2006-01-01"),as.Date("2099-12-31"),by="day")

for(i in 1:length(projfilelist)){
  ptm = proc.time()
  message("Starting work on file ",projfilelist[i])
  if(projfilebreakdown$GCM[i]=="MPI-ESM-LR"){
    noleap=FALSE
  }else {
    noleap=TRUE
  }
  if(noleap==FALSE) datesin = dates
  if(noleap==TRUE) datesin = dates[-which(substr(dates,6,10)=="02-29")]
  
  if(appfunc!="heatwaves" & appfunc!="growing_season_length"){
  if(seasonin=="ann"){
    yearlyoutput = netcdftoyearlycalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(2006,2099),datesdataperiod=datesin,appliedfunction=appfunc)
  } else {
    if(seasonin=="DJF" | seasonin =="MAM" | seasonin =="JJA" | seasonin=="SON"){
      yearlyoutput = netcdftoseasonalcalcs(projfilelist[i],varname=varin,season=seasonin,dimnames=c("lon","lat","time"),threscalc=TC,thres=TH,condition=cond,yearlydataperiod=c(2006,2099),datesdataperiod=datesin,appliedfunction=appfunc)
    }
    if(seasonin=="daily" & (varname=="tasmax" | varname=="tasmin" | varname=="pr")){
      yearlyoutput = netcdftodailypointseries(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),datesdataperiod=datesin,loc_lat=loc_lat,loc_lon=loc_lon)
    }
    }
  }
  
  if(appfunc=="growing_season_length" | appfunc=="heatwaves"){
    yearlyoutput = netcdfcombocalcs(projfilelist[i],varname=varin,dimnames=c("lon","lat","time"),yearlydataperiod=c(2006,2099),datesdataperiod=datesin,combofunction=appfunc)
  }
  
  #print(range(yearlyoutput[[3]],na.rm=TRUE))
  
  if(i==1){
    if(seasonin!="daily"){
      projlist = matrix(NA,nrow=length(unique(substr(datesin,1,4))),ncol=length(projfilelist))
    } else {
      projlist = matrix(NA,nrow=length(datesin),ncol=length(projfilelist))
    }
  }
  
  if(seasonin!="daily"){
    projlist[,i] = apply(yearlyoutput[[3]][locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
  } else {
    if(length(yearlyoutput[[3]])>nrow(projlist)){
      yearlyoutput[[3]] = yearlyoutput[[3]][-which(substr(dates,6,10)=="02-29")]
    }
    projlist[,i] = yearlyoutput[[3]]
  }

  rm(yearlyoutput)
  gc()
  ptmend = proc.time()
  message("Finished with file ",i," / ",length(projfilelist))
  message("Time to complete gathering: ",ptmend[3]-ptm[3]," secs")
}
#print(range(projlist,na.rm=TRUE))

################
# 1-c time series compilation

tseriesout = matrix(NA,nrow=nrow(histlist)+nrow(projlist),ncol=ncol(projlist))
for(i in 1:length(projfilelist)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  DSin = projfilebreakdown$DS[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin & histfilebreakdown$DS==DSin)
  tseriesout[1:nrow(histlist),i]=histlist[,histidx]
  tseriesout[(nrow(histlist)+1):nrow(tseriesout),i] = projlist[,i]
}

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
if(seasonin!="daily"){
year = 1981:2099
tseriesout = data.frame(year,tseriesout)
} else {
  date = seq(as.Date("1981-01-01"),as.Date("2099-12-31"),by="day")
  date = date[-which(substr(date,6,10)=="02-29")]
  tseriesout = data.frame(year,tseriesout)
}
names(tseriesout)[2:ncol(tseriesout)] = projnotes

#######################
# 1-d Include observed values - optional

if(includeobs==TRUE & seasonin=="daily"){
  years= 1981:2005
  if(varin=="tasmax") varinobs = "tmmx"
  if(varin=="tasmin") varinobs = "tmmn"
  if(varin=="pr") varinobs = "pr"
  for(y in 1:length(years)){
    
    test = nc_open(paste("/data4/data/DS_proj/METDATA/",varinobs,"_",years[y],".nc",sep=""))
    
    if(y==1){
      lonobs = ncvar_get(test,"lon")
      latobs = ncvar_get(test,"lat")
      latobs = rev(latobs)
      
      LON = rep(lonobs,each=length(latobs))
      LAT = rep(latobs,length(lonobs))
      C = rep(1:length(latobs),each=length(lonobs))
      R = rep(1:length(lonobs),length(latobs))
      modelgrid = data.frame(R,C,LON,LAT)
      names(modelgrid) = c("R","C","lon","lat")
      if(all(modelgrid$lon>0)==TRUE) modelgrid$lon = modelgrid$lon-360
      
      ###
      # get cells to use
      loc_lon = as.numeric(loc_lon)
      loc_lat = as.numeric(loc_lat)
      if(loc_lon>0) loc_lon=loc_lon-360
      
      pointarea = distfunc(loc_lon,loc_lat,modelgrid)
      locstart = c(pointarea$R-2,pointarea$C-2)
      locend = c(pointarea$R+2,pointarea$C+2)
      
      if(seasonin!="daily"){
        yeararray = c()
      } else {
        datesobs = seq(as.Date("1981-01-01"),as.Date("2005-12-31"),by="day")
        yeararray = c()
      }
    }
    
    vardata = ncvar_get(test,test$var[[1]]$name)
      vardata2 = array(NA,dim=c(length(lonobs),length(latobs),dim(vardata)[3]))
      for(d in 1:dim(vardata)[3]){
        tmp = t(vardata[,,d])
        vardata2[,,d]=tmp[,length(latobs):1]
      }
    
    time = ncvar_get(test,"day")
    timeunits = test$var[[1]]$dim[[3]]$units
    nc_close(test)
    
    startdate = as.Date(substr(timeunits,12,21))
    times = startdate+time
    
    if(seasonin!="daily"){
    if(varname == "tasmax" | varname == "tasmin"){
      temp2 = mean(vardata2[locstart[1]:locend[1],locstart[2],locend[2]],na.rm=TRUE)
    }
    
    if(varname=="gsl"){ 
      test = nc_open(paste("/data4/data/DS_proj/METDATA/tmmn_",years[y],".nc",sep=""))
      vardata2 = ncvar_get(test,test$var[[1]]$name)
      nc_close(test)
      temp2 = matrix(NA,nrow=dim(vardata2)[1],ncol=dim(vardata2)[2])
      for(r in 1:length(lonobs)){
        for(c in 1:length(latobs)){
          temp2[r,c] = GSLcalc(vardata2[r,c,],vardata2[r,c,],inputtimes=time,startdate=startdate)
          #message("Finished Calculation R ",r," and C ",c)
        }
      }  
      temp2 = mean(vardata2[locstart[1]:locend[1],locstart[2],locend[2]],na.rm=TRUE)
    }
    
    
    if(varname == "pr"){
      temp2 = apply(vardata2[locstart[1]:locend[1],locstart[2],locend[2]],c(1,2),sum,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2)
    }
    
    if(varname == "tmax95" | varname == "tmax100" | varname == "tmin32" | varname == "tmin28" | varname == "pr25" | varname == "pr50" | varname == "mdrn"){
      if(varname == "tmax100")  temp = ifelse(vardata2>=310.928,1,0)
      if(varname == "tmax95")  temp = ifelse(vardata2>=308.15,1,0)
      if(varname == "tmin32")  temp = ifelse(vardata2<=273.15,1,0)
      if(varname == "tmin28")  temp = ifelse(vardata2<=270.928,1,0)
      if(varname == "pr25")  temp = ifelse(vardata2>=50.8,1,0)
      if(varname == "pr50")  temp = ifelse(vardata2>=25.4,1,0)
      if(varname == "mdrn")  temp = ifelse(vardata2>=0.254,1,0)
      temp2 = apply(temp[locstart[1]:locend[1],locstart[2],locend[2]],c(1,2),sum,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata2[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2)
      rm(temp)
    }
    
    if(varname=="rx1day"){
      temp2 = apply(vardata2,c(1,2),max,na.rm=TRUE)
      temp2 = ifelse(is.na(vardata2[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2[locstart[1]:locend[1],locstart[2],locend[2]])
    }
    if(varname=="rx5day"){
      temp2 = apply(vardata2,c(1,2),calcrollsum,size=5)
      temp2 = ifelse(is.na(vardata[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2[locstart[1]:locend[1],locstart[2],locend[2]])
    }
    if(varname=="frd"){
      #test1 = apply(tempdata,c(1,2),lastfreeze,startdate=startdate,inputtimes=times)
      ptm = proc.time()
      temp2 = matrix(NA,nrow=length(lonobs),ncol=length(latobs))
      for(r in 1:length(lonobs)){
        for(c in 1:length(latobs)){
          message("Working on r ",r," and c ",c)
          temp2[r,c] = lastfreeze(vardata2[r,c,],startdate=startdate,inputtimes=times)
        }
      }
      temp2 = ifelse(is.na(vardata2[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2[locstart[1]:locend[1],locstart[2],locend[2]])
      ptm.end = proc.time()
    }
      
    if(varname=="cdd"){
      temp2 = apply(vardata2,c(1,2),spell_length_calc,premasked=FALSE,cond="LT",spell_len=3,thold=0.254,outtype="max")
      temp2 = ifelse(is.na(vardata2[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2[locstart[1]:locend[1],locstart[2],locend[2]])
    }
    if(varname=="cwd"){
      temp2 = apply(vardata2,c(1,2),spell_length_calc,premasked=FALSE,cond="GE",spell_len=3,thold=0.254,outtype="max")
      temp2 = ifelse(is.na(vardata2[,,1])==TRUE,NA,temp2)
      temp2 = mean(temp2[locstart[1]:locend[1],locstart[2],locend[2]])
    }
    yeararray[y] = temp2
    } else{
      dateidx = which(as.numeric(substr(datesobs,1,4))==years[y])
      yeararray[dateidx]=apply(vardata2[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
    }
    rm(vardata)
    rm(vardata2)
    rm(temp2)
    message("Finished data grab for year: ",years[y])
  }

  if(seasonin=="daily"){
    leapdays = which(substr(datesobs,6,10))
    yeararray = yeararray[-leapdays]
  }
  
  tseriesout$METDATA = NA
  if(seasonin=="daily"){
    datesall= seq(as.Date("1981-01-01"),as.Date("2099-12-31"),by="day")
    datesall = datesall[-which(substr(datesall,6,10)=="02-29")]
    inputidx = which(as.numeric(substr(datesall,1,4))>=1981 & as.numeric(substr(datesall,1,4))<=2005)
    tseriesout$METDATA[inputidx]=yeararray
  } else {
    yearsall = 1981:2099
    inputidx = which(yearsall>=1981 & yearsall<=2005)
    tseriesout$METDATA[inputidx]=yeararray
  }
}

# notes on the order of the ensemble members for the output netcdf file.

#######################
# 1-e Write out csv time series

step8_filename = paste(varname,"_allmem_timeseries_",loc_name,"_",seasonin,".csv",sep="")

write.table(tseriesout,paste("/data2/3to5/I35/timeseries/",step8_filename,sep=""),sep=",",row.names=FALSE)

#######
# return output

returnlist = list(step8_filename,projnotes,histnotes)
return(returnlist)
}

###
# Argument parser

library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

ParseArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c('-v', "--var"), action="store", default='na',
                dest='varname',
                help=paste("The variable of interest (tasmax, tasmin, pr,tmax95,tmin32,pr25) in the input file.", 
                           "If not present, an error will be thrown.")),
    make_option(c("-i", "--histlist"), action="store", default="",
                dest='histlist',
                help=paste("Comma delimited list of filenames for the model historical input files;", 
                           "an error will be thrown if this argument is not present")),
    make_option(c("-p", "--projlist"), action="store", default="",
                dest='projlist',
                help=paste("Comma delimited list of filenames for the model future input files;", 
                           "an error will be thrown if this argument is not present")),
    make_option(c("-a","--appfunc"), action="store", default='mean',dest='appfunc',
                help=paste("The function to be applied prior to determining climatology. ", 
                           "Includes the following options: mean (default), sum,",
                           " max, rx5day (5 day maximum), lastfreeze (Date of last freeze),",
                           "maxdryspell (longest dry spell), and maxwetspell (longest wet spell). If seasonin is daily this will be ignored.")),
    make_option(c("-T","--thres_calc"), action="store", default=FALSE,dest='TC',
                help=paste("Option to calculate certain variables that are threshold based. ", 
                           "Defaults to FALSE, but must be TRUE to calculate variables like ",
                           "tmax95, tmin32, mdrn, pr25. ",
                           "Should also be paired with appfunc option sum. ")),
    make_option(c("-H","--thres_val"), action="store", default='na',dest='TH',
                help=paste("The value of the threshold to be applied if -T is TRUE. ", 
                           "Defaults to a missing value otherwise. Will throw an error ",
                           "if TC is TRUE and -H is not specified")),
    make_option(c("-c","--condition"), action="store", default='na',dest='cond',
                help=paste("The condition that determines how the threshold will", 
                           "be applied to the input data. May be one of LT, LE,", 
                           "GT, GE")),
    make_option(c("-S","--season"), action="store", default='ann',dest='seasonin',
                help=paste("The season to use for calculations. Defaults to 'ann', ", 
                           "but can also be DJF, MAM, JJA, SON, or daily for a daily timeseries. ",
                           "If seasonin is daily, then the only variables that varnames that can be used are tasmax, tasmin, and pr. Will throw an error otherwise.")),
    make_option(c("-n", "--loc_name"), action="store",
                dest='loc_name',
                help=paste("Name of the Location you are grabbing information for. ", 
                           "Anything can be used, but it's recommended to use something short and simple. ",
                           "For example, using OKC instead of Oklahoma City. If not provided, the script will throw an error. ")),
    make_option(c("-x", "--loc_lon"), action="store",
                dest='loc_lon',
                help=paste("Longitude of the desired location. ", 
                           "If not provided, the script will throw an error. ")),
    make_option(c("-y", "--loc_lat"), action="store",
                dest='loc_lat',
                help=paste("Latitude of the desired location. If not provided, the script will throw an error. ")),
    make_option(c("-O", "--useobs"), action="store",default=FALSE,
                dest='useobs',
                help=paste("Include the historical observations from the METDATA dataset for the same location. ",
                           "Given METDATA has a finer resolution, a 5 by 5 block around the location is taken compared to a 3x3 block in 3^5."))
  )
  
  description = paste('Given the desired variable name, the list of files (with paths) ', 
                      "and associated options for thresholds, future period, applied function, temporal resolution, and location,  ",
                      "pull a time series for a location of interest ",
                      "from all available members from a collection of downscaled projections.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog --var varname --histlist histlist --projlist projlist ", 
                "[-a appfunc] [-T calcthres] [-H thresvalue] [-c condition] [-S seasonin]", 
                "-n loc_name -x loc_lon -y loc_lat [-O useobs]", 
                " [-h --help]")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

ts_pull(varname=parsed.args$varname,histfilelist=parsed.args$histlist,projfilelist=parsed.args$projlist,
      appfunc=parsed.args$appfunc,TC=parsed.args$TC,TH=parsed.args$TH,cond=parsed.args$cond,
      seasonin=parsed.args$seasonin,useobs=parsed.args$useobs,loc_name=parsed.args$loc_name,loc_lon=parsed.args$loc_lon,loc_lat=parsed.args$loc_lat)

