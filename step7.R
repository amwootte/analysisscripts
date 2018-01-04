##################
#
# Step 7: Location grab and calculation

source("analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

# must pass in the following: loc_lon,loc_lat,loc_name,step1_filename, varname, futureperiod,tempperiod,varunits,changeunits,histnotes,projnotes

test = nc_open(step1_filename)
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
histlist = ncvar_get(test,"histmean")
projlist = ncvar_get(test,"projmean")
diffs = ncvar_get(test,"projmeandiff")
nc_close(test)

histfilebreakdown = do.call(rbind,strsplit(histnotes,",",fixed=TRUE))
names(histfilebreakdown) = c("GCM","DS","obs")

projfilebreakdown = do.call(rbind,strsplit(projnotes,",",fixed=TRUE))
names(projfilebreakdown) = c("GCM","DS","obs","scen")

###
# create model grid

LON = rep(lon,each=length(lat))
LAT = rep(lat,length(lon))
R = rep(1:length(lon),each=length(lat))
C = rep(1:length(lat),length(lon))
modelgrid = data.frame(R,C,LON,LAT)
names(modelgrid) = c("R","C","lon","lat")
modelgrid$lon = modelgrid$lon-360

###
# get cells to use

pointarea = distfunc(loc_lon,loc_lat,modelgrid)
locstart = c(pointarea$R-1,pointarea$C-1)
locend = c(pointarea$R+1,pointarea$C+1)

###
# get point values

histvals = apply(histlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
projvals = apply(projlist[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)
diffvals = apply(diffs[locstart[1]:locend[1],locstart[2]:locend[2],],3,mean,na.rm=TRUE)

####
# Create output table

projchangedat = projfilebreakdown
projchangedat$modfut = projvals
projchangedat$projchange = diffvals

projchangedat$modhist = NA

for(i in 1:nrow(projchangedat)){
  GCMin = projfilebreakdown$GCM[i]
  obsin = projfilebreakdown$obs[i]
  DSin = projfilebreakdown$DS[i]
  histidx = which(histfilebreakdown$GCM==GCMin & histfilebreakdown$obs==obsin & histfilebreakdown$DS==DSin)
  projchangedat$modhist[i] = histvals[histidx]
}

projchangedat = projchangedat[,c(3,5,7:11)]


####
# Write out file

filename = paste("/data2/3to5/I35/point_output/",varname,"_",futureperiod,"_",location_name,"_",difftype,".csv",sep="")
write.table(projchangedat,file=filename,row.names=FALSE,sep=",")

print("Location outfile ",filename," written!")







