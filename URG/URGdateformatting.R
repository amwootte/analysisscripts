######################
#
# Attach dates to files

library(graphics)
setwd("/home/woot0002/Blanco/")
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(ncdf4)
library(mailR)

#filestoprocess = system("ls /data2/3to5/I35/Navajo/*/*/*/*.txt",intern=TRUE)
filestoprocess = system("ls /data2/3to5/I35/Blanco/tasmin/*.txt",intern=TRUE)

for(i in 1:length(filestoprocess)){
  
  filedat = read.table(filestoprocess[i],sep=" ",header=TRUE)
  filesplit = strsplit(filestoprocess[i],"/")
  filesplit2 = strsplit(filesplit[[1]][7],"_")[[1]]
  varname = filesplit2[1]
  period = c(as.numeric(substr(filesplit2[7],1,4)),as.numeric(substr(filesplit2[7],6,9)))
  
  dates = seq(as.Date(paste(period[1],"-01-01",sep="")),as.Date(paste(period[2],"-12-31",sep="")),by="day")
  if(length(dates)>nrow(filedat)){
    datesin = dates[-which(substr(dates,6,10)=="02-29")]
  } else {
    datesin = dates
  }
  
  year = as.numeric(substr(datesin,1,4))
  month = as.numeric(substr(datesin,6,7))
  day = as.numeric(substr(datesin,9,10))
  hr = min = sec = rep(0,length(datesin))
  
  dateframe = data.frame(year,month,day,hr,min,sec)
  
  filedat = cbind(dateframe,filedat)
  fileout=  filesplit[[1]][7]
  message("Writing out file")
  ptm = proc.time()
  write.table(filedat,file=fileout,col.names=FALSE,row.names=FALSE,sep=" ")
  ptmend = proc.time()-ptm
  message("Finished writing file after ",ptmend[3]," secs")
  message("Finished data formatting for file ",i," / ",length(filestoprocess))
}

send.mail(from = "amwootte@ou.edu",
          to = "amwootte@ou.edu",
          subject = "message from R on climatedata",
          body = "URGdateformatting.R has finished running", 
          authenticate = TRUE,
          smtp = list(host.name = "smtp.office365.com", port = 587,
                      user.name = "amwootte@ou.edu", passwd = "D0wnSc2l!ng", tls = TRUE))



