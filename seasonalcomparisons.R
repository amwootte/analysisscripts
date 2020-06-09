########
# Step 3 analysis 3^5
# Plot individual members

#plot_indiv = function(step1_filename,projnotes,colorchoicediff,BINLIMIT,diffbartype,use_fixed_scale,fixed_scale){

###################################################################################################
## must pass in the following: step1_filename ,projnotes, colorchoicediff, BINLIMIT, diffbartype ##
###################################################################################################
filename_p0_ann = "/data2/3to5/I35/all_mems/pr_allmem_percent_2006-2035_ann.nc"
filename_p1_ann = "/data2/3to5/I35/all_mems/pr_allmem_percent_2036-2065_ann.nc"
filename_p2_ann = "/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_ann.nc"

filename_p0_DJF = "/data2/3to5/I35/all_mems/pr_allmem_percent_2006-2035_DJF.nc"
filename_p1_DJF = "/data2/3to5/I35/all_mems/pr_allmem_percent_2036-2065_DJF.nc"
filename_p2_DJF = "/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_DJF.nc"

filename_p0_MAM = "/data2/3to5/I35/all_mems/pr_allmem_percent_2006-2035_MAM.nc"
filename_p1_MAM = "/data2/3to5/I35/all_mems/pr_allmem_percent_2036-2065_MAM.nc"
filename_p2_MAM = "/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_MAM.nc"

filename_p0_JJA = "/data2/3to5/I35/all_mems/pr_allmem_percent_2006-2035_JJA.nc"
filename_p1_JJA = "/data2/3to5/I35/all_mems/pr_allmem_percent_2036-2065_JJA.nc"
filename_p2_JJA = "/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_JJA.nc"

filename_p0_SON = "/data2/3to5/I35/all_mems/pr_allmem_percent_2006-2035_SON.nc"
filename_p1_SON = "/data2/3to5/I35/all_mems/pr_allmem_percent_2036-2065_SON.nc"
filename_p2_SON = "/data2/3to5/I35/all_mems/pr_allmem_percent_2070-2099_SON.nc"

colorchoicehist = "whitetogreen"
colorchoicediff = "browntogreen"

BINLIMIT = 30
diffbartype = "difference"
use_fixed_scale = TRUE
fixed_scale=c(-70,110)
varin="pr"

library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

check.integer <- function(N){
  !grepl("[^[:digit:]]", format(N,  digits = 20, scientific = FALSE))
}

# google doc color bar options: https://docs.google.com/spreadsheets/d/1tZjLyKVxdvBW1z73a_8LEEdHBR_0Sp296oT4YeENExY/edit?usp=sharing
# colorbrewer options for color issues: http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=11
## tasmax, tasmin, tmax95
#colorchoicehist = "yellowtored"
#colorchoicediff = "bluetored"
# ## tmin32, Tasmin 

histfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*00_historical*.nc",sep=""),intern=T)
projfilelist = system(paste("ls /data2/3to5/I35/",varin,"/*/",varin,"_day_*rcp*.nc",sep=""),intern=T)

filebreakdown = do.call(rbind,strsplit(projfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
projfilebreakdown = filebreakdown3
rm(filebreakdown3)

filebreakdown = do.call(rbind,strsplit(histfilelist,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,3],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,4:7])
filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,4])))
filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,4])))
filebreakdown3 = filebreakdown3[,-c(3,8,9)]
names(filebreakdown3) = c("var","tempres","DS","code","scen","experiment","GCM","obs")
histfilebreakdown = filebreakdown3
rm(filebreakdown3)
rm(filebreakdown2)
rm(filebreakdown)

projnotes = paste(projfilebreakdown$GCM,projfilebreakdown$DS,projfilebreakdown$obs,projfilebreakdown$scen,sep="_")
histnotes = paste(histfilebreakdown$GCM,histfilebreakdown$DS,histfilebreakdown$obs,sep="_")

#projnotes = do.call("c",strsplit(projnotes,",",fixed=TRUE))

split1 = strsplit(filename_p1_ann,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod1 = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
seasonin = substr(split2[5],1,3)

split1 = strsplit(filename_p2_ann,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod2 = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
seasonin = substr(split2[5],1,3)

split1 = strsplit(filename_p0_ann,"/",fixed=TRUE)[[1]]
split2 = strsplit(split1[length(split1)],"_",fixed=TRUE)[[1]]

varname = split2[1]
futureperiod0 = c(as.numeric(substr(split2[4],1,4)),as.numeric(substr(split2[4],6,9)))
difftype = split2[3]
seasonin = substr(split2[5],1,3)

#####
# annual data
test = nc_open(filename_p1_ann) #data for period 1
lon = ncvar_get(test,"lon")
lat = ncvar_get(test,"lat")
histvals_p1_ann = ncvar_get(test,"histmean")
diffs_p1_ann = ncvar_get(test,"projmeandiff")
varunits = test$var[[1]]$units
changeunits = test$var[[3]]$units
nc_close(test)

test = nc_open(filename_p2_ann) # data for period 2
histvals_p2_ann = ncvar_get(test,"histmean")
diffs_p2_ann = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p0_ann) # data for period 2
histvals_p0_ann = ncvar_get(test,"histmean")
diffs_p0_ann = ncvar_get(test,"projmeandiff")
nc_close(test)

#####
# DJF data
test = nc_open(filename_p1_DJF) #data for period 1
histvals_p1_DJF = ncvar_get(test,"histmean")
diffs_p1_DJF = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p2_DJF) # data for period 2
histvals_p2_DJF = ncvar_get(test,"histmean")
diffs_p2_DJF = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p0_DJF) # data for period 2
histvals_p0_DJF = ncvar_get(test,"histmean")
diffs_p0_DJF = ncvar_get(test,"projmeandiff")
nc_close(test)

#####
# MAM data
test = nc_open(filename_p1_MAM) #data for period 1
histvals_p1_MAM = ncvar_get(test,"histmean")
diffs_p1_MAM = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p2_MAM) # data for period 2
histvals_p2_MAM = ncvar_get(test,"histmean")
diffs_p2_MAM = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p0_MAM) # data for period 2
histvals_p0_MAM = ncvar_get(test,"histmean")
diffs_p0_MAM = ncvar_get(test,"projmeandiff")
nc_close(test)

#####
# JJA data
test = nc_open(filename_p1_JJA) #data for period 1
histvals_p1_JJA = ncvar_get(test,"histmean")
diffs_p1_JJA = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p2_JJA) # data for period 2
histvals_p2_JJA = ncvar_get(test,"histmean")
diffs_p2_JJA = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p0_JJA) # data for period 2
histvals_p0_JJA = ncvar_get(test,"histmean")
diffs_p0_JJA = ncvar_get(test,"projmeandiff")
nc_close(test)

#####
# SON data
test = nc_open(filename_p1_SON) #data for period 1
histvals_p1_SON = ncvar_get(test,"histmean")
diffs_p1_SON = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p2_SON) # data for period 2
histvals_p2_SON = ncvar_get(test,"histmean")
diffs_p2_SON = ncvar_get(test,"projmeandiff")
nc_close(test)

test = nc_open(filename_p0_SON) # data for period 2
histvals_p0_SON = ncvar_get(test,"histmean")
diffs_p0_SON = ncvar_get(test,"projmeandiff")
nc_close(test)

#####

Daymetidx = which(projfilebreakdown$obs=="Daymet")

dvec = c(diffs_p0_SON[,,Daymetidx],diffs_p0_JJA[,,Daymetidx],diffs_p0_MAM[,,Daymetidx],diffs_p0_DJF[,,Daymetidx],diffs_p0_ann[,,Daymetidx],
         diffs_p1_SON[,,Daymetidx],diffs_p1_JJA[,,Daymetidx],diffs_p1_MAM[,,Daymetidx],diffs_p1_DJF[,,Daymetidx],diffs_p1_ann[,,Daymetidx],
         diffs_p2_SON[,,Daymetidx],diffs_p2_JJA[,,Daymetidx],diffs_p2_MAM[,,Daymetidx],diffs_p2_DJF[,,Daymetidx],diffs_p2_ann[,,Daymetidx])
diffcolorbar = colorramp(dvec,colorchoice=colorchoicediff,Blimit=BINLIMIT,type=diffbartype,use_fixed_scale = use_fixed_scale,fixed_scale=fixed_scale)

projfilebreakdown$DS = ifelse(projfilebreakdown$DS=="QDM","EDQM",as.character(projfilebreakdown$DS))

scens = unique(as.character(projfilebreakdown$scen))


######
# 2006-2035
pdf(paste("/data2/3to5/I35/plots/all_mems/",varname,"_2006-2035_Daymetonly.pdf",sep=""),onefile=TRUE,width=15,height=15)
par(mfrow=c(3,3))
for(i in 1:length(scens)){
    
    idxin = which(projfilebreakdown$obs=="Daymet" & as.character(projfilebreakdown$scen)==scens[i])
  
    for(s in 1:5){
      
      if(s==1){
        tmpin = diffs_p0_ann
        seasonin = "annual"
      }
      
      if(s==2){
        tmpin = diffs_p0_DJF
        seasonin = "DJF"
      }
      
      if(s==3){
        tmpin = diffs_p0_MAM
        seasonin = "MAM"
      }
      
      if(s==4){
        tmpin = diffs_p0_JJA
        seasonin = "JJA"
      }
      
      if(s==5){
        tmpin = diffs_p0_SON
        seasonin = "SON"
      }
      
      for(j in 1:length(idxin)){
        GCM = projfilebreakdown$GCM[idxin[j]]
        scen = projfilebreakdown$scen[idxin[j]]
        obs = projfilebreakdown$obs[idxin[j]]
        DS = projfilebreakdown$DS[idxin[j]]
        testsfc1 = list(x=lon-360,y=lat,z=tmpin[,,idxin[j]])
        surface(testsfc1,type="I",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=legadd)
        map("state",add=TRUE) # main=paste("Scen: ",scen," GCM: ",GCM,"\n DS: ", DS," Obs: ",obs,sep="")
        text(-108.85,28.45,labels=GCM,cex=2,pos = 4)
        text(-108.85,27.45,labels=DS,cex=2,pos = 4)
        text(-108.85,26.45,labels=obs,cex=2,pos = 4)
        text(-90.05,26.45,labels=scen,cex=2,pos = 2)
        text(-90.05,27.45,labels=paste("Season: ",seasonin,sep=""),cex=1.75,pos = 2)
      }
    }  
}
dev.off()


######
# 2036-2065
pdf(paste("/data2/3to5/I35/plots/all_mems/",varname,"_2036-2065_Daymetonly.pdf",sep=""),onefile=TRUE,width=15,height=15)
par(mfrow=c(3,3))
for(i in 1:length(scens)){
  
  idxin = which(projfilebreakdown$obs=="Daymet" & as.character(projfilebreakdown$scen)==scens[i])
  
  for(s in 1:5){
    
    if(s==1){
      tmpin = diffs_p1_ann
      seasonin = "annual"
    }
    
    if(s==2){
      tmpin = diffs_p1_DJF
      seasonin = "DJF"
    }
    
    if(s==3){
      tmpin = diffs_p1_MAM
      seasonin = "MAM"
    }
    
    if(s==4){
      tmpin = diffs_p1_JJA
      seasonin = "JJA"
    }
    
    if(s==5){
      tmpin = diffs_p1_SON
      seasonin = "SON"
    }
    
    for(j in 1:length(idxin)){
      GCM = projfilebreakdown$GCM[idxin[j]]
      scen = projfilebreakdown$scen[idxin[j]]
      obs = projfilebreakdown$obs[idxin[j]]
      DS = projfilebreakdown$DS[idxin[j]]
      testsfc1 = list(x=lon-360,y=lat,z=tmpin[,,idxin[j]])
      surface(testsfc1,type="I",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=legadd)
      map("state",add=TRUE) # main=paste("Scen: ",scen," GCM: ",GCM,"\n DS: ", DS," Obs: ",obs,sep="")
      text(-108.85,28.45,labels=GCM,cex=2,pos = 4)
      text(-108.85,27.45,labels=DS,cex=2,pos = 4)
      text(-108.85,26.45,labels=obs,cex=2,pos = 4)
      text(-90.05,26.45,labels=scen,cex=2,pos = 2)
      text(-90.05,27.45,labels=paste("Season: ",seasonin,sep=""),cex=1.75,pos = 2)
    }
  }  
}
dev.off()


######
# 2070-2099
pdf(paste("/data2/3to5/I35/plots/all_mems/",varname,"_2070-2099_Daymetonly.pdf",sep=""),onefile=TRUE,width=15,height=15)
par(mfrow=c(3,3))
for(i in 1:length(scens)){
  
  idxin = which(projfilebreakdown$obs=="Daymet" & as.character(projfilebreakdown$scen)==scens[i])
  
  for(s in 1:5){
    
    if(s==1){
      tmpin = diffs_p2_ann
      seasonin = "annual"
    }
    
    if(s==2){
      tmpin = diffs_p2_DJF
      seasonin = "DJF"
    }
    
    if(s==3){
      tmpin = diffs_p2_MAM
      seasonin = "MAM"
    }
    
    if(s==4){
      tmpin = diffs_p2_JJA
      seasonin = "JJA"
    }
    
    if(s==5){
      tmpin = diffs_p2_SON
      seasonin = "SON"
    }
    
    for(j in 1:length(idxin)){
      GCM = projfilebreakdown$GCM[idxin[j]]
      scen = projfilebreakdown$scen[idxin[j]]
      obs = projfilebreakdown$obs[idxin[j]]
      DS = projfilebreakdown$DS[idxin[j]]
      testsfc1 = list(x=lon-360,y=lat,z=tmpin[,,idxin[j]])
      surface(testsfc1,type="I",zlim=diffcolorbar[[1]],col=diffcolorbar[[3]],breaks=diffcolorbar[[2]],xlab="",ylab="",xaxt='n', yaxt='n',ann=FALSE,add.legend=legadd)
      map("state",add=TRUE) # main=paste("Scen: ",scen," GCM: ",GCM,"\n DS: ", DS," Obs: ",obs,sep="")
      text(-108.85,28.45,labels=GCM,cex=2,pos = 4)
      text(-108.85,27.45,labels=DS,cex=2,pos = 4)
      text(-108.85,26.45,labels=obs,cex=2,pos = 4)
      text(-90.05,26.45,labels=scen,cex=2,pos = 2)
      text(-90.05,27.45,labels=paste("Season: ",seasonin,sep=""),cex=1.75,pos = 2)
    }
  }  
}
dev.off()
