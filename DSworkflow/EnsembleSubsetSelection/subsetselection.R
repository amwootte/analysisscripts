source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)
library(modi)
library(parallel)
library(foreach)
library(doParallel)

#####
# Prep filenames and read.

enssize = 5
ensgroup = "CMIP6"
mask = "SGP-NCA"
#climotype="annual"

RMSEfiles = system(paste("ls /home/woot0002/RCMES/*_",ensgroup,"_RMSEdat_enssize",enssize,"_dom",mask,"_*.Rdata",sep=""),intern=TRUE)
FRCfiles = system(paste("ls /home/woot0002/RCMES/*_",ensgroup,"_FRCdat_enssize",enssize,"_dom",mask,"_*.Rdata",sep=""),intern=TRUE)

RMSEfilestab1 = do.call("rbind",strsplit(RMSEfiles,"/",fixed=TRUE))
RMSEfilestab2 = do.call("rbind",strsplit(RMSEfilestab1[,ncol(RMSEfilestab1)],"_",fixed=TRUE))
RMSEfilestab = data.frame(RMSEfilestab2)
RMSEfilestab = RMSEfilestab[,1:6]
names(RMSEfilestab) = c("variable","ensgroup","metric","ensemblesize","domain","climotype")
#RMSEfilestab$climotype = substr(as.character(RMSEfilestab$climotype),1,nchar(as.character(RMSEfilestab$climotype))-6)

FRCfilestab1 = do.call("rbind",strsplit(FRCfiles,"/",fixed=TRUE))
FRCfilestab2 = do.call("rbind",strsplit(FRCfilestab1[,ncol(FRCfilestab1)],"_",fixed=TRUE))
FRCfilestab = data.frame(FRCfilestab2)
names(FRCfilestab) = c("variable","ensgroup","metric","ensemblesize","domain","climotype")
FRCfilestab = FRCfilestab[,1:6]
#FRCfilestab$climotype = substr(as.character(FRCfilestab$climotype),1,nchar(as.character(FRCfilestab$climotype))-6)

RMSEfilestab$varclimo = paste(RMSEfilestab$variable,RMSEfilestab$climotype,sep="-")
FRCfilestab$varclimo = paste(FRCfilestab$variable,FRCfilestab$climotype,sep="-")

idxin = which(RMSEfilestab$varclimo %in% FRCfilestab$varclimo)

RMSEfilesin = RMSEfiles[idxin]
FRCfilesin = FRCfiles

#RMSEfilesin = RMSEfilesin[-grep("r1mm",RMSEfilesin)]
#FRCfilesin = FRCfilesin[-grep("r1mm",FRCfilesin)]

#####
# load RMSEs and FRCss for each variable

alldat = list()
for(i in 1:length(RMSEfilesin)){
  load(RMSEfilesin[i])
  load(FRCfilesin[i])
  tmpdat = data.frame(idxs,FRCs,RMSEs,NRMSEs)
  names(tmpdat) = c("idx","frc","nfrc","rmse","nrmse")
  alldat[[i]] = tmpdat
}

######
# calculate Euclidean Distance values for multivariate (or univariate calcs)

alldatmv = alldat[[1]]

###
# calculate distance for individual variables
indivvar = FALSE
if(indivvar == TRUE){
for(i in 1:length(alldat)){
  alldat[[i]]$optvals = optimfunc(alldat[[i]]$nrmse,alldat[[i]]$frc)
  if(i==1){
    alldatmv$optvals = alldat[[i]]$optvals*(1/length(RMSEfilesin))
  } else {
    alldatmv$optvals = alldatmv$optvals+alldat[[i]]$optvals*(1/length(RMSEfilesin))
  }
}
}

####
# Determine optimal combination

mindist = min(alldatmv$optvals)
alldatmvmin = alldatmv[which(alldatmv$optvals==mindist),]


pdf(paste("/home/woot0002/RCMES/allvars_",ensgroup,"_optselect_enssize",enssize,"_dom",mask,"_images.pdf",sep=""),onefile=TRUE,height=6,width=6)

for(i in 1:length(RMSEfilesin)){
  
  tmp1 = do.call("rbind",strsplit(RMSEfilesin[i],"/",fixed=TRUE))
  tmp2 = do.call("rbind",strsplit(tmp1[,ncol(tmp1)],"_",fixed=TRUE))
  varname = tmp2[1,1]
  ctype = substr(tmp2[1,ncol(tmp2)],1,nchar(tmp2[1,ncol(tmp2)])-6)
  
  plot(frc~nrmse,data=alldat[[i]],ylim=c(0,1),xlim=c(0,max(nrmse)),xlab="NRMSE",ylab="FRC",main=paste("FRC vs. NRMSE - ",varname,"\nEnsemble size: ",enssize," Domain: ",mask," ClimoType: ",ctype,sep=""))
  maxidxonevar = which(alldat[[i]]$optvals == min(alldat[[i]]$optvals,na.rm=TRUE))
  maxidxmv = which(alldatmv$optvals == min(alldatmv$optvals,na.rm=TRUE))
  maxidxmv_cold = alldatmv_cold$idx[which(alldatmv_cold$optvals == min(alldatmv_cold$optvals,na.rm=TRUE))]
  points(alldat[[i]]$frc[maxidxmv]~alldat[[i]]$nrmse[maxidxmv],col="purple",pch=19)
  points(alldat[[i]]$frc[maxidxonevar]~alldat[[i]]$nrmse[maxidxonevar],col="red",pch=19)
  legend("bottomright",legend=c("All options",paste(varname," Optimal option: ",maxidxonevar,sep=""),paste("mv Optimal option: ",maxidxmv,sep="")),pch=c(1,19,19),col=c("black","red","purple"))
}

dev.off()

#####

# Determine GCMs selected in multivariate and for other variables - all models considered

GCMmat = matrix(NA,nrow=length(alldat)+1,ncol=enssize)
minidxvals = c()
varsin = c()
optvalin = c()
ctype = c()
for(i in 1:length(alldat)){
  tmp = alldat[[i]][order(alldat[[i]]$optvals,decreasing=FALSE),]
  GCMmat[i,]=varnames[(combosused[,tmp$idx[1]]-1)]
  minidxvals[i] = tmp$idx[1]
  
  tmp1 = do.call("rbind",strsplit(RMSEfilesin[i],"/",fixed=TRUE))
  tmp2 = do.call("rbind",strsplit(tmp1[,ncol(tmp1)],"_",fixed=TRUE))
  varsin[i] = tmp2[1,1]
  optvalin[i] = tmp$optvals[1]
  ctype[i] = tmp2[1,(ncol(tmp2)-1)] #substr(tmp2[1,(ncol(tmp2)-1)],1,nchar(tmp2[1,(ncol(tmp2)-1)])-6)
}

mindist = min(alldatmv$optvals)
alldatmvmin = alldatmv[which(alldatmv$optvals==mindist),]
mvidx = which(alldatmv$optvals==mindist)
GCMmat[nrow(GCMmat),]=varnames[(combosused[,mvidx]-1)]

varsin = c(varsin,"multivariate")
minidxvals = c(minidxvals,mvidx)
optvalin = c(optvalin,alldatmvmin$optvals[1])
ctype = c(ctype,NA)

GCMmat = data.frame(GCMmat)
names(GCMmat) = paste("GCM",1:enssize,sep="-")
for(i in 1:enssize) GCMmat[,i] = as.character(GCMmat[,i])

GCMmat$variable = varsin
GCMmat$minidx = minidxvals
GCMmat$optvals = optvalin
GCMmat$ctype = ctype

####
# save results - all models

fileout = paste("/home/woot0002/RCMES/allvars_",ensgroup,"_optselect_enssize",enssize,"_dom",mask,".Rdata",sep="")
save(list=c("varnames","combosused","alldat","alldatmv","enssize","GCMmat"),file=fileout)


