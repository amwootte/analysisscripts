
locname = c("WichitaKS","FortSmithAR","AlbuquerqueNM","SantaFeNM","FarmingtonNM","ShreveportLA","BatonRougeLA","AmarilloTX","SanAntonioTX","DallasTX","McAllenTX","MidlandTX","SpringfieldMO","OklahomaCityOK","TulsaOK","WoodwardOK")

for(i in 1:length(locname)){
  
  command1 = paste("ls /data2/3to5/I35/ts_output/",locname[i],"_pr_*DeltaSD*.csv",sep="")
  command2 = paste("ls /data2/3to5/I35/ts_output/",locname[i],"_pr_*QDM*.csv",sep="")
  command3 = paste("ls /data2/3to5/I35/ts_output/",locname[i],"_pr_*PARM*.csv",sep="")
  filesin = c(system(command1,intern=TRUE),system(command2,intern=TRUE),system(command3,intern=TRUE))
  
  histfiles = filesin[grep("historical",filesin)]
  projfiles = filesin[grep("rcp",filesin)]
  
  filebreakdown = do.call(rbind,strsplit(histfiles,"_",fixed=TRUE))
  filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,7],".",fixed=TRUE))
  filebreakdown3 = data.frame(filebreakdown[,4:6],filebreakdown2[,1])
  names(filebreakdown3) = c("GCM","DS","obs","scen")
  histfilebreakdown = filebreakdown3
  rm(filebreakdown3)
  rm(filebreakdown2)
  rm(filebreakdown)
  
  filebreakdown = do.call(rbind,strsplit(projfiles,"_",fixed=TRUE))
  filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,7],".",fixed=TRUE))
  filebreakdown3 = data.frame(filebreakdown[,4:6],filebreakdown2[,1])
  names(filebreakdown3) = c("GCM","DS","obs","scen")
  projfilebreakdown = filebreakdown3
  rm(filebreakdown3)
  rm(filebreakdown2)
  rm(filebreakdown)
  # histfiles
  histdata = NULL
  for(h in 1:length(histfiles)){
    hdat = read.table(histfiles[h],sep=",",header=TRUE)
    hdat[,2] = hdat[,2]*86400
    hdat[,2] = hdat[,2]/25.4
    hmon = aggregate(hdat[,2],by=list(month = substr(hdat[,1],1,7)),sum,na.rm=TRUE)
    names(hmon)[2] = names(hdat)[2]
    if(h==1){
      histdata = hmon
    } else {
      histdata = merge(histdata,hmon,by="month")
    }
  }
  
  #projfiles
  projdata = NULL
  for(h in 1:length(projfiles)){
    hdat = read.table(projfiles[h],sep=",",header=TRUE)
    hdat[,2] = hdat[,2]*86400
    hdat[,2] = hdat[,2]/25.4
    hmon = aggregate(hdat[,2],by=list(month = substr(hdat[,1],1,7)),sum,na.rm=TRUE)
    names(hmon)[2] = paste(names(hdat)[2],projfilebreakdown$scen[h],sep="_")
    if(h==1){
      projdata = hmon
    } else {
      projdata = merge(projdata,hmon,by="month")
    }
  }
  
  histmean = apply(histdata[,2:ncol(histdata)],1,mean,na.rm=TRUE)
  rcp85idx = grep("rcp85",names(projdata))
  rcp45idx = grep("rcp45",names(projdata))
  rcp26idx = grep("rcp26",names(projdata))
  rcp85mean = apply(projdata[,rcp85idx],1,mean,na.rm=TRUE)
  rcp45mean = apply(projdata[,rcp45idx],1,mean,na.rm=TRUE)
  rcp26mean = apply(projdata[,rcp26idx],1,mean,na.rm=TRUE)
  
  projdata$rcp85mean = rcp85mean
  projdata$rcp45mean = rcp45mean
  projdata$rcp26mean = rcp26mean
  histdata$histmean = histmean
  
  histfilename = paste("/home/woot0002/GPM/",locname[i],"_pr_monthly_historical.csv",sep="")
  projfilename = paste("/home/woot0002/GPM/",locname[i],"_pr_monthly_future.csv",sep="")
  
  write.table(histdata,histfilename,sep=",",row.names=FALSE)
  write.table(projdata,projfilename,sep=",",row.names=FALSE)
  
  message("Finished for location: ",locname[i])
}


