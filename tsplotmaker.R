####
# ts plotter script

tsplotdata = function(location,varname,unitsout,changetype="absolute",filesep=","){
  
  location = "OKC"
  varname = "tasmax"
  filesep=","
  unitsout = "degrees_F"
  changetype="absolute"
  
  if(varname=="tasmax" | varname=="tasmin" | varname=="pr"){
    position="topleft"
  }
  
  histfiles = system(paste("ls /data2/3to5/I35/ts_output/",location,"_",varname,"*historical.csv",sep=""),intern=TRUE)
  projfiles = system(paste("ls /data2/3to5/I35/ts_output/",location,"_",varname,"*rcp*.csv",sep=""),intern=TRUE)
  
  filebreakdown = do.call(rbind,strsplit(histfiles,"_",fixed=TRUE))
  filebreakdown[,2] = do.call("rbind",strsplit(as.character(filebreakdown[,2]),"/",fixed=TRUE))[,2]
  filebreakdown[,7] = do.call("rbind",strsplit(as.character(filebreakdown[,7]),".",fixed=TRUE))[,1]
  histfilebreakdown = data.frame(filebreakdown[,-1])
  names(histfilebreakdown) = c("city","varname","GCM","DS","obs","period")
  
  filebreakdown = do.call(rbind,strsplit(projfiles,"_",fixed=TRUE))
  filebreakdown[,2] = do.call("rbind",strsplit(as.character(filebreakdown[,2]),"/",fixed=TRUE))[,2]
  filebreakdown[,7] = do.call("rbind",strsplit(as.character(filebreakdown[,7]),".",fixed=TRUE))[,1]
  projfilebreakdown = data.frame(filebreakdown[,-1])
  names(projfilebreakdown) = c("city","varname","GCM","DS","obs","period")
  
  if(varname=="pr"){
    units="mm"
  }
  if(varname=="tasmax" | varname=="tasmin"){
    units="degrees_K"
  }
  
  projfilebreakdown$idx = 1:nrow(projfilebreakdown)
  histfilebreakdown$idx = 1:nrow(histfilebreakdown)
  
  projfilebreakdown_rcp26 = subset(projfilebreakdown,period=="rcp26")
  projfilebreakdown_rcp45 = subset(projfilebreakdown,period=="rcp45")
  projfilebreakdown_rcp85 = subset(projfilebreakdown,period=="rcp85")
  
  histtsdat = NULL
  proj26tsdat = proj45tsdat = proj85tsdat = NULL
  
  for(i in 1:nrow(histfilebreakdown)){
    
    projidx26 = which(projfilebreakdown_rcp26$GCM==histfilebreakdown$GCM[i] & projfilebreakdown_rcp26$DS==histfilebreakdown$DS[i] & projfilebreakdown_rcp26$obs==histfilebreakdown$obs[i])
    projidx45 = which(projfilebreakdown_rcp45$GCM==histfilebreakdown$GCM[i] & projfilebreakdown_rcp45$DS==histfilebreakdown$DS[i] & projfilebreakdown_rcp45$obs==histfilebreakdown$obs[i])
    projidx85 = which(projfilebreakdown_rcp85$GCM==histfilebreakdown$GCM[i] & projfilebreakdown_rcp85$DS==histfilebreakdown$DS[i] & projfilebreakdown_rcp85$obs==histfilebreakdown$obs[i])
    histtmp = read.table(histfiles[i],sep=filesep,header=TRUE)
    projtmp26 = read.table(projfiles[projfilebreakdown_rcp26$idx[projidx26]],sep=filesep,header=TRUE)
    projtmp45 = read.table(projfiles[projfilebreakdown_rcp45$idx[projidx45]],sep=filesep,header=TRUE)
    projtmp85 = read.table(projfiles[projfilebreakdown_rcp85$idx[projidx85]],sep=filesep,header=TRUE)
    
    if(varname=="pr"){
      histtmp[,2] = histtmp[,2]*86400
      projtmp26[,2] = projtmp26[,2]*86400
      projtmp45[,2] = projtmp45[,2]*86400
      projtmp85[,2] = projtmp85[,2]*86400
    }
    
    if(i==1){
      histtsdat = histtmp
      proj26tsdat = projtmp26
      proj45tsdat = projtmp45
      proj85tsdat = projtmp85
    } else {
      histtsdat = merge(histtsdat,histtmp,by="datesin",all=TRUE)
      proj26tsdat = merge(proj26tsdat,projtmp26,by="datesin",all=TRUE)
      proj45tsdat = merge(proj45tsdat,projtmp45,by="datesin",all=TRUE)
      proj85tsdat = merge(proj85tsdat,projtmp85,by="datesin",all=TRUE)
    }
  }
  
  if(varname=="tasmax" | varname=="tasmin"){
    histagg = aggregate(histtsdat[,2:ncol(histtsdat)],by=list(year=as.numeric(substr(histtsdat[,1],1,4))),mean,na.rm=TRUE)
    proj26agg = aggregate(proj26tsdat[,2:ncol(proj26tsdat)],by=list(year=as.numeric(substr(proj26tsdat[,1],1,4))),mean,na.rm=TRUE)
    proj45agg = aggregate(proj45tsdat[,2:ncol(proj45tsdat)],by=list(year=as.numeric(substr(proj45tsdat[,1],1,4))),mean,na.rm=TRUE)
    proj85agg = aggregate(proj85tsdat[,2:ncol(proj85tsdat)],by=list(year=as.numeric(substr(proj85tsdat[,1],1,4))),mean,na.rm=TRUE)
    
    if(units!=unitsout){
      if(unitsout=="degrees_F"){
        histagg[,2:ncol(histagg)] = ((histagg[,2:ncol(histagg)]-273.15)*9/5)+32
        proj26agg[,2:ncol(proj26agg)] = ((proj26agg[,2:ncol(proj26agg)]-273.15)*9/5)+32
        proj45agg[,2:ncol(proj45agg)] = ((proj45agg[,2:ncol(proj45agg)]-273.15)*9/5)+32
        proj85agg[,2:ncol(proj85agg)] = ((proj85agg[,2:ncol(proj85agg)]-273.15)*9/5)+32
      }
      if(unitsout=="degrees_C"){
        histagg[,2:ncol(histagg)] = (histagg[,2:ncol(histagg)]-273.15)
        proj26agg[,2:ncol(proj26agg)] = (proj26agg[,2:ncol(proj26agg)]-273.15)
        proj45agg[,2:ncol(proj45agg)] = (proj45agg[,2:ncol(proj45agg)]-273.15)
        proj85agg[,2:ncol(proj85agg)] = (proj85agg[,2:ncol(proj85agg)]-273.15)
      }
    }
  }
  
  if(varname=="pr"){
    histagg = aggregate(histtsdat[,2:ncol(histtsdat)],by=list(year=as.numeric(substr(histtsdat[,1],1,4))),sum,na.rm=TRUE)
    proj26agg = aggregate(proj26tsdat[,2:ncol(proj26tsdat)],by=list(year=as.numeric(substr(proj26tsdat[,1],1,4))),sum,na.rm=TRUE)
    proj45agg = aggregate(proj45tsdat[,2:ncol(proj45tsdat)],by=list(year=as.numeric(substr(proj45tsdat[,1],1,4))),sum,na.rm=TRUE)
    proj85agg = aggregate(proj85tsdat[,2:ncol(proj85tsdat)],by=list(year=as.numeric(substr(proj85tsdat[,1],1,4))),sum,na.rm=TRUE)
    
    if(units!=unitsout){
      if(unitsout=="in" & var=="pr"){
        histagg[,2:ncol(histagg)] = (histagg[,2:ncol(histagg)]/25.4)
        proj26agg[,2:ncol(proj26agg)] = (proj26agg[,2:ncol(proj26agg)]/25.4)
        proj45agg[,2:ncol(proj45agg)] = (proj45agg[,2:ncol(proj45agg)]/25.4)
        proj85agg[,2:ncol(proj85agg)] = (proj85agg[,2:ncol(proj85agg)]/25.4)
      }
      
    }
  }
    histclimo = apply(histagg[,2:ncol(histagg)],2,mean,na.rm=TRUE)
    
    histanom = histagg
    proj26anom = proj26agg
    proj45anom = proj45agg
    proj85anom = proj85agg
    
    if(changetype=="absolute"){
      histanom[,2:ncol(histanom)] = histanom[,2:ncol(histanom)]-histclimo
      proj26anom[,2:ncol(proj26anom)] = proj26anom[,2:ncol(proj26anom)]-histclimo
      proj45anom[,2:ncol(proj45anom)] = proj45anom[,2:ncol(proj45anom)]-histclimo
      proj85anom[,2:ncol(proj85anom)] = proj85anom[,2:ncol(proj85anom)]-histclimo
    } else {
      histanom[,2:ncol(histanom)] = ((histanom[,2:ncol(histanom)]-histclimo)/histclimo)*100
      proj26anom[,2:ncol(proj26anom)] = ((proj26anom[,2:ncol(proj26anom)]-histclimo)/histclimo)*100
      proj45anom[,2:ncol(proj45anom)] = ((proj45anom[,2:ncol(proj45anom)]-histclimo)/histclimo)*100
      proj85anom[,2:ncol(proj85anom)] = ((proj85anom[,2:ncol(proj85anom)]-histclimo)/histclimo)*100
      unitsout="%"
    }
    
    histyears = 1981:2005
    projyears = 2006:2099
    histanommean = apply(histanom[,2:ncol(histanom)],1,mean,na.rm=TRUE)
    proj26anommean = apply(proj26anom[,2:ncol(proj26anom)],1,mean,na.rm=TRUE)
    proj45anommean = apply(proj45anom[,2:ncol(proj45anom)],1,mean,na.rm=TRUE)
    proj85anommean = apply(proj85anom[,2:ncol(proj85anom)],1,mean,na.rm=TRUE)
    projmeandat = data.frame(proj26anommean,proj45anommean,proj85anommean)
    projmean = apply(projmeandat,1,mean,na.rm=TRUE)
    
    yrange = range(range(histanom[,2:ncol(histanom)],na.rm=TRUE),range(proj26anom[,2:ncol(proj26anom)],na.rm=TRUE),range(proj45anom[,2:ncol(proj45anom)],na.rm=TRUE),range(proj85anom[,2:ncol(proj85anom)],na.rm=TRUE))
    
    ###
    
    pdf(paste("/data2/3to5/I35/plots/tsplots/",location,"_",varname,"_",unitsout,".pdf",sep=""),onefile=TRUE,width=10,height=6)
    
    plot(projmean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l")
    abline(h=0,lty=2)
    
    plot(proj26anommean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l",col="#2166ac")
    lines(proj45anommean~projyears,lwd=2,type="l",col="#92c5de")
    lines(proj85anommean~projyears,lwd=2,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    plot(proj26anommean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l",col="#2166ac")
    lines(proj45anommean~projyears,lwd=2,type="l",col="#92c5de")
    lines(proj85anommean~projyears,lwd=2,type="l",col="#b2182b")
    lines(proj85anom[,6]~projyears,lwd=0.5,type="l",col="#b2182b")
    lines(proj85anom[,12]~projyears,lwd=0.5,type="l",col="#b2182b")
    lines(proj85anom[,18]~projyears,lwd=0.5,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    plot(proj85anommean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l",col="#b2182b")
    lines(proj85anom[,6]~projyears,lwd=0.5,type="l",col="#b2182b")
    lines(proj85anom[,12]~projyears,lwd=0.5,type="l",col="#b2182b")
    lines(proj85anom[,18]~projyears,lwd=0.5,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    plot(proj26anommean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l",col="#2166ac")
    lines(proj45anommean~projyears,lwd=2,type="l",col="#92c5de")
    lines(proj85anommean~projyears,lwd=2,type="l",col="#b2182b")
    for(i in 2:19) lines(proj85anom[,i]~projyears,lwd=0.5,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    plot(proj85anommean~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=2,type="l",col="#b2182b")
    for(i in 2:19) lines(proj85anom[,i]~projyears,lwd=0.5,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
     if(i==1){
       plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=0.5,type="l",col="gray")
       lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
       lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
       } else {
         lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
         lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
         lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
     }
  }
    lines(projmean~projyears,lwd=2,type="l")
    abline(h=0,lty=2)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
      if(i==1){
        plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      } else {
        lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      }
    }
    lines(proj26anommean~projyears,lwd=2,type="l",col="#2166ac")
    lines(proj45anommean~projyears,lwd=2,type="l",col="#92c5de")
    lines(proj85anommean~projyears,lwd=2,type="l",col="#b2182b")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
      if(i==1){
        plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      } else {
        lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      }
    }
    lines(projmean~projyears,lwd=2,type="l")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
      if(i==1){
        plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      } else {
        lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
      }
    }
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
      if(i==1){
        plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,xlim=range(c(histyears,projyears)),lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
        lines(histanom[,(i+1)]~histyears,lwd=0.5,type="l",col="gray")
      } else {
        lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#2166ac")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#92c5de")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="#b2182b")
        lines(histanom[,(i+1)]~histyears,lwd=0.5,type="l",col="gray")
      }
    }
    lines(proj26anommean~projyears,lwd=2,type="l",col="#2166ac")
    lines(proj45anommean~projyears,lwd=2,type="l",col="#92c5de")
    lines(proj85anommean~projyears,lwd=2,type="l",col="#b2182b")
    lines(histanommean~histyears,lwd=2,type="l")
    abline(h=0,lty=2)
    legend(position,c("RCP 2.6","RCP 4.5","RCP 8.5"),col=c("#2166ac","#92c5de","#b2182b"),lwd=4)
    
    ###
    
    for(i in 1:(ncol(proj26anom)-1)){
      if(i==1){
        plot(proj26anom[,(i+1)]~projyears,main = paste("Projected Change in ",varname," for ",location,sep=""),xlab="year",ylab=unitsout,ylim=yrange,xlim=range(c(histyears,projyears)),lwd=0.5,type="l",col="gray")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
      } else {
        lines(proj26anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
        lines(proj45anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
        lines(proj85anom[,(i+1)]~projyears,lwd=0.5,type="l",col="gray")
      }
      lines(histanom[,(i+1)]~histyears,lwd=0.5,type="l",col="gray")
    }
    lines(histanommean~histyears,lwd=2,type="l")
    lines(projmean~projyears,lwd=2,type="l")
    abline(h=0,lty=2)
    dev.off()
}

  
  
  



