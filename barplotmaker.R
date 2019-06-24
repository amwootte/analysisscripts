####
# bar plotter script

barplotdata = function(datafiles,varname,filesep=","){

varname = "mdrn"
datafiles = paste("/data2/3to5/I35/point_output/",c("mdrn_Albuquerque_absolute_2041-2070_ann.csv","mdrn_Chama_absolute_2041-2070_ann.csv","mdrn_Espanola_absolute_2041-2070_ann.csv","mdrn_Socorro_absolute_2041-2070_ann.csv"),sep="")
filesep=","

filebreakdown = do.call(rbind,strsplit(datafiles,"_",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,5],"-",fixed=TRUE))
filebreakdown3 = data.frame(filebreakdown[,1:2],filebreakdown2,filebreakdown[,3:6])
filebreakdown3[,2] = do.call("rbind",strsplit(as.character(filebreakdown3[,2]),"/",fixed=TRUE))[,2]
filebreakdown3[,8] = do.call("rbind",strsplit(as.character(filebreakdown3[,8]),".",fixed=TRUE))[,1]
datafilebreakdown = filebreakdown3[,-c(1,7)]
names(datafilebreakdown) = c("varname","startyear","endyear","city","units","period")

if(varname=="pr50" | varname=="tmax100" |varname=="gsl" | varname=="frd" | varname=="mdrn" | varname=="tmax95" | varname=="cdd"){
  datafilebreakdown$units="days"
}

if(varname=="rx1day" | varname=="rx5day"){
  datafilebreakdown$units="mm"
}


if(varname=="tasmax" | varname=="tasmin"){
  datafilebreakdown$units="degrees_K"
}


library(ggplot2)

if(length(datafiles)==1){

} else {
  
  message("Reading multiple input files")
  message("Working with ",length(datafiles)," total")
  dataranges = c()
  for(i in 1:length(datafiles)){
    message("working on file ",datafiles[i])
    datain = read.table(datafiles[i],header=TRUE,sep=filesep)
    datainmean = aggregate(datain$projchange,by=list(scen=datain$scen),mean,na.rm=TRUE)
    datainsd = aggregate(datain$projchange,by=list(scen=datain$scen),sd,na.rm=TRUE)
    dataindat = merge(datainmean,datainsd,by="scen")
    names(dataindat) = c("scen","mean","sd")
    dataindat$ci = (dataindat$sd/sqrt(9))*1.96
    dataindat$var = varname
    assign(paste("databounds.",i,sep=""),dataindat)
    dataranges = c(dataranges,dataindat$mean-dataindat$ci,dataindat$mean+dataindat$ci)
  }
  
  datarange= range(dataranges)
  if(all(datarange>0)==TRUE){
    datarange[1]=0
  }
  if(all(datarange<0)==TRUE){
    datarange[2]=0
  }
  
  
  for(i in 1:length(datafiles)){
    tmp = eval(parse(text=paste("databounds.",i,sep="")))
    ggplot(tmp, aes(x=var, y=mean, fill=scen)) + 
      geom_bar(position=position_dodge(), stat="identity")+
      geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9))+
      scale_fill_manual(values=c("#2166ac", "#92c5de", "#b2182b"))+ylim(datarange)+
      ggtitle(paste(datafilebreakdown$city[i]," ",varname," Projected Change ",datafilebreakdown$startyear[i],"-",datafilebreakdown$endyear[i],sep=""))+xlab(varname)+ylab(paste("projected change",datafilebreakdown$units[i],sep=" "))
    outfile = paste("/data2/3to5/I35/plots/barplots/",varname,"_barplot_",datafilebreakdown$city[i],"_",datafilebreakdown$units[i],"_",datafilebreakdown$startyear[i],"-",datafilebreakdown$endyear[i],"_",datafilebreakdown$period[i],".pdf",sep="")
    ggsave(outfile,width=5,height=10)
  }

  
}
  
}

######


