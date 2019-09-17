####
# bar plotter script

barplotdata = function(datafiles,filesep=","){

#varname = "heatwaves"
#datafiles = "/data2/3to5/I35/area_output/heatwaves_PottowotomieCounty_absolute_2041-2070_ann.csv"
#datafiles = paste("/data2/3to5/I35/point_output/",c("pr50_Anadarko_absolute_2071-2099_ann.csv","pr50_Bartlesville_absolute_2071-2099_ann.csv","pr50_Miami_absolute_2071-2099_ann.csv","pr50_Guymon_absolute_2071-2099_ann.csv"),sep="")
#filesep=","

  #datafiles = "/data2/3to5/I35/point_output/pr50_Anadarko_absolute_2071-2099_ann.csv,/data2/3to5/I35/point_output/pr50_Bartlesville_absolute_2071-2099_ann.csv,/data2/3to5/I35/point_output/pr50_Miami_absolute_2071-2099_ann.csv,/data2/3to5/I35/point_output/pr50_Guymon_absolute_2071-2099_ann.csv"
  
  datafiles = do.call("c",strsplit(datafiles,",",fixed=TRUE))
filebreakdown = do.call(rbind,strsplit(datafiles,"/",fixed=TRUE))
filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,ncol(filebreakdown)],"_",fixed=TRUE))
filebreakdown3 = do.call(rbind,strsplit(filebreakdown2[,4],"-",fixed=TRUE))
filebreakdown4 = do.call(rbind,strsplit(filebreakdown2[,ncol(filebreakdown2)],".",fixed=TRUE))


if(length(datafiles)>1){
  message("running for multiple files")
  datafilebreakdown = data.frame(filebreakdown2[,1],filebreakdown3,filebreakdown2[,2:3],filebreakdown4[,1])
  message("NUmber of columns ",ncol(datafilebreakdown))
  names(datafilebreakdown) = c("varname","startyear","endyear","area","units","period")
  datafilebreakdown[,1]=as.character(datafilebreakdown[,1])
  datafilebreakdown[,2]=as.character(datafilebreakdown[,2])
  datafilebreakdown[,3]=as.character(datafilebreakdown[,3])
  datafilebreakdown[,4]=as.character(datafilebreakdown[,4])
  datafilebreakdown[,5]=as.character(datafilebreakdown[,5])
  datafilebreakdown[,6]=as.character(datafilebreakdown[,6])
} else {
  message("running for one file")
  #filebreakdown3 = matrix(c(varname,filebreakdown2,filebreakdown[,3:5]),ncol=6,nrow=1)
  datafilebreakdown = data.frame(filebreakdown2[,1],filebreakdown3,filebreakdown2[,2:3],filebreakdown4[,1])
  names(datafilebreakdown) = c("varname","startyear","endyear","area","units","period")
  datafilebreakdown[,1]=as.character(datafilebreakdown[,1])
  datafilebreakdown[,2]=as.character(datafilebreakdown[,2])
  datafilebreakdown[,3]=as.character(datafilebreakdown[,3])
  datafilebreakdown[,4]=as.character(datafilebreakdown[,4])
  datafilebreakdown[,5]=as.character(datafilebreakdown[,5])
  datafilebreakdown[,6]=as.character(datafilebreakdown[,6])
}

for(r in 1:nrow(datafilebreakdown)){
if(datafilebreakdown$varname[r]=="pr50" | datafilebreakdown$varname[r]=="pr25" | datafilebreakdown$varname[r]=="tmax100" | datafilebreakdown$varname[r]=="gsl" | datafilebreakdown$varname[r]=="frd" | datafilebreakdown$varname[r]=="mdrn" | datafilebreakdown$varname[r]=="tmax95" | datafilebreakdown$varname[r]=="cdd"){
  datafilebreakdown$units[r]="days"
}

if(datafilebreakdown$varname[r]=="rx1day" | datafilebreakdown$varname[r]=="rx5day"){
  datafilebreakdown$units[r]="mm"
}


if(datafilebreakdown$varname[r]=="tasmax" | datafilebreakdown$varname[r]=="tasmin"){
  datafilebreakdown$units[r]="degrees_K"
}

if(datafilebreakdown$varname[r]=="heatwaves"){
  datafilebreakdown$units[r]="events per year"
}
}

library(ggplot2)

if(length(datafiles)==1){
  message("Reading one input file")
  datain = read.table(datafiles[1],header=TRUE,sep=filesep)
  datainmean = aggregate(datain$projchange,by=list(scen=datain$scen),mean,na.rm=TRUE)
  datainsd = aggregate(datain$projchange,by=list(scen=datain$scen),sd,na.rm=TRUE)
  dataindat = merge(datainmean,datainsd,by="scen")
  names(dataindat) = c("scen","mean","sd")
  dataindat$ci = (dataindat$sd/sqrt(9))*1.96
  dataindat$var = datafilebreakdown$varname[1]
  
  datarange= range(c(dataindat$mean-dataindat$ci,dataindat$mean+dataindat$ci))
  if(all(datarange>0)==TRUE){
    datarange[1]=0
  }
  if(all(datarange<0)==TRUE){
    datarange[2]=0
  }
  
  ggplot(dataindat, aes(x=var, y=mean, fill=scen)) + 
    geom_bar(position=position_dodge(), stat="identity")+
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    scale_fill_manual(values=c("#2166ac", "#92c5de", "#b2182b"))+ylim(datarange)+
    ggtitle(paste(datafilebreakdown$area[1]," ",datafilebreakdown$varname[1],"\n Projected Change ",datafilebreakdown$startyear[1],"-",datafilebreakdown$endyear[1],sep=""))+xlab(datafilebreakdown$varname[1])+ylab(paste("projected change",datafilebreakdown$units[1],sep=" "))
  outfile = paste("/data2/3to5/I35/plots/barplots/",datafilebreakdown$varname[1],"_barplot_",datafilebreakdown$area[1],"_",datafilebreakdown$units[1],"_",datafilebreakdown$startyear[1],"-",datafilebreakdown$endyear[1],".pdf",sep="")
  ggsave(outfile,width=5,height=10)
  
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
    dataindat$var = datafilebreakdown$varname[i]
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
      ggtitle(paste(datafilebreakdown$area[i]," ",datafilebreakdown$varname[i]," Projected Change ",datafilebreakdown$startyear[i],"-",datafilebreakdown$endyear[i],sep=""))+xlab(datafilebreakdown$varname[i])+ylab(paste("projected change",datafilebreakdown$units[i],sep=" "))
    outfile = paste("/data2/3to5/I35/plots/barplots/",datafilebreakdown$varname[i],"_barplot_",datafilebreakdown$area[i],"_",datafilebreakdown$units[i],"_",datafilebreakdown$startyear[i],"-",datafilebreakdown$endyear[i],"_",datafilebreakdown$period[i],".pdf",sep="")
    ggsave(outfile,width=5,height=10)
  }
  
}
  
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
    make_option(c('-d', "--datafiles"), action="store", default='na',
                dest='datafiles',
                help=paste("The variable of interest (tasmax, tasmin, pr,tmax95,tmin32,pr25) in the input file.", 
                           "If not present, an error will be thrown.")),
    make_option(c("-f", "--filesep"), action="store", default="comma",
                dest='filesep',
                help=paste("Delimiter for the list of files passed to -d"))
  )
  
  description = paste('Given the files provided, this script will create barplots. ', 
                      "The files should be of a format (both file name and structure) ",
                      "similar to those resulting from point_calcs.R. ")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                   "not required to specify strings.")
  usage = paste("usage: %prog -d datafiles", 
                "[-f filesep]","[-h --help]")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

#####

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)

parsed.args <- ParseArgs(args)

if(parsed.args$filesep=="comma") filesep=","
if(parsed.args$filesep=="space") filesep=" "

barplotdata(datafiles = parsed.args$datafiles,filesep=filesep)


