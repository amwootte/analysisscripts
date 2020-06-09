###
# Attribute adder

DS = "PARM"
var = "pr"
period = "future"

filelist = system(paste("ls /data/static_web/3to5_GDP/",DS,"/",period,"/",var,"/",var,"_day*.nc",sep=""),intern=TRUE)
#setwd(paste("/data/static_web/3to5_GDP/",DS,"/",period,"/",var,"/",sep=""))

if(period=="historical"){
  filebreakdown = do.call(rbind,strsplit(filelist,"_",fixed=TRUE))
  filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,9],"-",fixed=TRUE))
  filebreakdown3 = data.frame(var,filebreakdown2,filebreakdown[,4:7])
  filebreakdown3$GCM =rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=3),length(unique(filebreakdown3[,6])))
  filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),3),length(unique(filebreakdown3[,6])))
  filebreakdown3 = filebreakdown3[,-c(2,3)]
  names(filebreakdown3) = c("var","tempres","code","scen","experiment","GCM","obs")
  rm(filebreakdown2)
  rm(filebreakdown)
} else {
  filebreakdown = do.call(rbind,strsplit(filelist,"_",fixed=TRUE))
  filebreakdown2 = do.call(rbind,strsplit(filebreakdown[,9],"-",fixed=TRUE))
  filebreakdown3 = data.frame(var,filebreakdown2,filebreakdown[,4:7])
  filebreakdown3$GCM = rep(rep(c("CCSM4","MIROC5","MPI-ESM-LR"),each=9),length(unique(filebreakdown3[,4])))
  filebreakdown3$obs = rep(rep(c("Daymet","Livneh","PRISM"),9),length(unique(filebreakdown3[,4])))
  filebreakdown3 = filebreakdown3[,-c(2,3)]
  names(filebreakdown3) = c("var","tempres","code","scen","experiment","GCM","obs")
  rm(filebreakdown2)
  rm(filebreakdown)
}

commandbase = "ncatted -O -a"

if(DS=="PARM") DSpart = "Piecewise Asychronous Regression Model (PARM)"
if(DS=="EDQM" & (var=="tasmax" | var=="tasmin")) DSpart = "Equi-Distant Quantile Mapping (EDQM)"
if(DS=="EDQM" & var=="pr") DSpart = "Equi-Ratio Quantile Mapping (ERQM)"
if(DS=="DeltaSD" & (var=="tasmax" | var=="tasmin")) DSpart = "Delta Method - Additive (DeltaSD)"
if(DS=="DeltaSD" & var=="pr") DSpart = "Delta Method - Ratio (DeltaSD)"

commandsout = c()

for(i in 1:length(filelist)){
   
  message("Working on file ",i," / ",length(filelist))
  
  # set training data line
  if(filebreakdown3$obs[i]=="Daymet") tdpart = "Daymet v. 2.1"
  if(filebreakdown3$obs[i]=="Livneh") tdpart = "Livneh v. 1.2"
  if(filebreakdown3$obs[i]=="PRISM") tdpart = "PRISM AN81d"
  
  # training data addition 
  command1 = paste(commandbase," training_data,global,c,c,'",tdpart,"' ",filelist[i],sep="")
  
  # Downscaling addition
  command2 = paste(commandbase," DS,global,c,c,'",DSpart,"' ",filelist[i],sep="")
  
  # GCM addition
  command3 = paste(commandbase," GCM,global,c,c,'",filebreakdown3$GCM[i],"' ",filelist[i],sep="")
  
  if(DS=="DeltaSD"){
    extracommand = paste(commandbase," DS_notes,global,c,c,'To compute anomalies with DeltaSD use observations from the training data for this file as the historical baseline' ",filelist[i],sep="")
  }
  
  # scenario addition
  command4 = paste(commandbase," scenario,global,c,c,'",filebreakdown3$scen[i],"' ",filelist[i],sep="")
  
  #message("Running commands")
  #system(command1)
  #system(command2)
  #system(command3)
  #system(command4)
  #if(DS=="DeltaSD") system(extracommand)
  if(DS!="DeltaSD"){
    commands = c(command1,"wait",command2,"wait",command3,"wait",command4,"wait")
  } else {
    commands = c(command1,"wait",command2,"wait",command3,"wait",command4,"wait",extracommand,"wait")
  }
  
  commandsout = c(commandsout,commands)
}

write(commandsout,paste("/home/woot0002/commands_",var,"_",DS,"_",period,".csh",sep=""),sep=" ")
