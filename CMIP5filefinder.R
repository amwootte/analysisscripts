#####################
#
# CMIP5 file finder

olddirs = read.table("https://climatedata.oscer.ou.edu/climatedata_full_CMIP5_internal.csv",sep=",",header=TRUE)

SCCASC_climatedata_path_new = OSCER_schooner_path_new = filemoved = c()

for(i in 1:nrow(olddirs)){
  
  message("checking path for file ",olddirs$filename[i])
  
  filecheck = length(system(paste("ls ",olddirs$SCCASC_climatedata_path[i],sep=""),intern=TRUE))
  
  if(filecheck>=1){
    message(olddirs$filename[i]," was not moved")
    filemoved=0
    SCCASC_climatedata_path_new = as.character(olddirs$SCCASC_climatedata_path[i])
    OSCER_schooner_path_new = as.character(olddirs$OSCER_schooner_path[i])
  } else {
    message(olddirs$filename[i]," was moved")
    message("Looking for new file path")
    
    filemoved=1
    filepath_split = strsplit(as.character(olddirs$SCCASC_climatedata_path[i]),split="/",fixed=TRUE)[[1]]
    findcommand = paste("find /",filepath_split[2],"/synda/ -type f -name ",filepath_split[length(filepath_split)],sep="")
    newlocation = system(findcommand,intern=TRUE)
    SCCASC_climatedata_path_new = newlocation
    filepath_split2 = strsplit(newlocation,split="/",fixed=TRUE)[[1]]
    schoonerpath = paste("climate",filepath_split2[2],sep="")
    schoonerpath2 = paste(filepath_split2[3:length(filepath_split2)],collapse="/")
    schoonerpathall = paste("/condo/",schoonerpath,"/",schoonerpath2,sep="")
    OSCER_schooner_path_new = schoonerpathall
    message("new path found!")
  }
  message("Finished path check for file ",i," / ",nrow(olddirs))
  
  tmp = cbind(olddirs[i,],filemoved)
  tmp = cbind(tmp,SCCASC_climatedata_path_new)
  tmp = cbind(tmp,OSCER_schooner_path_new)
  
  if(i==1){
    write.table(tmp,"/home/woot0002/climatedata_full_CMIP5_internal_2022.csv",sep=",",row.names=FALSE,append=FALSE)
  } else {
    write.table(tmp,"/home/woot0002/climatedata_full_CMIP5_internal_2022.csv",sep=",",row.names=FALSE,append=TRUE,col.names=!file.exists("/home/woot0002/climatedata_full_CMIP5_internal_2022.csv"))
  }
  
}

#olddirs$filemoved = filemoved
#olddirs$SCCASC_climatedata_path_new = SCCASC_climatedata_path_new
#olddirs$OSCER_schooner_path_new = OSCER_schooner_path_new

#write.table(olddirs,"/home/woot0002/climatedata_full_CMIP5_internal_2022.csv",sep=",",row.names=FALSE)

