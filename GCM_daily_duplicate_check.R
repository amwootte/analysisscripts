filetable = read.table("duplicate_files.csv",sep=",",header=TRUE,colClasses="character")

filenamesplit = do.call(rbind,strsplit(filetable$filename,"_",fixed=TRUE))
monidx = which(filenamesplit[,2]=="day" | filenamesplit[,2]=="cfDay")

filenamesplitin = filenamesplit[monidx,]
checkfiletable = filetable[monidx,]

uniqueids = unique(checkfiletable$unique_tag)

Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/usr/local/anaconda/bin",sep=":"))

library(ncdf4)
diff_check = c()

for(i in 1:length(uniqueids)){
  unqidx = which(checkfiletable$unique_tag==uniqueids[i])
  tmptable = checkfiletable[unqidx,]
  tmpsplit = filenamesplitin[unqidx,]
  
  if(nrow(tmptable)==1) diff_check[i]=0
  
  if(nrow(tmptable)==2){
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[2]," filediff1.nc",sep="")
    system(command)
    test1 = nc_open("filediff1.nc")
    diff_check[i] = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff1.nc")
  }
  
  if(nrow(tmptable)==3){
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[2]," filediff1.nc",sep="")
    system(command)
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[3]," filediff2.nc",sep="")
    system(command)
    
    test1 = nc_open("filediff1.nc")
    tc1 = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff1.nc")
    
    test1 = nc_open("filediff2.nc")
    tc2 = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff2.nc")
    
    diff_check[i] = mean(c(tc1,tc2))
  }
  
  if(nrow(tmptable)==4){
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[2]," filediff1.nc",sep="")
    system(command)
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[3]," filediff2.nc",sep="")
    system(command)
    command = paste("ncdiff -v ",tmpsplit[1,1]," ",tmptable$local_file[1]," ",tmptable$local_file[4]," filediff3.nc",sep="")
    system(command)
    
    test1 = nc_open("filediff1.nc")
    tc1 = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff1.nc")
    
    test1 = nc_open("filediff2.nc")
    tc2 = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff2.nc")
    
    test1 = nc_open("filediff3.nc")
    tc3 = mean(ncvar_get(test1,tmpsplit[1,1]))
    nc_close(test1)
    system("rm filediff3.nc")
    
    diff_check[i] = mean(c(tc1,tc2,tc3))
  }
  
  gc()
  message("completed check for unique group ",i," / ",length(uniqueids))
}

difference = data.frame(uniqueids,diff_check)
write.table(difference,"Duplicate_daily_GCM_data_check.csv",sep=",",row.names=FALSE)





