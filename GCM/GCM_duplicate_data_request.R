filetable = read.table("duplicate_files.csv",sep=",",header=TRUE,colClasses="character")

filenamesplit = do.call(rbind,strsplit(filetable$filename,"_",fixed=TRUE))
monidx = which(filenamesplit[,2]=="Amon" | filenamesplit[,2]=="cfMon")

filenamesplitin = filenamesplit[monidx,]
checkfiletable = filetable[monidx,]

uniqueids = unique(checkfiletable$unique_tag)

library(ncdf4)
diff_check = c()

for(i in 1:length(uniqueids)){
  unqidx = which(checkfiletable$unique_tag==uniqueids[i])
  tmptable = checkfiletable[unqidx,]
  tmpsplit = filenamesplitin[unqidx,]
  
  datalist = list()
  if(nrow(tmptable)>1){
  for(j in 1:nrow(tmptable)){
    checkfile = system(paste("ls ",tmptable$local_file[j],sep=""),intern=TRUE)
    if(length(checkfile)>0){
    test1 = nc_open(tmptable$local_file[j])
    datalist[[j]]=ncvar_get(test1,tmpsplit[j,1])
    nc_close(test1)
    } else {
      datalist[[j]]=NA
    }
  }
  }
  
  if(nrow(tmptable)==1) diff_check[i]=0
  
  if(nrow(tmptable)==2){
    diff_check[i] = mean(datalist[[2]]-datalist[[1]],na.rm=TRUE)
  }
  
  if(nrow(tmptable)>=3){
    tcheck = c()
    for(k in 2:nrow(tmptable)){
      tcheck[(k-1)] = mean(datalist[[k]]-datalist[[1]],na.rm=TRUE)
    }
    diff_check[i] = mean(tcheck,na.rm=TRUE)
  }
  rm(datalist)
  gc()
  message("completed check for unique group ",i," / ",length(uniqueids))
}

difference = data.frame(uniqueids,diff_check)
write.table(difference,"Duplicate_monthly_GCM_data_check.csv",sep=",",row.names=FALSE)





