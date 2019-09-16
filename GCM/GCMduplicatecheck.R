# GCM duplicates check

filetable = read.table("duplicate_files.csv",sep=",",header=TRUE)

uniqueids = unique(filetable$unique_tag)

diffcheck = list()

for(i in 18464:length(uniqueids)){
  
  tmptable = subset(filetable,unique_tag==uniqueids[i])
  for(j in 1:nrow(tmptable)){
    system(paste("ncdump -h ",tmptable$local_file[j]," > tmpinfo_",j,".txt",sep=""))
  }
  
  if(nrow(tmptable)==2){
    diffcheck[[i]] = length(system("diff tmpinfo_1.txt tmpinfo_2.txt",intern=TRUE))
  }
  
  if(nrow(tmptable)>=3){
    dcheck = c()
    for(k in 2:nrow(tmptable)){
      dcheck= c(dcheck,length(system(paste("diff tmpinfo_1.txt tmpinfo_",k,".txt",sep=""),intern=TRUE)))
    }
    diffcheck[[i]] = dcheck
  }
  
  for(m in 1:nrow(tmptable)){
    system(paste("rm tmpinfo_",m,".txt",sep=""))
  }
  
  message("Finished checking uniqueid ",i," / ",length(uniqueids))
  
}

diffchecklengths = sapply(diffcheck,length)
checktheseids = which(diffchecklengths>1)

differentheaders = data.frame(uniqueids[checktheseids],diffchecklengths[checktheseids])
names(uniqueids) = c("uniqueids","number of differences found")

write.table(differentheaders,"GCM_duplicate_different_headers.csv",sep=",",row.names=FALSE)








