#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

vars = c("tmax","tmax95","tmin","tmin32","pr","r1mm","rx1day")

###
# get RMSEs

RMSEmatlist = NRMSEmatlist = metadatlist = fileslist = list()

for(i in 1:length(vars)){
  load(paste("/home/woot0002/RMSEfiles/",vars[i],"_RMSEmats.Rdata",sep=""))
  RMSEmatlist[[i]]=RMSEmat
  NRMSEmatlist[[i]]=normRMSEmat
  metadatlist[[i]]=metadat
  fileslist[[i]] = allfiles
}

###
# linear combination

weight = 1/length(vars)

NRMSEmatout = matrix(0,nrow=nrow(RMSEmatlist[[1]]),ncol=ncol(RMSEmatlist[[1]]))

for(i in 1:length(vars)){
  NRMSEmatout = NRMSEmatout+NRMSEmatlist[[i]]*weight
}

####
# heatmap with combined NRMSE

library(ggplot2)
library(reshape2)
get_upper_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

simnames = c()
for(i in 1:nrow(metadatlist[[1]])){
  if(i<=(nrow(metadatlist[[1]])-5)){
    simnames[i]= paste(metadatlist[[1]][i,1],metadatlist[[1]][i,2],metadatlist[[1]][i,3],metadatlist[[1]][i,4],sep="_")
  } else {
    simnames[i]= metadatlist[[1]][i,4]
  }
}

rownames(NRMSEmatout) = colnames(NRMSEmatout) = simnames
del.RMSE = get_upper_tri(NRMSEmatout)
melted_cormat1 <- melt(del.RMSE, na.rm = TRUE)
topend = ceiling(max(NRMSEmatout,na.rm=TRUE))

### Heatmap
pdf("NRMSEheatmap_combined.pdf",height=15,width=15)
ggplot(data = melted_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 5, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 5, hjust = 1))+
  coord_fixed()
dev.off()

###
# save out data
normRMSEmat = NRMSEmatout
metadat = metadatlist[[1]]
save(list=c("normRMSEmat","metadat","simnames"),file="/home/woot0002/RMSEmats_combo.Rdata")







