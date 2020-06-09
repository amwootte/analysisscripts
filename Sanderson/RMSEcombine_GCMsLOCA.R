#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

vars = c("tmax","tmax95","tmin","tmin32","pr","pr50")

###
# get RMSEs

RMSEmatlist_LOCA = NRMSEmatlist_LOCA = metadatlist_LOCA = fileslist_LOCA = list()

for(i in 1:length(vars)){
  load(paste("/home/woot0002/RMSEfiles/",vars[i],"_RMSEmats_LOCA.Rdata",sep=""))
  RMSEmatlist_LOCA[[i]]=LOCARMSEmat
  NRMSEmatlist_LOCA[[i]]=normLOCARMSEmat
  metadatlist_LOCA[[i]]=LOCAdat
  fileslist_LOCA[[i]] = LOCAgroup
}

RMSEmatlist_GCM = NRMSEmatlist_GCM = metadatlist_GCM = fileslist_GCM = list()

for(i in 1:length(vars)){
  load(paste("/home/woot0002/RMSEfiles/",vars[i],"_RMSEmats_GCM.Rdata",sep=""))
  RMSEmatlist_GCM[[i]]=GCMRMSEmat
  NRMSEmatlist_GCM[[i]]=normGCMRMSEmat
  metadatlist_GCM[[i]]=GCMdat
  fileslist_GCM[[i]] = GCMgroup
}


###
# linear combination

weight = 1/length(vars)

NRMSEmatout_LOCA = matrix(0,nrow=nrow(RMSEmatlist_LOCA[[1]]),ncol=ncol(RMSEmatlist_LOCA[[1]]))
NRMSEmatout_GCM = matrix(0,nrow=nrow(RMSEmatlist_GCM[[1]]),ncol=ncol(RMSEmatlist_GCM[[1]]))

for(i in 1:length(vars)){
  NRMSEmatout_LOCA = NRMSEmatout_LOCA+NRMSEmatlist_LOCA[[i]]*weight
  NRMSEmatout_GCM = NRMSEmatout_GCM+NRMSEmatlist_GCM[[i]]*weight
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

simnames_LOCA = simnames_GCM = c()
for(i in 1:nrow(metadatlist_LOCA[[1]])){
  if(i<=(nrow(metadatlist_LOCA[[1]])-1)){
    simnames_LOCA[i]= paste(metadatlist_LOCA[[1]][i,1],metadatlist_LOCA[[1]][i,2],metadatlist_LOCA[[1]][i,3],metadatlist_LOCA[[1]][i,4],sep="_")
    simnames_GCM[i]= paste(metadatlist_GCM[[1]][i,1],metadatlist_GCM[[1]][i,2],metadatlist_GCM[[1]][i,3],metadatlist_GCM[[1]][i,4],sep="_")
  } else {
    simnames_LOCA[i]= metadatlist_LOCA[[1]][i,4]
    simnames_GCM[i]= metadatlist_GCM[[1]][i,4]
  }
}

rownames(NRMSEmatout_LOCA) = colnames(NRMSEmatout_LOCA) = simnames_LOCA
del.RMSE_LOCA = get_upper_tri(NRMSEmatout_LOCA)
melted_cormat_LOCA <- melt(del.RMSE_LOCA, na.rm = TRUE)

rownames(NRMSEmatout_GCM) = colnames(NRMSEmatout_GCM) = simnames_GCM
del.RMSE_GCM = get_upper_tri(NRMSEmatout_GCM)
melted_cormat_GCM <- melt(del.RMSE_GCM, na.rm = TRUE)

topend = ceiling(max(c(NRMSEmatout_GCM,NRMSEmatout_LOCA),na.rm=TRUE))

### Heatmap
pdf("NRMSEheatmap_combined_GCM.pdf",height=18,width=18)
ggplot(data = melted_cormat_GCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf("NRMSEheatmap_combined_LOCA.pdf",height=18,width=18)
ggplot(data = melted_cormat_LOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = 1, limit = c(0,topend), breaks=seq(0,topend,by=0.25), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()


###
# save out data

save(list=c("NRMSEmatout_LOCA","LOCAdat","simnames_LOCA"),file="/home/woot0002/RMSEmats_combo_LOCA.Rdata")
save(list=c("NRMSEmatout_GCM","GCMdat","simnames_GCM"),file="/home/woot0002/RMSEmats_combo_GCM.Rdata")







