#####################
#
# RMSE calculator

source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)

#vars = c("tmax","tmax95","tmin","tmin32","pr","pr50")
vars = c("pr")
type="ann"
plotname = "pr"
statemask = "new mexico"
#vars=c("pr","pr50")
###
# get RMSEs

NRMSEhmatlist_LOCA = NRMSEcmatlist_LOCA = metadatlist_LOCA = fileslist_LOCA = list()

for(i in 1:length(vars)){
  load(paste("/home/woot0002/RMSEfiles/",vars[i],"_RMSEmats_combo_LOCA_v4_",type,"_",statemask,".Rdata",sep=""))
  NRMSEhmatlist_LOCA[[i]]=normLOCAhRMSEmat
  NRMSEcmatlist_LOCA[[i]]=normLOCAcRMSEmat
  metadatlist_LOCA[[i]]=LOCAhdat
  fileslist_LOCA[[i]] = LOCAgroup
}

NRMSEhmatlist_GCM = NRMSEcmatlist_GCM = metadatlist_GCM = fileslist_GCM = list()

for(i in 1:length(vars)){
  load(paste("/home/woot0002/RMSEfiles/",vars[i],"_RMSEmats_combo_GCM_v4_",type,"_",statemask,".Rdata",sep=""))
  NRMSEhmatlist_GCM[[i]]= normGCMhRMSEmat
  NRMSEcmatlist_GCM[[i]]=normGCMcRMSEmat
  metadatlist_GCM[[i]]=GCMhdat
  fileslist_GCM[[i]] = GCMgroup
}

###
# linear combination

weight = 1/length(vars)

NRMSEhmatout_LOCA = matrix(0,nrow=nrow(NRMSEhmatlist_LOCA[[1]]),ncol=ncol(NRMSEhmatlist_LOCA[[1]]))
NRMSEhmatout_GCM = matrix(0,nrow=nrow(NRMSEhmatlist_GCM[[1]]),ncol=ncol(NRMSEhmatlist_GCM[[1]]))

NRMSEcmatout_LOCA = matrix(0,nrow=nrow(NRMSEcmatlist_LOCA[[1]]),ncol=ncol(NRMSEcmatlist_LOCA[[1]]))
NRMSEcmatout_GCM = matrix(0,nrow=nrow(NRMSEcmatlist_GCM[[1]]),ncol=ncol(NRMSEcmatlist_GCM[[1]]))

for(i in 1:length(vars)){
  NRMSEhmatout_LOCA = NRMSEhmatout_LOCA+NRMSEhmatlist_LOCA[[i]]*weight
  NRMSEhmatout_GCM = NRMSEhmatout_GCM+NRMSEhmatlist_GCM[[i]]*weight
  
  NRMSEcmatout_LOCA = NRMSEcmatout_LOCA+NRMSEcmatlist_LOCA[[i]]*weight
  NRMSEcmatout_GCM = NRMSEcmatout_GCM+NRMSEcmatlist_GCM[[i]]*weight
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

rownames(NRMSEhmatout_LOCA) = colnames(NRMSEhmatout_LOCA) = simnames_LOCA
del.RMSE_LOCA = get_upper_tri(NRMSEhmatout_LOCA)
melted_cormat_LOCA <- melt(del.RMSE_LOCA, na.rm = TRUE)

rownames(NRMSEhmatout_GCM) = colnames(NRMSEhmatout_GCM) = simnames_GCM
del.RMSE_GCM = get_upper_tri(NRMSEhmatout_GCM)
melted_cormat_GCM <- melt(del.RMSE_GCM, na.rm = TRUE)

topend = ceiling(max(c(NRMSEhmatout_GCM,NRMSEhmatout_LOCA,NRMSEcmatout_GCM,NRMSEcmatout_LOCA),na.rm=TRUE))

### Heatmap
pdf(paste("/home/woot0002/DS_ind/NRMSEheatmap_",plotname,"_GCM_hist_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormat_GCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf(paste("/home/woot0002/DS_ind/NRMSEheatmap_",plotname,"_LOCA_hist_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormat_LOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

####

rownames(NRMSEcmatout_LOCA) = colnames(NRMSEcmatout_LOCA) = simnames_LOCA
del.RMSE_LOCA = get_upper_tri(NRMSEcmatout_LOCA)
melted_cormat_LOCA <- melt(del.RMSE_LOCA, na.rm = TRUE)
melted_cormat_LOCA$value = ifelse(melted_cormat_LOCA$value>1,1,melted_cormat_LOCA$value)

rownames(NRMSEcmatout_GCM) = colnames(NRMSEcmatout_GCM) = simnames_GCM
del.RMSE_GCM = get_upper_tri(NRMSEcmatout_GCM)
melted_cormat_GCM <- melt(del.RMSE_GCM, na.rm = TRUE)
melted_cormat_GCM$value = ifelse(melted_cormat_GCM$value>1,1,melted_cormat_GCM$value)

topend = ceiling(max(c(NRMSEcmatout_GCM,NRMSEcmatout_LOCA),na.rm=TRUE))
#topend=1
### Heatmap
pdf(paste("/home/woot0002/DS_ind/NRMSEheatmap_",plotname,"_GCM_change_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormat_GCM, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
                       name="NRMSE") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
dev.off()

pdf(paste("/home/woot0002/DS_ind/NRMSEheatmap_",plotname,"_LOCA_change_",type,"_",statemask,".pdf",sep=""),height=18,width=18)
ggplot(data = melted_cormat_LOCA, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "gray96", 
                       midpoint = topend/2, limit = c(0,topend), breaks=seq(0,topend,by=0.1), space = "Lab", 
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

save(list=c("NRMSEhmatout_LOCA","NRMSEcmatout_LOCA","LOCAhdat","simnames_LOCA"),file=paste("/home/woot0002/DS_ind/RMSEmats_combohc_LOCA_",plotname,"_",type,"_",statemask,".Rdata",sep=""))
save(list=c("NRMSEhmatout_GCM","NRMSEcmatout_GCM","GCMhdat","simnames_GCM"),file=paste("/home/woot0002/DS_ind/RMSEmats_combohc_GCM_",plotname,"_",type,"_",statemask,".Rdata",sep=""))







