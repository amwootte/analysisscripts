source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(mapdata)
library(maptools)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(ggplot2)
library(modi)

########
# Focusing RMSE and Bias dot plots - one variable at a time

varin = "pr"
var = "pr"
group = "GCM"
group2 = "CMIP5"

changevalsfile = paste("/home/woot0002/DS_ind/",var,"_datavals_projchange.Rdata",sep="")
interannvarfile = paste("/home/woot0002/DS_ind/intanntab_",group,"_",varin,".csv",sep="")

load(changevalsfile)

intanntab = read.table(interannvarfile,header=TRUE,sep=",")
intanntabagg = aggregate(intanntab[,4:5],by=list(region=intanntab$region),mean,na.rm=TRUE)
intanntabaggsd = aggregate(intanntab[,4:5],by=list(region=intanntab$region),sd,na.rm=TRUE)

intmodsd = aggregate(memsdat[which(memsdat$DS==group2),6],by=list(region=memsdat$region[which(memsdat$DS==group2)]),sd,na.rm=TRUE)

meansdat = meansdat_all[which(meansdat_all$DS==group2),]

regions = as.character(unique(meansdat$region))

meansdat2 = NULL

for(r in 1:length(regions)){
  
  if(regions[r]=="louisiana"){
    r2 = "LA"
  }
  if(regions[r]=="new mexico"){
    r2 = "NM"
  }
  if(regions[r]=="full"){
    r2 = "full"
  }
  
  tmpdat = meansdat[which(as.character(meansdat$region)==regions[r]),]
  
  idx1 = which(as.character(intanntabagg$region)==r2)
  idx2 = which(as.character(intmodsd$region)==regions[r])
  
  tmpdat$intmodsd = intmodsd$x[idx2]
  tmpdat$histintann = intanntabagg$histintann[idx1]
  tmpdat$projintann = intanntabagg$projintann[idx1]
  
  meansdat2 = rbind(meansdat2,tmpdat)
  
}

####
# rearrange table for plotting

wgs = c("pr","tmax")
plotdat = NULL

for(r in 1:3){
  for(s in 1:3){
    for(w in 1:length(wgs)){
      tmp = meansdat2[which(meansdat2$region==regions[r] & meansdat2$WUdomain==regions[s] & meansdat2$wg==wgs[w]),]
      meandat = c(abs(tmp$changemeans),tmp$intmodsd[1],tmp$histintann[1],tmp$projintann[1])
      group = c(as.character(tmp$group),"intmodsd","histintann","projintann")
      tmpframe = data.frame(group,meandat)
      tmpframe$wg = wgs[w]
      tmpframe$region = regions[r]
      tmpframe$SAdomain = regions[r]
      tmpframe$WUdomain = regions[s]
      tmpframe$var = var
      plotdat = rbind(plotdat,tmpframe)
    }    
  }
}

plotdat$region <- factor(plotdat$region,levels = c('full','new mexico','louisiana'),ordered = TRUE)
plotdat$WUdomain <- factor(plotdat$WUdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)
plotdat$SAdomain <- factor(plotdat$SAdomain,levels = c('full','new mexico','louisiana'),ordered = TRUE)
plotdat$wg <- factor(plotdat$wg,levels = c('tmax','pr'),ordered = TRUE)
plotdat$group <- factor(plotdat$group,levels = c("BMA","SI-c","SI-h","skill","unweighted","intmodsd","histintann","projintann"),ordered = TRUE)


##
# barplots

ggplot(plotdat, aes(x=SAdomain, y=meandat, fill=group)) + 
  geom_bar(position=position_dodge(),width=0.7, stat="identity")+scale_fill_manual(values=c("#f8766d", "#a3a500", "#00bf7d","#00b0f6", "#e76bf3", "#45965c", "#f4b653", "#0520ed"))+
  facet_grid(cols=vars(WUdomain),rows=vars(wg))+geom_hline(yintercept=0,linetype="dashed")+
  ggtitle(paste(group2," Mean Change vs. interann vars - ",var,sep=""))+xlab("region")+ylab("mm")
ggsave(paste("/home/woot0002/DS_ind/interannvars_",group2,"_",varin,".pdf",sep=""),device = "pdf",width=10,height=7)



