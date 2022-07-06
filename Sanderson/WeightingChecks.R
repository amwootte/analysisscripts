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

weighted.var2 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =na.rm)
}

weighted.var3 <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
}

setwd("/home/woot0002/DS_ind/")

var = varin = "pr"
type="ann"
weightingused = "new mexico"
stateapplied = "new mexico"

if(weightingused=="full"){
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,".Rdata",sep=""))
  BMAweights_GCM = read.table(paste("best_BMA_combo_",var,".txt",sep=""))
  BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,".txt",sep=""))
  load(paste("/home/woot0002/DS_ind/BMAposterior_meansandvars_",var,"_WU",weightingused,".Rdata",sep=""))
  BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,".txt",sep=""))
  BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,".txt",sep=""))
} else {
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,"_",weightingused,".Rdata",sep=""))
  BMAweights_GCM = read.table(paste("best_BMA_combo_",var,"_",weightingused,".txt",sep=""))
  BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
  load(paste("/home/woot0002/DS_ind/BMAposterior_meansandvars_",var,"_WU",weightingused,".Rdata",sep=""))
  BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,"_",weightingused,".txt",sep=""))
  BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
}

GCMhdat$BMA = t(BMAweights_GCM)[,1]
LOCAhdat$BMA = t(BMAweights_LOCA)[,1]
GCMhdatpr = GCMhdat
LOCAhdatpr = LOCAhdat
BMAweights_GCMpr = BMAweights_GCM
BMAweights_LOCApr = BMAweights_LOCA
postweights_GCMpr = BMAweightsGCM 
postweights_LOCApr = BMAweightsLOCA

checksGCMpr = c()
diffGCMpr = list()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_GCMpr[j])==as.numeric(postweights_GCMpr[i,j]),1,0)
  }
  checksGCMpr[i] = sum(tmp,na.rm=TRUE)
  diffGCMpr[[i]] = as.numeric(BMAweights_GCMpr)-as.numeric(postweights_GCMpr[i,])
}

checksLOCApr = c()
diffLOCApr = list()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_LOCApr[j])==as.numeric(postweights_LOCApr[i,j]),1,0)
  }
  checksLOCApr[i] = sum(tmp,na.rm=TRUE)
  diffLOCApr[[i]] = as.numeric(BMAweights_LOCApr)-as.numeric(postweights_LOCApr[i,])
}

#####

var = varin = "tmax"

if(weightingused=="full"){
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,".Rdata",sep=""))
  BMAweights_GCM = read.table(paste("best_BMA_combo_",var,".txt",sep=""))
  BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,".txt",sep=""))
  load(paste("/home/woot0002/DS_ind/BMAposterior_meansandvars_",var,"_WU",weightingused,".Rdata",sep=""))
  BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,".txt",sep=""))
  BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,".txt",sep=""))
} else {
  load(file=paste("Sanderson_EnsembleWeights_",var,"_",type,"_",weightingused,".Rdata",sep=""))
  BMAweights_GCM = read.table(paste("best_BMA_combo_",var,"_",weightingused,".txt",sep=""))
  BMAweights_LOCA = read.table(paste("best_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
  load(paste("/home/woot0002/DS_ind/BMAposterior_meansandvars_",var,"_WU",weightingused,".Rdata",sep=""))
  BMAweightsGCM = read.table(paste("posterior_BMA_combo_",var,"_",weightingused,".txt",sep=""))
  BMAweightsLOCA = read.table(paste("posterior_BMA_combo_LOCA_",var,"_",weightingused,".txt",sep=""))
}

GCMhdat$BMA = t(BMAweights_GCM)[,1]
LOCAhdat$BMA = t(BMAweights_LOCA)[,1]
GCMhdattmax = GCMhdat
LOCAhdattmax = LOCAhdat
BMAweights_GCMtmax = BMAweights_GCM
BMAweights_LOCAtmax = BMAweights_LOCA
postweights_GCMtmax = BMAweightsGCM 
postweights_LOCAtmax = BMAweightsLOCA


checksGCMtmax = c()
diffGCMtmax = list()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_GCMtmax[j])==as.numeric(postweights_GCMtmax[i,j]),1,0)
  }
  checksGCMtmax[i] = sum(tmp,na.rm=TRUE)
  diffGCMtmax[[i]] = as.numeric(BMAweights_GCMtmax)-as.numeric(postweights_GCMtmax[i,])
}

checksLOCAtmax = c()
diffLOCAtmax = list()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_LOCAtmax[j])==as.numeric(postweights_LOCAtmax[i,j]),1,0)
  }
  checksLOCAtmax[i] = sum(tmp,na.rm=TRUE)
  diffLOCAtmax[[i]] = as.numeric(BMAweights_LOCAtmax)-as.numeric(postweights_LOCAtmax[i,])
}

################
# Does one set of weights exist in the other? BMA best weights identical between pr vs. tmax weighting?

BMAweights_GCMtmax==BMAweights_GCMpr
BMAweights_LOCAtmax==BMAweights_LOCApr

################
# Does one set of weights exist in the other? BMA posterior weights identical between pr and tmax weighting?

any(postweights_GCMtmax == postweights_GCMpr)
any(postweights_LOCAtmax == postweights_LOCApr)

################
# Does one set of weights exist in the other? BMA best weights of tmax or pr exist in posteriors of pr or tmax?


checksGCMpr_in_tmax = c()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_GCMpr[j])==as.numeric(postweights_GCMtmax[i,j]),1,0)
  }
  checksGCMpr_in_tmax[i] = sum(tmp,na.rm=TRUE)
}

checksLOCApr_in_tmax = c()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_LOCApr[j])==as.numeric(postweights_LOCAtmax[i,j]),1,0)
  }
  checksLOCApr_in_tmax[i] = sum(tmp,na.rm=TRUE)
}

checksGCMtmax_in_pr = c()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_GCMtmax[j])==as.numeric(postweights_GCMpr[i,j]),1,0)
  }
  checksGCMtmax_in_pr[i] = sum(tmp,na.rm=TRUE)
}

checksLOCAtmax_in_pr = c()
for(i in 1:100){
  tmp=c()
  for(j in 1:26){
    tmp[j] = ifelse(as.numeric(BMAweights_LOCApr[j])==as.numeric(postweights_LOCAtmax[i,j]),1,0)
  }
  checksLOCAtmax_in_pr[i] = sum(tmp,na.rm=TRUE)
}



