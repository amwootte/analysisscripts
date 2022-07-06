
library(optparse)
source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("/data2/3to5/I35/scripts/colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(raster)
library(rasterVis)
library(maptools)
library(ggplot2)
library(zoo)
library(lars)
library(mailR)

load("/home/woot0002/streamflowtests_complexresults_5varsleftout.Rdata")

maxcors = c()
maxcoridx = c()
maxcorRMSE = c()
maxcornames = list()

for(i in 1:14){
  tmp = do.call("c",corswithin)
  maxcors[i] = max(tmp)
  maxcoridx[i] = which(tmp==max(tmp))
  maxcornames[[i]] = varswithin[[i]][[which(tmp==max(tmp))]]
  maxcorRMSE[i]= RMSEswithin[[i]][[which(tmp==max(tmp))]]
}








