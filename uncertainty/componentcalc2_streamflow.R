################
#
# Gridded component calculator

############
# library and function load

###########
# Just in case random things got saved, clean the workspace first

rm(list=ls(all=TRUE))

############
# load libraries

library(ncdf4)
library(fields)
library(maps)

############
# load calculator functions

source("/home/woot0002/scripts/uncertainty/componentfunctions_2021.R")

###
# set variable name

varname = "TAVG_12M"
datasetname = "08144500"
regfilepath = "regfits"
anfilepath = "analysis"

############
# Load fitted values and residuals from all GCMs

# determine files
load(paste("/home/woot0002/uncertainty/",regfilepath,"/CPREP_",varname,"_",datasetname,"_regfit.Rdata",sep=""))
load(paste("/home/woot0002/uncertainty/",regfilepath,"/",varname,"_",datasetname,"_regfit.Rdata",sep=""))


############
# Load GCM and DS weights

load(paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_weights.Rdata",sep=""))

#############
# Get unique GCM, DS, and scenario names from filenames

memnames = names(fittedvals)[2:ncol(fittedvals)]
filesplit = do.call(rbind,strsplit(names(fittedvals)[2:ncol(fittedvals)],"_",fixed=TRUE))
scens = unique(filesplit[,4])

##############
# Calculate V - Internal Variability

V = intvar_spts(resvals,weights=list(GCMweights,DSweights),GCMs,DSs,files=memnames)

##############
# Calculate M - GCM Uncertainty

M = GCMuncts(fittedvals,weights=GCMweights,scens,DSs,GCMs,files=memnames)

##############
# Calculate D - DS Uncertainty

D = DSuncts(fittedvals,weights=DSweights,scens,DSs,GCMs,files=memnames)

##############
# Calculate S - Scenario Uncertainty

S = scenuncts(fittedvals,scens,DSs,GCMs,weightsGCM=GCMweights,weightsDS=DSweights,files=memnames)

##############
# Smoothed values of M, D, and S

M_smooth = M #aperm(apply(M,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
D_smooth = D #aperm(apply(D,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))
S_smooth = S #aperm(apply(S,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

#####
# map check

#testsfc1=list(x=lon,y=lat,z=V)
#testsfc2=list(x=lon,y=lat,z=M_smooth[,,100])
#testsfc3=list(x=lon,y=lat,z=D_smooth[,,100])
#testsfc4=list(x=lon,y=lat,z=S_smooth[,,100])

#par(mfrow=c(2,2))
#surface(testsfc1,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc2,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc3,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)
#surface(testsfc4,type="I",zlim=c(0,0.25))
#map("state",add=TRUE)

###############
# Calculate mean change

G=meanchangets(fittedvals,scens,GCMs,DSs,GCMweights,DSweights,files=memnames)
G_smooth = G #aperm(apply(G,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

################
# Calculate Total Variance

T = c()
for(t in 1:length(M)) T[t] = M[t]+D[t]+S[t]+V
T_smooth = T #aperm(apply(T,c(1,2),filter,sides=2,filter=rep(1,10)/10),c(2,3,1))

years = fittedvals[,1]

###############
# Save out files
fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_components_fix2.Rdata",sep="")

save(list=c("years","V","M_smooth","D_smooth","S_smooth","T_smooth","G_smooth","GCMs","DSs"),file=fileout)

