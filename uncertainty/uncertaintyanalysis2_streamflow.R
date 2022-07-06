################
#
# Uncertainty analysis gridded

############
# library and function load

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
anfilepath = "analysis"

###
# set extra options

printPDF = TRUE # if you want some plots printed to PDF set this to TRUE

##############
# Read in components from netcdf file

compfile = paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_components_fix2.Rdata",sep="")

# Gather fitted data and residuals for all files
# rather than a list, the fitted values and residuals are ultimately four dimensional arrays
# lon,lat,time,file

load(compfile)

################
# Subset to future period only

timeidx = which(years>=2006)

M=M_smooth[timeidx]
D=D_smooth[timeidx]
S=S_smooth[timeidx]
G=G_smooth[timeidx]
T=T_smooth[timeidx]

################
# Fractional Uncertainty Calcs

FM = FD = FS = FV = FT = M

for(i in 1:length(M)){
FM[i] = (1.65*sqrt(M[i]))/G[i]
FD[i] = (1.65*sqrt(D[i]))/G[i]
FS[i] = (1.65*sqrt(S[i]))/G[i]
FV[i] = (1.65*sqrt(V))/G[i]
FT[i] = (1.65*sqrt(T[i]))/G[i]
}

#################
# Signal to Noise Ratio Calcs

SM = 1/FM
SD = 1/FD
SS = 1/FS
SV = 1/FV
ST = 1/FT

################
# Contribution to the Total Uncertainty per component calc

CM = CD = CS = CV = M

for(i in 1:length(M)){
CM[i] = (M[i]/T[i])*100
CD[i] = (D[i]/T[i])*100
CS[i] = (S[i]/T[i])*100
CV[i] = (V/T[i])*100
}

years = years[timeidx]

###############
# Write out analysis output

fileout = paste("/home/woot0002/uncertainty/",anfilepath,"/CPREP_",varname,"_",datasetname,"_analysis_fix2.Rdata",sep="")

save(list=c("years","FM","FD","FS","FV","FT","SM","SD","SS","SV","ST","CM","CD","CS","CV"),file=fileout)



