
source("/data2/3to5/I35/scripts/analysisfunctions.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
library(PCICt)

###
# historical GCM files

GCM_hist_tasmax = system("ls /home/woot0002/GCMs/regrid/tasmax*histclimo*",intern=TRUE)
GCM_hist_tasmin = system("ls /home/woot0002/GCMs/regrid/tasmin*histclimo*",intern=TRUE)
GCM_hist_pr = system("ls /home/woot0002/GCMs/regrid/pr*histclimo*",intern=TRUE)

###
# future GCM files

GCM_proj_tasmax = system("ls /home/woot0002/GCMs/regrid/tasmax*projclimo*",intern=TRUE)
GCM_proj_tasmin = system("ls /home/woot0002/GCMs/regrid/tasmin*projclimo*",intern=TRUE)
GCM_proj_pr = system("ls /home/woot0002/GCMs/regrid/pr*projclimo*",intern=TRUE)

###
# historical LOCA files

LOCA_hist_tasmax = system("ls /home/woot0002/LOCA/regrid/tasmax*histclimo*",intern=TRUE)
LOCA_hist_tasmin = system("ls /home/woot0002/LOCA/regrid/tasmin*histclimo*",intern=TRUE)
LOCA_hist_pr = system("ls /home/woot0002/LOCA/regrid/pr*histclimo*",intern=TRUE)

###
# future LOCA files

LOCA_proj_tasmax = system("ls /home/woot0002/LOCA/regrid/tasmax*projclimo*",intern=TRUE)
LOCA_proj_tasmin = system("ls /home/woot0002/LOCA/regrid/tasmin*projclimo*",intern=TRUE)
LOCA_proj_pr = system("ls /home/woot0002/LOCA/regrid/pr*projclimo*",intern=TRUE)

######

filebreakdown = do.call(rbind,strsplit(GCM_hist_pr,"_",fixed=TRUE))
GCMs_pr = filebreakdown[,2]

filebreakdown = do.call(rbind,strsplit(GCM_hist_tasmax,"_",fixed=TRUE))
GCMs_tasmax = filebreakdown[,2]

filebreakdown = do.call(rbind,strsplit(GCM_hist_tasmin,"_",fixed=TRUE))
GCMs_tasmin = filebreakdown[,2]

GCMlist = GCMs_pr

GCM_proj_pr_check = GCM_proj_tasmax_check = GCM_proj_tasmin_check = c()
LOCA_proj_pr_check = LOCA_proj_tasmax_check = LOCA_proj_tasmin_check = c()
LOCA_hist_pr_check = LOCA_hist_tasmax_check = LOCA_hist_tasmin_check = c()

for(i in 1:length(GCMlist)){
  
  GCM_proj_pr_check[i] = length(grep(GCMlist[i],GCM_proj_pr))
  GCM_proj_tasmax_check[i] = length(grep(GCMlist[i],GCM_proj_tasmax))
  GCM_proj_tasmin_check[i] = length(grep(GCMlist[i],GCM_proj_tasmin))
  
  LOCA_proj_pr_check[i] = length(grep(GCMlist[i],LOCA_proj_pr))
  LOCA_proj_tasmax_check[i] = length(grep(GCMlist[i],LOCA_proj_tasmax))
  LOCA_proj_tasmin_check[i] = length(grep(GCMlist[i],LOCA_proj_tasmin))
  
  LOCA_hist_pr_check[i] = length(grep(GCMlist[i],LOCA_hist_pr))
  LOCA_hist_tasmax_check[i] = length(grep(GCMlist[i],LOCA_hist_tasmax))
  LOCA_hist_tasmin_check[i] = length(grep(GCMlist[i],LOCA_hist_tasmin))
  
}

save(list=c("GCMlist"),file="/home/woot0002/GCMlist.Rdata")




