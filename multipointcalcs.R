##################
#
# Multiple Location Point calculation
setwd("/data2/3to5/I35/scripts/")
location = "Guymon"
lon = 258.5185
lat = 36.6828
period = "2041-2070"

vars = c("tasmax","tasmin","tmax95","tmax100","tmin32","tmin28","heatwaves","gsl","frd","mdrn","pr25","pr50","rx1day","rx5day","cwd","cdd","pr")
type = c(rep("absolute",16),"percent")

commandstart = "Rscript /data2/3to5/I35/scripts/point_calcs.R -i /data2/3to5/I35/all_mems/"
commandmiddle = " -s CCSM4_DeltaSD_Daymet,CCSM4_DeltaSD_Livneh,CCSM4_DeltaSD_PRISM,MIROC5_DeltaSD_Daymet,MIROC5_DeltaSD_Livneh,MIROC5_DeltaSD_PRISM,MPI-ESM-LR_DeltaSD_Daymet,MPI-ESM-LR_DeltaSD_Livneh,MPI-ESM-LR_DeltaSD_PRISM,CCSM4_QDM_Daymet,CCSM4_QDM_Livneh,CCSM4_QDM_PRISM,MIROC5_QDM_Daymet,MIROC5_QDM_Livneh,MIROC5_QDM_PRISM,MPI-ESM-LR_QDM_Daymet,MPI-ESM-LR_QDM_Livneh,MPI-ESM-LR_QDM_PRISM -p CCSM4_DeltaSD_Daymet_rcp26,CCSM4_DeltaSD_Livneh_rcp26,CCSM4_DeltaSD_PRISM_rcp26,CCSM4_DeltaSD_Daymet_rcp45,CCSM4_DeltaSD_Livneh_rcp45,CCSM4_DeltaSD_PRISM_rcp45,CCSM4_DeltaSD_Daymet_rcp85,CCSM4_DeltaSD_Livneh_rcp85,CCSM4_DeltaSD_PRISM_rcp85,MIROC5_DeltaSD_Daymet_rcp26,MIROC5_DeltaSD_Livneh_rcp26,MIROC5_DeltaSD_PRISM_rcp26,MIROC5_DeltaSD_Daymet_rcp45,MIROC5_DeltaSD_Livneh_rcp45,MIROC5_DeltaSD_PRISM_rcp45,MIROC5_DeltaSD_Daymet_rcp85,MIROC5_DeltaSD_Livneh_rcp85,MIROC5_DeltaSD_PRISM_rcp85,MPI-ESM-LR_DeltaSD_Daymet_rcp26,MPI-ESM-LR_DeltaSD_Livneh_rcp26,MPI-ESM-LR_DeltaSD_PRISM_rcp26,MPI-ESM-LR_DeltaSD_Daymet_rcp45,MPI-ESM-LR_DeltaSD_Livneh_rcp45,MPI-ESM-LR_DeltaSD_PRISM_rcp45,MPI-ESM-LR_DeltaSD_Daymet_rcp85,MPI-ESM-LR_DeltaSD_Livneh_rcp85,MPI-ESM-LR_DeltaSD_PRISM_rcp85,CCSM4_QDM_Daymet_rcp26,CCSM4_QDM_Livneh_rcp26,CCSM4_QDM_PRISM_rcp26,CCSM4_QDM_Daymet_rcp45,CCSM4_QDM_Livneh_rcp45,CCSM4_QDM_PRISM_rcp45,CCSM4_QDM_Daymet_rcp85,CCSM4_QDM_Livneh_rcp85,CCSM4_QDM_PRISM_rcp85,MIROC5_QDM_Daymet_rcp26,MIROC5_QDM_Livneh_rcp26,MIROC5_QDM_PRISM_rcp26,MIROC5_QDM_Daymet_rcp45,MIROC5_QDM_Livneh_rcp45,MIROC5_QDM_PRISM_rcp45,MIROC5_QDM_Daymet_rcp85,MIROC5_QDM_Livneh_rcp85,MIROC5_QDM_PRISM_rcp85,MPI-ESM-LR_QDM_Daymet_rcp26,MPI-ESM-LR_QDM_Livneh_rcp26,MPI-ESM-LR_QDM_PRISM_rcp26,MPI-ESM-LR_QDM_Daymet_rcp45,MPI-ESM-LR_QDM_Livneh_rcp45,MPI-ESM-LR_QDM_PRISM_rcp45,MPI-ESM-LR_QDM_Daymet_rcp85,MPI-ESM-LR_QDM_Livneh_rcp85,MPI-ESM-LR_QDM_PRISM_rcp85 -n "

for(i in 1:length(vars)){
  #"Rscript point_calcs.R -i /data2/3to5/I35/all_mems/tmin28_allmem_absolute_2041-2070_ann.nc -s CCSM4_DeltaSD_Daymet,CCSM4_DeltaSD_Livneh,CCSM4_DeltaSD_PRISM,MIROC5_DeltaSD_Daymet,MIROC5_DeltaSD_Livneh,MIROC5_DeltaSD_PRISM,MPI-ESM-LR_DeltaSD_Daymet,MPI-ESM-LR_DeltaSD_Livneh,MPI-ESM-LR_DeltaSD_PRISM,CCSM4_QDM_Daymet,CCSM4_QDM_Livneh,CCSM4_QDM_PRISM,MIROC5_QDM_Daymet,MIROC5_QDM_Livneh,MIROC5_QDM_PRISM,MPI-ESM-LR_QDM_Daymet,MPI-ESM-LR_QDM_Livneh,MPI-ESM-LR_QDM_PRISM -p CCSM4_DeltaSD_Daymet_rcp26,CCSM4_DeltaSD_Livneh_rcp26,CCSM4_DeltaSD_PRISM_rcp26,CCSM4_DeltaSD_Daymet_rcp45,CCSM4_DeltaSD_Livneh_rcp45,CCSM4_DeltaSD_PRISM_rcp45,CCSM4_DeltaSD_Daymet_rcp85,CCSM4_DeltaSD_Livneh_rcp85,CCSM4_DeltaSD_PRISM_rcp85,MIROC5_DeltaSD_Daymet_rcp26,MIROC5_DeltaSD_Livneh_rcp26,MIROC5_DeltaSD_PRISM_rcp26,MIROC5_DeltaSD_Daymet_rcp45,MIROC5_DeltaSD_Livneh_rcp45,MIROC5_DeltaSD_PRISM_rcp45,MIROC5_DeltaSD_Daymet_rcp85,MIROC5_DeltaSD_Livneh_rcp85,MIROC5_DeltaSD_PRISM_rcp85,MPI-ESM-LR_DeltaSD_Daymet_rcp26,MPI-ESM-LR_DeltaSD_Livneh_rcp26,MPI-ESM-LR_DeltaSD_PRISM_rcp26,MPI-ESM-LR_DeltaSD_Daymet_rcp45,MPI-ESM-LR_DeltaSD_Livneh_rcp45,MPI-ESM-LR_DeltaSD_PRISM_rcp45,MPI-ESM-LR_DeltaSD_Daymet_rcp85,MPI-ESM-LR_DeltaSD_Livneh_rcp85,MPI-ESM-LR_DeltaSD_PRISM_rcp85,CCSM4_QDM_Daymet_rcp26,CCSM4_QDM_Livneh_rcp26,CCSM4_QDM_PRISM_rcp26,CCSM4_QDM_Daymet_rcp45,CCSM4_QDM_Livneh_rcp45,CCSM4_QDM_PRISM_rcp45,CCSM4_QDM_Daymet_rcp85,CCSM4_QDM_Livneh_rcp85,CCSM4_QDM_PRISM_rcp85,MIROC5_QDM_Daymet_rcp26,MIROC5_QDM_Livneh_rcp26,MIROC5_QDM_PRISM_rcp26,MIROC5_QDM_Daymet_rcp45,MIROC5_QDM_Livneh_rcp45,MIROC5_QDM_PRISM_rcp45,MIROC5_QDM_Daymet_rcp85,MIROC5_QDM_Livneh_rcp85,MIROC5_QDM_PRISM_rcp85,MPI-ESM-LR_QDM_Daymet_rcp26,MPI-ESM-LR_QDM_Livneh_rcp26,MPI-ESM-LR_QDM_PRISM_rcp26,MPI-ESM-LR_QDM_Daymet_rcp45,MPI-ESM-LR_QDM_Livneh_rcp45,MPI-ESM-LR_QDM_PRISM_rcp45,MPI-ESM-LR_QDM_Daymet_rcp85,MPI-ESM-LR_QDM_Livneh_rcp85,MPI-ESM-LR_QDM_PRISM_rcp85 -n Anadarko -x 261.7563 -y 35.0726"
  command = paste(commandstart,vars[i],"_allmem_",type[i],"_",period,"_ann.nc",commandmiddle,location," -x ",lon," -y ",lat,sep="")
  system(command,wait=TRUE)
  
  message("Finished calcs for var ",vars[i])
}


