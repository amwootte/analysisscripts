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
# grid basics for later

GCMfiles_pr = system("ls /home/woot0002/GCMs/regrid/pr_*histclimo*.nc",intern=TRUE)
  nctest = nc_open(GCMfiles_pr[1])
  lat = ncvar_get(nctest,"lat")
  lon=ncvar_get(nctest,"lon")
  nc_close(nctest)

########

setwd("/home/woot0002/DS_ind/")
load("WeightedMeansVars_pr_pr.Rdata")
load("WeightedMeansVars_tmax_pr.Rdata")
meansdat_pr_Wpr = meansdat_pr
meansdat_tmax_Wpr = meansdat_tmax
meansdat_pr_Wpr$wg ="pr"
meansdat_pr_Wpr$var ="pr"
meansdat_tmax_Wpr$wg="pr"
meansdat_tmax_Wpr$var="tmax"

#####

load("Biasesbasedon_pr_weighting.Rdata")

GCMunweightedmean_bias_pr_Wpr=GCMunweightedmean_bias_pr
GCMunweightedmean_bias_tmax_Wpr=GCMunweightedmean_bias_tmax
LOCAunweightedmean_bias_pr_Wpr=LOCAunweightedmean_bias_pr
LOCAunweightedmean_bias_tmax_Wpr=LOCAunweightedmean_bias_tmax

GCMskillmean_bias_pr_Wpr=GCMskillmean_bias_pr
LOCAskillmean_bias_pr_Wpr=LOCAskillmean_bias_pr
GCMskillmean_bias_tmax_Wpr=GCMskillmean_bias_tmax
LOCAskillmean_bias_tmax_Wpr=LOCAskillmean_bias_tmax

GCMSIhmean_bias_pr_Wpr=GCMSIhmean_bias_pr
GCMSIhmean_bias_tmax_Wpr=GCMSIhmean_bias_tmax
LOCASIhmean_bias_pr_Wpr=LOCASIhmean_bias_pr
LOCASIhmean_bias_tmax_Wpr=LOCASIhmean_bias_tmax

GCMSIcmean_bias_pr_Wpr=GCMSIcmean_bias_pr
GCMSIcmean_bias_tmax_Wpr=GCMSIcmean_bias_tmax
LOCASIcmean_bias_pr_Wpr=LOCASIcmean_bias_pr
LOCASIcmean_bias_tmax_Wpr=LOCASIcmean_bias_tmax

#####

load("Changesbasedon_pr_weighting.Rdata")

GCMunweightedmean_change_pr_Wpr=GCMunweightedmean_change_pr
GCMunweightedmean_change_tmax_Wpr=GCMunweightedmean_change_tmax
LOCAunweightedmean_change_pr_Wpr=LOCAunweightedmean_change_pr
LOCAunweightedmean_change_tmax_Wpr=LOCAunweightedmean_change_tmax

GCMskillmean_change_pr_Wpr=GCMskillmean_change_pr
LOCAskillmean_change_pr_Wpr=LOCAskillmean_change_pr
GCMskillmean_change_tmax_Wpr=GCMskillmean_change_tmax
LOCAskillmean_change_tmax_Wpr=LOCAskillmean_change_tmax

GCMSIhmean_change_pr_Wpr=GCMSIhmean_change_pr
GCMSIhmean_change_tmax_Wpr=GCMSIhmean_change_tmax
LOCASIhmean_change_pr_Wpr=LOCASIhmean_change_pr
LOCASIhmean_change_tmax_Wpr=LOCASIhmean_change_tmax

GCMSIcmean_change_pr_Wpr=GCMSIcmean_change_pr
GCMSIcmean_change_tmax_Wpr=GCMSIcmean_change_tmax
LOCASIcmean_change_pr_Wpr=LOCASIcmean_change_pr
LOCASIcmean_change_tmax_Wpr=LOCASIcmean_change_tmax

#####

load("ChangeVarsbasedon_pr_weighting.Rdata")

GCMunweightedvar_change_pr_Wpr=GCMunweightedvar_change_pr
GCMunweightedvar_change_tmax_Wpr=GCMunweightedvar_change_tmax
LOCAunweightedvar_change_pr_Wpr=LOCAunweightedvar_change_pr
LOCAunweightedvar_change_tmax_Wpr=LOCAunweightedvar_change_tmax

GCMskillvar_change_pr_Wpr=GCMskillvar_change_pr
LOCAskillvar_change_pr_Wpr=LOCAskillvar_change_pr
GCMskillvar_change_tmax_Wpr=GCMskillvar_change_tmax
LOCAskillvar_change_tmax_Wpr=LOCAskillvar_change_tmax

GCMSIhvar_change_pr_Wpr=GCMSIhvar_change_pr
GCMSIhvar_change_tmax_Wpr=GCMSIhvar_change_tmax
LOCASIhvar_change_pr_Wpr=LOCASIhvar_change_pr
LOCASIhvar_change_tmax_Wpr=LOCASIhvar_change_tmax

GCMSIcvar_change_pr_Wpr=GCMSIcvar_change_pr
GCMSIcvar_change_tmax_Wpr=GCMSIcvar_change_tmax
LOCASIcvar_change_pr_Wpr=LOCASIcvar_change_pr
LOCASIcvar_change_tmax_Wpr=LOCASIcvar_change_tmax


######
load("WeightedMeansVars_pr_tmax.Rdata")
load("WeightedMeansVars_tmax_tmax.Rdata")
meansdat_pr_Wtx = meansdat_pr
meansdat_tmax_Wtx = meansdat_tmax
meansdat_pr_Wtx$wg ="tmax"
meansdat_pr_Wtx$var ="pr"
meansdat_tmax_Wtx$wg="tmax"
meansdat_tmax_Wtx$var="tmax"

#####

load("Biasesbasedon_tmax_weighting.Rdata")

GCMunweightedmean_bias_pr_Wtx=GCMunweightedmean_bias_pr
GCMunweightedmean_bias_tmax_Wtx=GCMunweightedmean_bias_tmax
LOCAunweightedmean_bias_pr_Wtx=LOCAunweightedmean_bias_pr
LOCAunweightedmean_bias_tmax_Wtx=LOCAunweightedmean_bias_tmax

GCMskillmean_bias_pr_Wtx=GCMskillmean_bias_pr
LOCAskillmean_bias_pr_Wtx=LOCAskillmean_bias_pr
GCMskillmean_bias_tmax_Wtx=GCMskillmean_bias_tmax
LOCAskillmean_bias_tmax_Wtx=LOCAskillmean_bias_tmax

GCMSIhmean_bias_pr_Wtx=GCMSIhmean_bias_pr
GCMSIhmean_bias_tmax_Wtx=GCMSIhmean_bias_tmax
LOCASIhmean_bias_pr_Wtx=LOCASIhmean_bias_pr
LOCASIhmean_bias_tmax_Wtx=LOCASIhmean_bias_tmax

GCMSIcmean_bias_pr_Wtx=GCMSIcmean_bias_pr
GCMSIcmean_bias_tmax_Wtx=GCMSIcmean_bias_tmax
LOCASIcmean_bias_pr_Wtx=LOCASIcmean_bias_pr
LOCASIcmean_bias_tmax_Wtx=LOCASIcmean_bias_tmax

######
load("Changesbasedon_tmax_weighting.Rdata")

GCMunweightedmean_change_pr_Wtx=GCMunweightedmean_change_pr
GCMunweightedmean_change_tmax_Wtx=GCMunweightedmean_change_tmax
LOCAunweightedmean_change_pr_Wtx=LOCAunweightedmean_change_pr
LOCAunweightedmean_change_tmax_Wtx=LOCAunweightedmean_change_tmax

GCMskillmean_change_pr_Wtx=GCMskillmean_change_pr
LOCAskillmean_change_pr_Wtx=LOCAskillmean_change_pr
GCMskillmean_change_tmax_Wtx=GCMskillmean_change_tmax
LOCAskillmean_change_tmax_Wtx=LOCAskillmean_change_tmax

GCMSIhmean_change_pr_Wtx=GCMSIhmean_change_pr
GCMSIhmean_change_tmax_Wtx=GCMSIhmean_change_tmax
LOCASIhmean_change_pr_Wtx=LOCASIhmean_change_pr
LOCASIhmean_change_tmax_Wtx=LOCASIhmean_change_tmax

GCMSIcmean_change_pr_Wtx=GCMSIcmean_change_pr
GCMSIcmean_change_tmax_Wtx=GCMSIcmean_change_tmax
LOCASIcmean_change_pr_Wtx=LOCASIcmean_change_pr
LOCASIcmean_change_tmax_Wtx=LOCASIcmean_change_tmax

######
load("ChangeVarsbasedon_tmax_weighting.Rdata")

GCMunweightedvar_change_pr_Wtx=GCMunweightedvar_change_pr
GCMunweightedvar_change_tmax_Wtx=GCMunweightedvar_change_tmax
LOCAunweightedvar_change_pr_Wtx=LOCAunweightedvar_change_pr
LOCAunweightedvar_change_tmax_Wtx=LOCAunweightedvar_change_tmax

GCMskillvar_change_pr_Wtx=GCMskillvar_change_pr
LOCAskillvar_change_pr_Wtx=LOCAskillvar_change_pr
GCMskillvar_change_tmax_Wtx=GCMskillvar_change_tmax
LOCAskillvar_change_tmax_Wtx=LOCAskillvar_change_tmax

GCMSIhvar_change_pr_Wtx=GCMSIhvar_change_pr
GCMSIhvar_change_tmax_Wtx=GCMSIhvar_change_tmax
LOCASIhvar_change_pr_Wtx=LOCASIhvar_change_pr
LOCASIhvar_change_tmax_Wtx=LOCASIhvar_change_tmax

GCMSIcvar_change_pr_Wtx=GCMSIcvar_change_pr
GCMSIcvar_change_tmax_Wtx=GCMSIcvar_change_tmax
LOCASIcvar_change_pr_Wtx=LOCASIcvar_change_pr
LOCASIcvar_change_tmax_Wtx=LOCASIcvar_change_tmax

######
load("WeightedMeansVars_pr_mv.Rdata")
load("WeightedMeansVars_tmax_mv.Rdata")
meansdat_pr_Wmv = meansdat_pr
meansdat_tmax_Wmv = meansdat_tmax
meansdat_pr_Wmv$wg ="mv"
meansdat_pr_Wmv$var ="pr"
meansdat_tmax_Wmv$wg="mv"
meansdat_tmax_Wmv$var="tmax"

#####
load("Biasesbasedon_mv_weighting.Rdata")

GCMunweightedmean_bias_pr_Wmv=GCMunweightedmean_bias_pr
GCMunweightedmean_bias_tmax_Wmv=GCMunweightedmean_bias_tmax
LOCAunweightedmean_bias_pr_Wmv=LOCAunweightedmean_bias_pr
LOCAunweightedmean_bias_tmax_Wmv=LOCAunweightedmean_bias_tmax

GCMskillmean_bias_pr_Wmv=GCMskillmean_bias_pr
LOCAskillmean_bias_pr_Wmv=LOCAskillmean_bias_pr
GCMskillmean_bias_tmax_Wmv=GCMskillmean_bias_tmax
LOCAskillmean_bias_tmax_Wmv=LOCAskillmean_bias_tmax

GCMSIhmean_bias_pr_Wmv=GCMSIhmean_bias_pr
GCMSIhmean_bias_tmax_Wmv=GCMSIhmean_bias_tmax
LOCASIhmean_bias_pr_Wmv=LOCASIhmean_bias_pr
LOCASIhmean_bias_tmax_Wmv=LOCASIhmean_bias_tmax

GCMSIcmean_bias_pr_Wmv=GCMSIcmean_bias_pr
GCMSIcmean_bias_tmax_Wmv=GCMSIcmean_bias_tmax
LOCASIcmean_bias_pr_Wmv=LOCASIcmean_bias_pr
LOCASIcmean_bias_tmax_Wmv=LOCASIcmean_bias_tmax

######

load("Changesbasedon_mv_weighting.Rdata")

GCMunweightedmean_change_pr_Wmv=GCMunweightedmean_change_pr
GCMunweightedmean_change_tmax_Wmv=GCMunweightedmean_change_tmax
LOCAunweightedmean_change_pr_Wmv=LOCAunweightedmean_change_pr
LOCAunweightedmean_change_tmax_Wmv=LOCAunweightedmean_change_tmax

GCMskillmean_change_pr_Wmv=GCMskillmean_change_pr
LOCAskillmean_change_pr_Wmv=LOCAskillmean_change_pr
GCMskillmean_change_tmax_Wmv=GCMskillmean_change_tmax
LOCAskillmean_change_tmax_Wmv=LOCAskillmean_change_tmax

GCMSIhmean_change_pr_Wmv=GCMSIhmean_change_pr
GCMSIhmean_change_tmax_Wmv=GCMSIhmean_change_tmax
LOCASIhmean_change_pr_Wmv=LOCASIhmean_change_pr
LOCASIhmean_change_tmax_Wmv=LOCASIhmean_change_tmax

GCMSIcmean_change_pr_Wmv=GCMSIcmean_change_pr
GCMSIcmean_change_tmax_Wmv=GCMSIcmean_change_tmax
LOCASIcmean_change_pr_Wmv=LOCASIcmean_change_pr
LOCASIcmean_change_tmax_Wmv=LOCASIcmean_change_tmax

######
load("ChangeVarsbasedon_mv_weighting.Rdata")

GCMunweightedvar_change_pr_Wmv=GCMunweightedvar_change_pr
GCMunweightedvar_change_tmax_Wmv=GCMunweightedvar_change_tmax
LOCAunweightedvar_change_pr_Wmv=LOCAunweightedvar_change_pr
LOCAunweightedvar_change_tmax_Wmv=LOCAunweightedvar_change_tmax

GCMskillvar_change_pr_Wmv=GCMskillvar_change_pr
LOCAskillvar_change_pr_Wmv=LOCAskillvar_change_pr
GCMskillvar_change_tmax_Wmv=GCMskillvar_change_tmax
LOCAskillvar_change_tmax_Wmv=LOCAskillvar_change_tmax

GCMSIhvar_change_pr_Wmv=GCMSIhvar_change_pr
GCMSIhvar_change_tmax_Wmv=GCMSIhvar_change_tmax
LOCASIhvar_change_pr_Wmv=LOCASIhvar_change_pr
LOCASIhvar_change_tmax_Wmv=LOCASIhvar_change_tmax

GCMSIcvar_change_pr_Wmv=GCMSIcvar_change_pr
GCMSIcvar_change_tmax_Wmv=GCMSIcvar_change_tmax
LOCASIcvar_change_pr_Wmv=LOCASIcvar_change_pr
LOCASIcvar_change_tmax_Wmv=LOCASIcvar_change_tmax
#####

meansdat_pr = rbind(meansdat_pr_Wpr,meansdat_pr_Wtx)
meansdat_pr = rbind(meansdat_pr,meansdat_pr_Wmv)
meansdat_tmax = rbind(meansdat_tmax_Wpr,meansdat_tmax_Wtx)
meansdat_tmax = rbind(meansdat_tmax,meansdat_tmax_Wmv)


######
meansdat_pr$region <- factor(meansdat_pr$region,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)
meansdat_pr$wg <- factor(meansdat_pr$wg,levels = c('tmax','pr','mv'),ordered = TRUE)
meansdat_tmax$region <- factor(meansdat_tmax$region,levels = c('full','new mexico','texas','oklahoma','louisiana'),ordered = TRUE)
meansdat_tmax$wg <- factor(meansdat_tmax$wg,levels = c('tmax','pr','mv'),ordered = TRUE)
ggplot(meansdat_pr, aes(x=region, y=bias))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("Bias of Ens. Means against Livneh")+xlab("Region")+ylab("Bias (mm)")+facet_grid(cols=vars(wg))
ggplot(meansdat_pr, aes(x=region, y=rmse)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region")+ylab("RMSE (mm)")+facet_grid(cols=vars(wg))
ggplot(meansdat_pr, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Change (mm)")+facet_grid(cols=vars(wg))
ggplot(meansdat_pr, aes(x=region, y=changevars)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Variance Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Var (mm^2)")+facet_grid(cols=vars(wg))

ggplot(meansdat_tmax, aes(x=region, y=bias))+geom_point(aes(colour = factor(group),shape=factor(DS)),size=5)+geom_hline(yintercept=0,linetype="dashed")+ggtitle("Bias of Ens. Means against Livneh")+xlab("Region")+ylab("Bias (deg C)")+facet_grid(cols=vars(wg))
ggplot(meansdat_tmax, aes(x=region, y=rmse)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("RMSE of Ens. Means against Livneh")+xlab("Region")+ylab("RMSE (deg C)")+facet_grid(cols=vars(wg))
ggplot(meansdat_tmax, aes(x=region, y=changemeans)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Mean Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Change (deg C)")+facet_grid(cols=vars(wg))
ggplot(meansdat_tmax, aes(x=region, y=changevars)) + geom_point(aes(colour = factor(group),shape=factor(DS)),size=5) +geom_hline(yintercept=0,linetype="dashed") + ggtitle("Ens. Variance Projected Changes - RCP8.5 2070-2099")+xlab("Region")+ylab("Var (degC^2)")+facet_grid(cols=vars(wg))

####
# maps

biasvec_pr = c(LOCAunweightedmean_bias_pr_Wpr,LOCAskillmean_bias_pr_Wpr,LOCASIhmean_bias_pr_Wpr,LOCASIcmean_bias_pr_Wpr,
            GCMunweightedmean_bias_pr_Wpr,GCMskillmean_bias_pr_Wpr,GCMSIhmean_bias_pr_Wpr,GCMSIcmean_bias_pr_Wpr,
            LOCAunweightedmean_bias_pr_Wtx,LOCAskillmean_bias_pr_Wtx,LOCASIhmean_bias_pr_Wtx,LOCASIcmean_bias_pr_Wtx,
            GCMunweightedmean_bias_pr_Wtx,GCMskillmean_bias_pr_Wtx,GCMSIhmean_bias_pr_Wtx,GCMSIcmean_bias_pr_Wtx,
            LOCAunweightedmean_bias_pr_Wmv,LOCAskillmean_bias_pr_Wmv,LOCASIhmean_bias_pr_Wmv,LOCASIcmean_bias_pr_Wmv,
            GCMunweightedmean_bias_pr_Wmv,GCMskillmean_bias_pr_Wmv,GCMSIhmean_bias_pr_Wmv,GCMSIcmean_bias_pr_Wmv)

biasvec_tmax = c(LOCAunweightedmean_bias_tmax_Wpr,LOCAskillmean_bias_tmax_Wpr,LOCASIhmean_bias_tmax_Wpr,LOCASIcmean_bias_tmax_Wpr,
               GCMunweightedmean_bias_tmax_Wpr,GCMskillmean_bias_tmax_Wpr,GCMSIhmean_bias_tmax_Wpr,GCMSIcmean_bias_tmax_Wpr,
               LOCAunweightedmean_bias_tmax_Wtx,LOCAskillmean_bias_tmax_Wtx,LOCASIhmean_bias_tmax_Wtx,LOCASIcmean_bias_tmax_Wtx,
               GCMunweightedmean_bias_tmax_Wtx,GCMskillmean_bias_tmax_Wtx,GCMSIhmean_bias_tmax_Wtx,GCMSIcmean_bias_tmax_Wtx,
               LOCAunweightedmean_bias_tmax_Wmv,LOCAskillmean_bias_tmax_Wmv,LOCASIhmean_bias_tmax_Wmv,LOCASIcmean_bias_tmax_Wmv,
               GCMunweightedmean_bias_tmax_Wmv,GCMskillmean_bias_tmax_Wmv,GCMSIhmean_bias_tmax_Wmv,GCMSIcmean_bias_tmax_Wmv)

changevec_pr = c(LOCAunweightedmean_change_pr_Wpr,LOCAskillmean_change_pr_Wpr,LOCASIhmean_change_pr_Wpr,LOCASIcmean_change_pr_Wpr,
              GCMunweightedmean_change_pr_Wpr,GCMskillmean_change_pr_Wpr,GCMSIhmean_change_pr_Wpr,GCMSIcmean_change_pr_Wpr,
              GCMunweightedmean_change_pr_Wtx,GCMskillmean_change_pr_Wtx,GCMSIhmean_change_pr_Wtx,GCMSIcmean_change_pr_Wtx,
              LOCAunweightedmean_change_pr_Wtx,LOCAskillmean_change_pr_Wtx,LOCASIhmean_change_pr_Wtx,LOCASIcmean_change_pr_Wtx,
              LOCAunweightedmean_change_pr_Wmv,LOCAskillmean_change_pr_Wmv,LOCASIhmean_change_pr_Wmv,LOCASIcmean_change_pr_Wmv,
              GCMunweightedmean_change_pr_Wmv,GCMskillmean_change_pr_Wmv,GCMSIhmean_change_pr_Wmv,GCMSIcmean_change_pr_Wmv)

changevec_tmax = c(LOCAunweightedmean_change_tmax_Wpr,LOCAskillmean_change_tmax_Wpr,LOCASIhmean_change_tmax_Wpr,LOCASIcmean_change_tmax_Wpr,
                 GCMunweightedmean_change_tmax_Wpr,GCMskillmean_change_tmax_Wpr,GCMSIhmean_change_tmax_Wpr,GCMSIcmean_change_tmax_Wpr,
                 GCMunweightedmean_change_tmax_Wtx,GCMskillmean_change_tmax_Wtx,GCMSIhmean_change_tmax_Wtx,GCMSIcmean_change_tmax_Wtx,
                 LOCAunweightedmean_change_tmax_Wtx,LOCAskillmean_change_tmax_Wtx,LOCASIhmean_change_tmax_Wtx,LOCASIcmean_change_tmax_Wtx,
                 LOCAunweightedmean_change_tmax_Wmv,LOCAskillmean_change_tmax_Wmv,LOCASIhmean_change_tmax_Wmv,LOCASIcmean_change_tmax_Wmv,
                 GCMunweightedmean_change_tmax_Wmv,GCMskillmean_change_tmax_Wmv,GCMSIhmean_change_tmax_Wmv,GCMSIcmean_change_tmax_Wmv)


changevarvec_pr = c(LOCAunweightedvar_change_pr_Wpr,LOCAskillvar_change_pr_Wpr,LOCASIhvar_change_pr_Wpr,LOCASIcvar_change_pr_Wpr,
                 GCMunweightedvar_change_pr_Wpr,GCMskillvar_change_pr_Wpr,GCMSIhvar_change_pr_Wpr,GCMSIcvar_change_pr_Wpr,
                 LOCAunweightedvar_change_pr_Wtx,LOCAskillvar_change_pr_Wtx,LOCASIhvar_change_pr_Wtx,LOCASIcvar_change_pr_Wtx,
                 GCMunweightedvar_change_pr_Wtx,GCMskillvar_change_pr_Wtx,GCMSIhvar_change_pr_Wtx,GCMSIcvar_change_pr_Wtx,
                 LOCAunweightedvar_change_pr_Wmv,LOCAskillvar_change_pr_Wmv,LOCASIhvar_change_pr_Wmv,LOCASIcvar_change_pr_Wmv,
                 GCMunweightedvar_change_pr_Wmv,GCMskillvar_change_pr_Wmv,GCMSIhvar_change_pr_Wmv,GCMSIcvar_change_pr_Wmv)

changevarvec_tmax = c(LOCAunweightedvar_change_tmax_Wpr,LOCAskillvar_change_tmax_Wpr,LOCASIhvar_change_tmax_Wpr,LOCASIcvar_change_tmax_Wpr,
                    GCMunweightedvar_change_tmax_Wpr,GCMskillvar_change_tmax_Wpr,GCMSIhvar_change_tmax_Wpr,GCMSIcvar_change_tmax_Wpr,
                    LOCAunweightedvar_change_tmax_Wtx,LOCAskillvar_change_tmax_Wtx,LOCASIhvar_change_tmax_Wtx,LOCASIcvar_change_tmax_Wtx,
                    GCMunweightedvar_change_tmax_Wtx,GCMskillvar_change_tmax_Wtx,GCMSIhvar_change_tmax_Wtx,GCMSIcvar_change_tmax_Wtx,
                    LOCAunweightedvar_change_tmax_Wmv,LOCAskillvar_change_tmax_Wmv,LOCASIhvar_change_tmax_Wmv,LOCASIcvar_change_tmax_Wmv,
                    GCMunweightedvar_change_tmax_Wmv,GCMskillvar_change_tmax_Wmv,GCMSIhvar_change_tmax_Wmv,GCMSIcvar_change_tmax_Wmv)


biascolorbar_pr = colorramp(biasvec_pr,colorchoice="bluetored",Blimit=30,use_fixed_scale = TRUE, fixed_scale = c(-1000,1000))
biascolorbar_tmax = colorramp(biasvec_tmax,colorchoice="bluetored",Blimit=30,use_fixed_scale = TRUE, fixed_scale = c(-10,10))

changecolorbar_pr = colorramp(changevec_pr,colorchoice="browntogreen",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(-130,130))
changecolorbar_tmax = colorramp(changevec_tmax,colorchoice="bluetored",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(-6,6))

changevarcolorbar_pr = colorramp(sqrt(changevarvec_pr),colorchoice="whitetogreen",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(0,250))
changevarcolorbar_tmax = colorramp(sqrt(changevarvec_tmax),colorchoice="whitetored",Blimit=30,use_fixed_scale = TRUE,fixed_scale = c(0,3))


pdf("/home/woot0002/DS_ind/biasmaps_pr.pdf",width=20,height=15,onefile=TRUE)
par(mfrow=c(3,4))

for(d in 1:2){
for(i in 1:12){
  
  if(i==1){
    namepart = "Unweighted"
    weightname = ""
    if(d==1) tmp1 = GCMunweightedmean_bias_pr_Wpr;
    if(d==2) tmp1 = LOCAunweightedmean_bias_pr_Wpr;
  }
  if(d==1){
    groupname="GCM"
  }
  if(d==2){
    groupname="LOCA"
  }
  if(i>=2 & i<=4){
    weightname="pr"
  }
  if(i>=6 & i<=8){
    weightname="tmax"
  }
  if(i>=10 & i<=12){
    weightname="mv"
  }
  
  if(i==2 | i==6 | i==10){
    namepart = "Skill Weighted"
    if(i==2 & d==1){
      tmp1 = GCMskillmean_bias_pr_Wpr
    }
    if(i==2 & d==2){
      tmp1 = LOCAskillmean_bias_pr_Wpr
    }
    if(i==6 & d==1){
      tmp1 = GCMskillmean_bias_pr_Wtx
    }
    if(i==6 & d==2){
      tmp1 = LOCAskillmean_bias_pr_Wtx
    }
    if(i==10 & d==1){
      tmp1 = GCMskillmean_bias_pr_Wmv
    }
    if(i==10 & d==2){
      tmp1 = LOCAskillmean_bias_pr_Wmv
    }
  }
  
  if(i==3 | i==7 | i==11){
    namepart = "SI-h Weighted"
    if(i==3 & d==1){
      tmp1 = GCMSIhmean_bias_pr_Wpr
    }
    if(i==3 & d==2){
      tmp1 = LOCASIhmean_bias_pr_Wpr
    }
    if(i==7 & d==1){
      tmp1 = GCMSIhmean_bias_pr_Wtx
    }
    if(i==7 & d==2){
      tmp1 = LOCASIhmean_bias_pr_Wtx
    }
    if(i==11 & d==1){
      tmp1 = GCMSIhmean_bias_pr_Wmv
    }
    if(i==11 & d==2){
      tmp1 = LOCASIhmean_bias_pr_Wmv
    }
  }
  
  if(i==4 | i==8 | i==12){
    namepart = "SI-c Weighted"
    if(i==4 & d==1){
      tmp1 = GCMSIcmean_bias_pr_Wpr
    }
    if(i==4 & d==2){
      tmp1 = LOCASIcmean_bias_pr_Wpr
    }
    if(i==8 & d==1){
      tmp1 = GCMSIcmean_bias_pr_Wtx
    }
    if(i==8 & d==2){
      tmp1 = LOCASIcmean_bias_pr_Wtx
    }
    if(i==12 & d==1){
      tmp1 = GCMSIcmean_bias_pr_Wmv
    }
    if(i==12 & d==2){
      tmp1 = LOCASIcmean_bias_pr_Wmv
    }
  }
  
  if(i==5 | i==9){
    frame()
  } else {
    testsfc = list(x=lon,y=lat,z=tmp1)
    surface(testsfc,type="I",zlim=biascolorbar_pr[[1]],col=biascolorbar_pr[[3]],breaks=biascolorbar_pr[[2]],xlab="Longitude",ylab="Latitude",main=paste(groupname,weightname,namepart,sep=" "))
    map("state",add=TRUE)
    text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  }

}
}
dev.off()


pdf("/home/woot0002/DS_ind/biasmaps_tmax.pdf",width=20,height=15,onefile=TRUE)
par(mfrow=c(3,4))

for(d in 1:2){
  for(i in 1:12){
    
    if(i==1){
      namepart = "Unweighted"
      weightname = ""
      if(d==1) tmp1 = GCMunweightedmean_bias_tmax_Wpr;
      if(d==2) tmp1 = LOCAunweightedmean_bias_tmax_Wpr;
    }
    if(d==1){
      groupname="GCM"
    }
    if(d==2){
      groupname="LOCA"
    }
    if(i>=2 & i<=4){
      weightname="pr"
    }
    if(i>=6 & i<=8){
      weightname="tmax"
    }
    if(i>=10 & i<=12){
      weightname="mv"
    }
    
    if(i==2 | i==6 | i==10){
      namepart = "Skill Weighted"
      if(i==2 & d==1){
        tmp1 = GCMskillmean_bias_tmax_Wpr
      }
      if(i==2 & d==2){
        tmp1 = LOCAskillmean_bias_tmax_Wpr
      }
      if(i==6 & d==1){
        tmp1 = GCMskillmean_bias_tmax_Wtx
      }
      if(i==6 & d==2){
        tmp1 = LOCAskillmean_bias_tmax_Wtx
      }
      if(i==10 & d==1){
        tmp1 = GCMskillmean_bias_tmax_Wmv
      }
      if(i==10 & d==2){
        tmp1 = LOCAskillmean_bias_tmax_Wmv
      }
    }
    
    if(i==3 | i==7 | i==11){
      namepart = "SI-h Weighted"
      if(i==3 & d==1){
        tmp1 = GCMSIhmean_bias_tmax_Wpr
      }
      if(i==3 & d==2){
        tmp1 = LOCASIhmean_bias_tmax_Wpr
      }
      if(i==7 & d==1){
        tmp1 = GCMSIhmean_bias_tmax_Wtx
      }
      if(i==7 & d==2){
        tmp1 = LOCASIhmean_bias_tmax_Wtx
      }
      if(i==11 & d==1){
        tmp1 = GCMSIhmean_bias_tmax_Wmv
      }
      if(i==11 & d==2){
        tmp1 = LOCASIhmean_bias_tmax_Wmv
      }
    }
    
    if(i==4 | i==8 | i==12){
      namepart = "SI-c Weighted"
      if(i==4 & d==1){
        tmp1 = GCMSIcmean_bias_tmax_Wpr
      }
      if(i==4 & d==2){
        tmp1 = LOCASIcmean_bias_tmax_Wpr
      }
      if(i==8 & d==1){
        tmp1 = GCMSIcmean_bias_tmax_Wtx
      }
      if(i==8 & d==2){
        tmp1 = LOCASIcmean_bias_tmax_Wtx
      }
      if(i==12 & d==1){
        tmp1 = GCMSIcmean_bias_tmax_Wmv
      }
      if(i==12 & d==2){
        tmp1 = LOCASIcmean_bias_tmax_Wmv
      }
    }
    
    if(i==5 | i==9){
      frame()
    } else {
      testsfc = list(x=lon,y=lat,z=tmp1)
      surface(testsfc,type="I",zlim=biascolorbar_tmax[[1]],col=biascolorbar_tmax[[3]],breaks=biascolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude",main=paste(groupname,weightname,namepart,sep=" "))
      map("state",add=TRUE)
      text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    }
    
  }
}
dev.off()


pdf("/home/woot0002/DS_ind/changemaps_pr.pdf",width=20,height=15,onefile=TRUE)
par(mfrow=c(3,4))

for(d in 1:2){
  for(i in 1:12){
    
    if(i==1){
      namepart = "Unweighted"
      weightname = ""
      if(d==1) tmp1 = GCMunweightedmean_change_pr_Wpr;
      if(d==2) tmp1 = LOCAunweightedmean_change_pr_Wpr;
    }
    if(d==1){
      groupname="GCM"
    }
    if(d==2){
      groupname="LOCA"
    }
    if(i>=2 & i<=4){
      weightname="pr"
    }
    if(i>=6 & i<=8){
      weightname="tmax"
    }
    if(i>=10 & i<=12){
      weightname="mv"
    }
    
    if(i==2 | i==6 | i==10){
      namepart = "Skill Weighted"
      if(i==2 & d==1){
        tmp1 = GCMskillmean_change_pr_Wpr
      }
      if(i==2 & d==2){
        tmp1 = LOCAskillmean_change_pr_Wpr
      }
      if(i==6 & d==1){
        tmp1 = GCMskillmean_change_pr_Wtx
      }
      if(i==6 & d==2){
        tmp1 = LOCAskillmean_change_pr_Wtx
      }
      if(i==10 & d==1){
        tmp1 = GCMskillmean_change_pr_Wmv
      }
      if(i==10 & d==2){
        tmp1 = LOCAskillmean_change_pr_Wmv
      }
    }
    
    if(i==3 | i==7 | i==11){
      namepart = "SI-h Weighted"
      if(i==3 & d==1){
        tmp1 = GCMSIhmean_change_pr_Wpr
      }
      if(i==3 & d==2){
        tmp1 = LOCASIhmean_change_pr_Wpr
      }
      if(i==7 & d==1){
        tmp1 = GCMSIhmean_change_pr_Wtx
      }
      if(i==7 & d==2){
        tmp1 = LOCASIhmean_change_pr_Wtx
      }
      if(i==11 & d==1){
        tmp1 = GCMSIhmean_change_pr_Wmv
      }
      if(i==11 & d==2){
        tmp1 = LOCASIhmean_change_pr_Wmv
      }
    }
    
    if(i==4 | i==8 | i==12){
      namepart = "SI-c Weighted"
      if(i==4 & d==1){
        tmp1 = GCMSIcmean_change_pr_Wpr
      }
      if(i==4 & d==2){
        tmp1 = LOCASIcmean_change_pr_Wpr
      }
      if(i==8 & d==1){
        tmp1 = GCMSIcmean_change_pr_Wtx
      }
      if(i==8 & d==2){
        tmp1 = LOCASIcmean_change_pr_Wtx
      }
      if(i==12 & d==1){
        tmp1 = GCMSIcmean_change_pr_Wmv
      }
      if(i==12 & d==2){
        tmp1 = LOCASIcmean_change_pr_Wmv
      }
    }
    
    if(i==5 | i==9){
      frame()
    } else {
      testsfc = list(x=lon,y=lat,z=tmp1)
      surface(testsfc,type="I",zlim=changecolorbar_pr[[1]],col=changecolorbar_pr[[3]],breaks=changecolorbar_pr[[2]],xlab="Longitude",ylab="Latitude",main=paste(groupname,weightname,namepart,sep=" "))
      map("state",add=TRUE)
      text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    }
    
  }
}
dev.off()

pdf("/home/woot0002/DS_ind/changemaps_tmax.pdf",width=20,height=15,onefile=TRUE)
par(mfrow=c(3,4))

for(d in 1:2){
  for(i in 1:12){
    
    if(i==1){
      namepart = "Unweighted"
      weightname = ""
      if(d==1) tmp1 = GCMunweightedmean_change_tmax_Wpr;
      if(d==2) tmp1 = LOCAunweightedmean_change_tmax_Wpr;
    }
    if(d==1){
      groupname="GCM"
    }
    if(d==2){
      groupname="LOCA"
    }
    if(i>=2 & i<=4){
      weightname="pr"
    }
    if(i>=6 & i<=8){
      weightname="tmax"
    }
    if(i>=10 & i<=12){
      weightname="mv"
    }
    
    if(i==2 | i==6 | i==10){
      namepart = "Skill Weighted"
      if(i==2 & d==1){
        tmp1 = GCMskillmean_change_tmax_Wpr
      }
      if(i==2 & d==2){
        tmp1 = LOCAskillmean_change_tmax_Wpr
      }
      if(i==6 & d==1){
        tmp1 = GCMskillmean_change_tmax_Wtx
      }
      if(i==6 & d==2){
        tmp1 = LOCAskillmean_change_tmax_Wtx
      }
      if(i==10 & d==1){
        tmp1 = GCMskillmean_change_tmax_Wmv
      }
      if(i==10 & d==2){
        tmp1 = LOCAskillmean_change_tmax_Wmv
      }
    }
    
    if(i==3 | i==7 | i==11){
      namepart = "SI-h Weighted"
      if(i==3 & d==1){
        tmp1 = GCMSIhmean_change_tmax_Wpr
      }
      if(i==3 & d==2){
        tmp1 = LOCASIhmean_change_tmax_Wpr
      }
      if(i==7 & d==1){
        tmp1 = GCMSIhmean_change_tmax_Wtx
      }
      if(i==7 & d==2){
        tmp1 = LOCASIhmean_change_tmax_Wtx
      }
      if(i==11 & d==1){
        tmp1 = GCMSIhmean_change_tmax_Wmv
      }
      if(i==11 & d==2){
        tmp1 = LOCASIhmean_change_tmax_Wmv
      }
    }
    
    if(i==4 | i==8 | i==12){
      namepart = "SI-c Weighted"
      if(i==4 & d==1){
        tmp1 = GCMSIcmean_change_tmax_Wpr
      }
      if(i==4 & d==2){
        tmp1 = LOCASIcmean_change_tmax_Wpr
      }
      if(i==8 & d==1){
        tmp1 = GCMSIcmean_change_tmax_Wtx
      }
      if(i==8 & d==2){
        tmp1 = LOCASIcmean_change_tmax_Wtx
      }
      if(i==12 & d==1){
        tmp1 = GCMSIcmean_change_tmax_Wmv
      }
      if(i==12 & d==2){
        tmp1 = LOCASIcmean_change_tmax_Wmv
      }
    }
    
    if(i==5 | i==9){
      frame()
    } else {
      testsfc = list(x=lon,y=lat,z=tmp1)
      surface(testsfc,type="I",zlim=changecolorbar_tmax[[1]],col=changecolorbar_tmax[[3]],breaks=changecolorbar_tmax[[2]],xlab="Longitude",ylab="Latitude",main=paste(groupname,weightname,namepart,sep=" "))
      map("state",add=TRUE)
      text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
      text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
    }
    
  }
}
dev.off()






pdf("/home/woot0002/DS_ind/changemaps.pdf",width=10,height=25)
par(mfrow=c(5,2))
for(i in 1:5){
  
  if(i==1){
    namepart = "Unweighted"
    tmp1 = GCMunweightedmean_change
    tmp2 = LOCAunweightedmean_change
  }
  if(i==2){
    namepart = "Skill Weighted"
    tmp1 = GCMskillmean_change
    tmp2 = LOCAskillmean_change
  }
  if(i==3){
    namepart = "SI-h Weighted"
    tmp1 = GCMSIhmean_change
    tmp2 = LOCASIhmean_change
  }
  if(i==4){
    namepart = "SI-c Weighted"
    tmp1 = GCMSIcmean_change
    tmp2 = LOCASIcmean_change
  }
  if(i==5){
    namepart = "BMA Weighted"
    tmp1 = GCMBMAmean_change
    tmp2 = LOCABMAmean_change
  }
  
  testsfc = list(x=lon-360,y=lat,z=tmp1)
  surface(testsfc,type="I",zlim=changecolorbar[[1]],col=changecolorbar[[3]],breaks=changecolorbar[[2]],xlab="Longitude",ylab="Latitude",main=paste("CMIP5: ",namepart,sep=""))
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc = list(x=lon-360,y=lat,z=tmp2)
  surface(testsfc,type="I",zlim=changecolorbar[[1]],col=changecolorbar[[3]],breaks=changecolorbar[[2]],xlab="Longitude",ylab="Latitude",main=paste("LOCA: ",namepart,sep=""))
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}
dev.off()


pdf("/home/woot0002/DS_ind/changesdmaps.pdf",width=10,height=25)
par(mfrow=c(5,2))
for(i in 1:5){
  
  if(i==1){
    namepart = "Unweighted"
    tmp1 = sqrt(GCMunweightedvar_change)
    tmp2 = sqrt(LOCAunweightedvar_change)
  }
  if(i==2){
    namepart = "Skill Weighted"
    tmp1 = sqrt(GCMskillvar_change)
    tmp2 = sqrt(LOCAskillvar_change)
  }
  if(i==3){
    namepart = "SI-h Weighted"
    tmp1 = sqrt(GCMSIhvar_change)
    tmp2 = sqrt(LOCASIhvar_change)
  }
  if(i==4){
    namepart = "SI-c Weighted"
    tmp1 = sqrt(GCMSIcvar_change)
    tmp2 = sqrt(LOCASIcvar_change)
  }
  if(i==5){
    namepart = "BMA Weighted"
    tmp1 = sqrt(GCMBMAvar_change)
    tmp2 = sqrt(LOCABMAvar_change)
  }
  
  testsfc = list(x=lon-360,y=lat,z=tmp1)
  surface(testsfc,type="I",zlim=changevarcolorbar[[1]],col=changevarcolorbar[[3]],breaks=changevarcolorbar[[2]],xlab="Longitude",ylab="Latitude",main=paste("CMIP5: ",namepart,sep=""))
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp1,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
  testsfc = list(x=lon-360,y=lat,z=tmp2)
  surface(testsfc,type="I",zlim=changevarcolorbar[[1]],col=changevarcolorbar[[3]],breaks=changevarcolorbar[[2]],xlab="Longitude",ylab="Latitude",main=paste("LOCA: ",namepart,sep=""))
  map("state",add=TRUE)
  text(-108.85,28.45,labels=paste("MAX = ",round(max(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,27.45,labels=paste("MEAN = ",round(mean(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  text(-108.85,26.45,labels=paste("MIN = ",round(min(tmp2,na.rm=TRUE),1),sep=""),cex=1,pos = 4)
  
}
dev.off()
