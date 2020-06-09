####
# Spring Phenology R test code - synthetic
# AMW 6/14/2019
#
# This script demonstrates appropriate use for calculations of first leaf, first bloom, damage index, and false springs using the functions
# in springpheno.R. Specific to vector analyses.
#
#  Questions on code should be emailed to amwootte@ou.edu

source("/data2/3to5/I35/scripts/analysisfunctions.R")
#source("colorramp.R")
library(ncdf4)
library(maps)
library(fields)
library(sp)
setwd("/home/woot0002") # edit or remove paths for scripts and files as needed
source("scripts/falsesprings/springpheno_v0.3.R")

load("lfs_tmin.Rdata")

RESULTS = calc_si(TMAX,TMIN,latval)
  
  

####










