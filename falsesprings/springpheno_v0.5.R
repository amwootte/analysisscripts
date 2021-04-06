####
# Spring Phenology R functions v. 0.5
# AMW 12/9/2020
#
# This script contains functions describing numerous spring onset indices and false spring
# exposure indices based on several papers. The following specifically are included here:
# 1) Extended Spring Indices (SI-x) and supporting functions (Ault et al. 2015; Schwarz et al. 2006; Schwartz et al. 2013)
# 2) Early and Late False Spring (EFS and LFS) - Allstadt et al. 2015
# 3) False Spring Exposure Index (FSEI) - Peterson and Abatzoglou, 2014
#    Early FSEI and Late FSEI are calculated using those indicators from the Allstadt et al. 2015 study.
#
# Notes by version
#
# v0.1
# - The chill hours and chill date functions are included based on SI-x, but it is my rough interpretation
#   given a lack of documentation in literature. Use with caution.
# - The SI-x are created based on translation from the original Matlab code from Ault et al. (2015).
#   However, there are a couple improvements made. First, the FLI and FBI calculations incorporate the bug
#   fix by Allstadt et al. (2015). Second, the original Matlab code had a function which read in precalculated solar
#   declination angle. In this R function, I've incorporated code for calculating solar declination from R's solrad package.
# - The FSEI formula in the functions is based on Peterson and Abatzoglou (2014), however, the definition
#   of a false spring is based on Allstadt et al. (2015). The FSEI function can be used to determine false
#   spring exposure with both EFS and LFS.
#
#
# v0.2 6/24/2019
# - incorporated the function soldec, which brings in the precalculated solar declination values from Ault et al (2015).
#   As a result, the declination function is no longer called.
# - Added the function calc_si, which (like the original Matlab) serves as a wrapper and calculates all SI-x and FSEI values for one or multiple years.
# - There are situations where the last freeze during the year occurs in the winter following first leaf. This is another no damage situation, as plants would be dormant.
#   As such, the damageindex function now contains an added conditional forcing the damage index to zero if the last freeze occurs in the winter following first leaf.
# - Tests with Livneh in the 3^5 domain demonstrated a situation where some locations will have years where tasmin is always >= the freeze value (frzval).
#   In this situation, there is no damage or false springs, so the values for both are set to zero (i.e. no freeze damage, no false spring occurred),
#   using an if statement at the end of the calc_si function.
#
#
# v0.3 8/16/2019
# - There is a chance that the a freeze will occur after the first bloom, but not occur between the first leaf and first bloom.
#   In the prior version, the falsesprings function would have counted this occurrence as both a late false spring and an early false spring.
#   Counting the event as both an early and late false spring would be logically inconsistent. A correction for this has been applied to the falsesprings function.
#   Essentially the correction ensures that the falsesprings function will only mark that those freeze days occurring one week after first leaf, but before first bloom are considered.
# - In preparation for implementation in a formal R package, the unused declination function has been deleted.
# - No other changes from prior versions.
# - Cleaned up the daylength function to drop the unused yearlength argument
#
#
# v0.4 6/9/2020
# - corrected bug in calc_si (length(tminNAidx)>0), where the if statement originally referred to tmaxNAidx.
#
# v0.5 12/9/2020
# - added a loess option for imputing missing data to the calc_si function. Also corrected an overwriting glitch
#   that was present in calc_si. TMAX was being overwritten by TMAXmat early in the function. Could lead to 
#   corrections for missing values being ignored.
#
#  Questions on code should be emailed to amwootte@ou.edu
#
#########################################
#########################################

# Calculate Spring Indices - Given appropriate inputs, this function will calculate all the spring indices
calc_si = function(TMAX,TMIN,lat,missingcalc="mean"){
  # INPUTS
  # TMAX, TMIN - Matrices of tasmax and tasmin values which are 366 rows x N columns (matching number of years).
  #              The matrix should always have 366 days, as the function can fill in gaps (such as missing leap days) on its own.
  #              Input units should also be degrees F, though it will warn and auto convert if TMAX and TMIN are provided in degrees K.
  # lat - scalar latitude of location
  # missingcalc - options for imputing missing values in TMAX and TMIN. Current options include "mean" (default)
  #               and "loess". The "mean" option is the original imputting put forward in the Matlab package.
  #               The "loess" option will direct calc_si to use loess to impute any missing values.
  #
  # OUTPUTS - in a list format variable
  # FLImat - Matrix of first leaf which is N rows (matching number of years) by 4 columns
  #          The 4 columns correspond to the mean first leaf and the first leaf for plants 1-3 [mean,plant 1, plant 2, plant 3]
  # FBImat - Alike to the FLImat, but for first bloom.
  # DMGmat - Alike to the FLImat, but for the damage index.
  # lastfreeze - Vector containing the day of year value (1-366) for the day of last freeze in spring. Should be length of N years.
  # FSmat - Matrix containing false spring indicators. Should be N rows (matching number of years) by 2 columns.
  #         The two columns represent each of the false spring indicators [early false spring, late false spring]
  #         A value of 1 in either column represents that a false spring occurred. A value of 0 represents that a false spring did not occur
  # FSEImat - False Spring Exposure Index (FSEI) values. This should be a vector of length 2.
  #           There is one value of FSEI each for the early and late false springs. [early FSEI, late FSEI]
  #
  # Note: Thresholds for base temperature (baset) and freezing temperature (frzval) are fixed in this code. This matches Ault et al. (2015),
  # but it could be changed to an extra argument for this function later.
  # Note: FSmat and FSEImat are current calculated using the first values of FLImat and FBImat. This is also consistent with prior code, but
  # could be altered to provide the FSmat and FSEImat based on all four values of FLImat and FBImat.
  
  #TMAX = TMAXmat # These three lines were not commented out in version 0.4
  #TMIN = TMINmat
  #lat = latin
  
  if(is.matrix(TMAX)==FALSE | is.matrix(TMIN)==FALSE){
    stop("TMAX or TMIN are not matrices",.call=TRUE)
  }
  
  if(is.numeric(TMAX)==FALSE | is.numeric(TMIN)==FALSE | is.numeric(lat)==FALSE){
    stop("TMAX, TMIN, or lat are not numeric",.call=TRUE)
  }
  
  if(all(TMAX[which(is.na(TMAX)==FALSE)]>200)==TRUE & all(TMIN[which(is.na(TMIN)==FALSE)]>200)==TRUE){
    warning("TMAX and TMIN appear to be in degrees K, so I converted them to degrees F")
    TMAX = (TMAX-273.15)*9/5+32
    TMIN = (TMIN-273.15)*9/5+32
  }
  
  baset = 31 # baset is 31F, 0C, 273.15K
  frzval = 28 # -2.2C, 28F - hard freeze temperature
  
  # initialize matrices for outputs
  FLImat = FBImat = DMGmat = matrix(NA,nrow=ncol(TMAX),ncol=4)
  FSmat = matrix(NA,nrow=ncol(TMAX),ncol=2)
  lastfreeze = c()
  
  for(i in 1:ncol(TMAX)){
    
    tasmax = as.vector(TMAX[,i]) # getting tasmin and tasmax into vectors
    tasmin = as.vector(TMIN[,i])
    
    
    if(missingcalc=="mean"){
      message("using monthly mean to impute missing values")
    # if values are missing from tasmax or tasmin, the average of the 30 days is used.
    # this section fills in missing values, similarly to the code of Ault et al. 2015
    for(x in seq(1,211,by=30)){
      start = x
      end = x+29
      flagm = 0
      
      tmaxNAidx = which(is.na(tasmax[start:end])==TRUE)
      tminNAidx = which(is.na(tasmin[start:end])==TRUE)
      
      cntmax = length(tasmax[start:end])-length(tmaxNAidx)
      cntmin = length(tasmin[start:end])-length(tminNAidx)
      
      tmax = tasmax[start:end]
      tmin = tasmin[start:end]
      idx = start:end
      if(length(tmaxNAidx)>0){
        avgmax = mean(tmax[-tmaxNAidx])
        tasmax[idx[tmaxNAidx]]=avgmax
      }
      if(length(tminNAidx)>0){
        avgmin = mean(tmin[-tminNAidx])
        tasmin[idx[tminNAidx]]=avgmin
      }
      if(cntmax < 20 | cntmin <20) flagm = 1
      if(flagm == 1){
        warning(paste("There is a lot of missing data in TMAX and TMIN for column ",i,sep=""))
      }
    }
    }
    
    if(missingcalc=="loess"){
      message("using loess to impute missing values")
      tstart = -(nrow(TMAX)/2):0
      tstartuse = 1:184
      t=1:nrow(TMAX)
      tend = seq(367,by=1,length.out=184)
      tenduse = seq(366,by=-1,length.out=184)
      tenduse = tenduse[order(tenduse)]
      tstring = c(tstart,t,tend)
      
      if(i==1){
        tmaxseries = c(TMAX[tenduse,i],tasmax,TMAX[tstartuse,(i+1)])
        tminseries = c(TMIN[tenduse,i],tasmin,TMIN[tstartuse,(i+1)])
      }
      if(i>1 & i<ncol(TMAX)){
        tmaxseries = c(TMAX[tenduse,(i-1)],tasmax,TMAX[tstartuse,(i+1)])
        tminseries = c(TMIN[tenduse,(i-1)],tasmin,TMIN[tstartuse,(i+1)])
      }
      if(i==ncol(TMAX)){
        tmaxseries = c(TMAX[tenduse,(i-1)],tasmax,TMAX[tstartuse,(i-1)])
        tminseries = c(TMIN[tenduse,(i-1)],tasmin,TMIN[tstartuse,(i-1)])
      }
      
      notNAidx = which(is.na(tmaxseries)==FALSE)
      isNAidx = which(is.na(tmaxseries)==TRUE)
      loessfit = loess(tmaxseries[notNAidx]~tstring[notNAidx],span=0.03)
      tmaxseries[isNAidx] = predict(loessfit,newdata=tstring[isNAidx])
      tasmax = tmaxseries[which(tstring>=1 & tstring<=366)]
      
      notNAidx = which(is.na(tminseries)==FALSE)
      isNAidx = which(is.na(tminseries)==TRUE)
      loessfit = loess(tminseries[notNAidx]~tstring[notNAidx],span=0.03)
      tminseries[isNAidx] = predict(loessfit,newdata=tstring[isNAidx])
      tasmin = tminseries[which(tstring>=1 & tstring<=366)]
      
      check1idx = which(is.na(tasmax[1:31])==TRUE)
      if(length(check1idx)>0){
        tasmax[check1idx] = mean(tasmax[1:31],na.rm=TRUE)
      }
      check1idx = which(is.na(tasmin[1:31])==TRUE)
      if(length(check1idx)>0){
        tasmin[check1idx] = mean(tasmin[1:31],na.rm=TRUE)
      }
      
      if(is.na(tasmax[length(tasmax)])==TRUE){
        tasmax[length(tasmax)] = tasmax[(length(tasmax)-1)]
      }
      
      if(is.na(tasmin[length(tasmin)])==TRUE){
        tasmin[length(tasmin)] = tasmin[(length(tasmin)-1)]
      }
      
      if(any(is.na(tasmax)==TRUE)==TRUE){
        message("Still has missing values in tasmax") 
        print(tasmax)
      }
      if(any(is.na(tasmin)==TRUE)==TRUE){
        message("Still has missing values in tasmin")  
        print(tasmin)
      }
      
      }
    
    
    DOY = 1:length(tasmax) # the example I have here is with a full year of data, so you may which to calculate with less than one full year.
    daystop = max(DOY) # for daylength, you can have daylength calculated until a particular DOY. I just have it run the whole year here.
    
    dlen = daylength(daystop,lat=lat) # Calculating daylength
    
    FLImat[i,2] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=1,type="leaf",plant=1) # Following with the Matlab code in Ault et al. 2015, the FLI numbers are ordered such that the mean FLI is the first number in the output vector
    FLImat[i,3] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=1,type="leaf",plant=2) # and the three following numbers are the FLI for each plant used to calculate first leaf.
    FLImat[i,4] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=1,type="leaf",plant=3)
    FLImat[i,1] = round(mean(FLImat[i,2:4],na.rm=TRUE))
    
    if(is.na(FLImat[i,2])==FALSE){FBImat[i,2] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=FLImat[i,2],type="bloom",plant=1)} # FBI has a similar vector out to FLI above. FBI also takes the FLI for each plant as it's reference date.
    if(is.na(FLImat[i,3])==FALSE){FBImat[i,3] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=FLImat[i,3],type="bloom",plant=2)} # This is shown in the Matlab code of Ault et al. 2015
    if(is.na(FLImat[i,4])==FALSE){FBImat[i,4] = leafindex(tasmax=tasmax,tasmin=tasmin,daylen=dlen,baset=baset,refdate=FLImat[i,4],type="bloom",plant=3)}
    FBImat[i,1] = round(mean(FBImat[i,2:4],na.rm=TRUE))
    
    freezeinfo = freezedates(tasmin=tasmin,frzval=frzval,DOY=DOY)
    lastfreeze[i] = freezeinfo$lastfreeze
    freezedata = subset(freezeinfo[[4]],DOY>=1 & DOY<=180) # table containing the tasmin, and day of year for the hard freeze.
    
    if(is.na(lastfreeze[i])==FALSE){
      if(is.na(FLImat[i,2])==FALSE){DMGmat[i,2] = damageindex(leafidx=FLImat[i,2],lastfreeze=freezeinfo$lastfreeze)} # per the Schwarz and Ault papers, the last freeze is all that is needed calculate damage.
      if(is.na(FLImat[i,3])==FALSE){DMGmat[i,3] = damageindex(leafidx=FLImat[i,3],lastfreeze=freezeinfo$lastfreeze)} # also with a similar vector output to FLI and FBI.
      if(is.na(FLImat[i,4])==FALSE){DMGmat[i,4] = damageindex(leafidx=FLImat[i,4],lastfreeze=freezeinfo$lastfreeze)}
      DMGmat[i,1] = round(mean(DMGmat[i,2:4],na.rm=TRUE))
      
      fs = falsesprings(SI=c(FLImat[i,1],FBImat[i,1]),freezedata = freezedata) # takes the freezedata table to determine if there's an early false spring or a late false spring.
      FSmat[i,1] = fs$EFS
      FSmat[i,2] = fs$LFS
    } else {
      message("Freeze didn't happen this year, so no damage or false springs")
      DMGmat[i,1:4] = 0
      FSmat[i,1:2] = 0
    }
    #message("Calcs finished for year: ",i," / ",ncol(TMAX))
  }
  
  FSEImat = apply(FSmat,2,FSEI)
  list(FLImat=FLImat,FBImat=FBImat,DMGmat=DMGmat,lastfreeze=lastfreeze,FSmat=FSmat,FSEImat=FSEImat)
}


####
# Solar declination - from the original matlab package, precalculated
soldec = function(DOY){
  # INPUTS
  # DOY - scalar day of year, from 1 to 366 or 1 to 365 depending on the calendar of the date
  #
  # OUTPUTS
  # spring climatological calendar day of year for the location and calendar day of year matching the appropriate solar declination
  
  decvals=c(-23.16,-23.08,-23,-22.92,-22.83,-22.73,-22.62,-22.5,-22.38,-22.24,-22.11,-21.96,-21.81,-21.65,-21.48,-21.31,-21.13,-20.94,
            -20.75,-20.55,-20.34,-20.13,-19.91,-19.68,-19.45,-19.21,-18.97,-18.72,-18.46,-18.2,-17.94,-17.68,-17.39,-17.11,-16.83,-16.53,-16.23,
            -15.93,-15.62,-15.31,-15,-14.68,-14.36,-14.03,-13.7,-13.36,-13.03,-12.68,-12.34,-11.99,-11.64,-11.28,-10.93,-10.57,-10.2,-9.84,-9.47,
            -9.1,-8.73,-8.35,-7.98,-7.6,-7.22,-6.83,-6.45,-6.06,-5.68,-5.29,-4.9,-4.51,-4.12,-3.72,-3.33,-2.93,-2.54,-2.15,-1.75,-1.36,-0.97,-0.57,
            -0.17,0.22,0.62,1.01,1.41,1.8,2.19,2.58,2.98,3.37,3.76,4.14,4.53,4.92,5.3,5.68,6.06,6.44,6.82,7.19,7.57,7.94,8.31,8.67,9.04,9.4,9.76,
            10.11,10.47,10.82,11.16,11.51,11.85,12.19,12.52,12.85,13.18,13.5,13.83,14.14,14.45,14.76,15.07,15.37,15.67,15.95,16.24,16.53,16.81,17.08,
            17.35,17.61,17.87,18.13,18.38,18.62,18.87,19.09,19.32,19.54,19.76,19.97,20.18,20.38,20.57,20.76,20.94,21.12,21.29,21.46,21.61,21.77,21.91,
            22.05,22.18,22.31,22.43,22.54,22.65,22.75,22.84,22.93,23.01,23.08,23.15,23.21,23.26,23.31,23.35,23.38,23.41,23.43,23.44,23.44,23.43,23.41,
            23.38,23.35,23.31,23.26,23.21,23.15,23.08,23.01,22.93,22.84,22.75,22.65,22.54,22.43,22.31,22.18,22.05,21.91,21.77,21.61,21.46,21.29,21.12,
            20.94,20.76,20.57,20.38,20.18,19.97,19.76,19.54,19.32,19.09,18.87,18.62,18.38,18.13,17.87,17.61,17.35,17.08,16.81,16.53,16.24,15.95,15.67,
            15.37,15.07,14.76,14.45,14.14,13.83,13.5,13.18,12.85,12.52,12.19,11.85,11.51,11.16,10.82,10.47,10.11,9.76,9.4,9.04,8.67,8.31,7.94,7.57,7.19,
            6.82,6.44,6.06,5.68,5.3,4.92,4.53,4.14,3.76,3.37,2.98,2.58,2.19,1.8,1.41,1.01,0.62,0.22,-0.17,-0.57,-0.97,-1.36,-1.75,-2.15,-2.54,-2.93,-3.33,
            -3.72,-4.12,-4.51,-4.9,-5.29,-5.68,-6.06,-6.45,-6.83,-7.22,-7.6,-7.98,-8.35,-8.73,-9.1,-9.47,-9.84,-10.2,-10.57,-10.93,-11.28,-11.64,-11.99,-12.34,
            -12.68,-13.03,-13.36,-13.7,-14.03,-14.36,-14.68,-15,-15.31,-15.62,-15.93,-16.23,-16.53,-16.83,-17.11,-17.39,-17.68,-17.94,-18.2,-18.46,-18.72,-18.97,
            -19.21,-19.45,-19.68,-19.91,-20.13,-20.34,-20.55,-20.75,-20.94,-21.13,-21.31,-21.48,-21.65,-21.81,-21.96,-22.11,-22.24,-22.38,-22.5,-22.62,-22.73,
            -22.83,-22.92,-23,-23.08,-23.16,-23.23,-23.29,-23.34,-23.38,-23.42,-23.45,-23.47,-23.48,-23.49,-23.5,-23.5,-23.49,-23.48,-23.47,-23.45,-23.42,
            -23.38,-23.34,-23.29,-23.23)
  
  dayvals = c(307:366,1:306)
  
  dayvals[DOY]
}


###

# Daylength calculator - supporting function, but is called directly for leaf index calculations
daylength = function(daystop,lat){
  # INPUTS
  # daystop - scalar day of year at which to end calculations
  # lat - latitude for location of interest in decimal degrees
  #
  # OUTPUTS
  # DAYLEN - vector of length equal to input daystop describing the total hours of daylight per day.
  
  DAYLEN = c()
  for(i in 1:daystop){
    if(lat < 40){ # high vs. low latitude calculations for daylength
      DLL = 12.14+3.34*tan(lat*pi/180)*cos(0.0172*soldec(i)-1.95)
    } else {
      DLL = 12.25 + (1.6164+1.7643*(tan(lat*pi/180))^2)*cos(0.0172*soldec(i)-1.95)
    }
    if(DLL < 1) DLL=1 # for very high latitudes, these ifs allow for having at least 1 hour of day or night.
    if(DLL> 23) DLL = 23
    DAYLEN[i]=DLL
  }
  DAYLEN
}

###

# Chill Hours calculator - supporting function for SI-x chill indices, not called directly.
chillh = function(tmax,tmin,daylen,baset){
  # INPUTS
  # tmax - scalar daily high temperature for a given day
  # tmin - scalar daily low temperature for a given day
  # daylen - scalar daylength for a given day (calculated from daylength function)
  # baset - base temperature for determining chill hours. Typically this is 7.2C, but can be plant specific.
  #
  # OUTPUTS
  # CHOUR - total chill hours for a given day. That is, number of hours the temperature fell below baset.
  
  DTR = tmax-tmin # diurnal temperature range
  IDL = floor(daylen) # rounds off daylen to an integer value
  
  ###
  # Estimates hourly temperature curves during daylight hours
  Thr = c()
  Thr[1] = tmin
  Thr[2:(IDL+1)]=DTR*sin(pi/(daylen+4)*(1:IDL))+tmin
  
  ###
  # Estimates hourly temperature curve during night hours and puts them together with daylight hours.
  TS1 = DTR*sin(pi/(daylen+4)*daylen)+tmin
  if(TS1 <=0) TS1 = 0.01
  Thr[(IDL+2):24] = TS1-(TS1-tmin)/(log(24-daylen))*log((1:(23-IDL)))
  
  ###
  # chill hours calculation
  CHRS = baset-Thr
  CHRSflag = ifelse(CHRS<0,0,1)
  CHOUR = sum(CHRSflag)
  #list(CHOUR,Thr) # you can uncomment this line and comment out the next to see the hourly temperature curve.
  CHOUR
}

###

# Chill date calculator - provides output to calculate the Composite Chill Date for SI-x.
chilldate = function(tasmax,tasmin,DOY,daylen,baset,plant=1){
  # INPUTS
  # tasmax - vector of daily high temperatures for a given year
  # tasmin - vector of daily low temperatures for a given year
  # daylen - vector of daylengths for a given year (calculated from daylength function).
  # DOY - vector of day of year values for a given year (1:365 or 1:366 depending on calendar)
  # *** tasmax, tasmin, daylen, and DOY must be the same length.
  # baset - base temperature for determining chill hours. Typically this is 7.2C, but can be plant specific.
  # plant - signifies plant type 1, 2, or 3. plant = 1 for lilac, plant = 2 for arnold red, plant =3 for zabelli
  #
  # OUTPUTS
  # chillDOY - scalar, estimated day of year the chill requirement for the specified plant is met.
  #
  # Many questions regarding how this is calculated, use at your own risk.
  
  ###
  # adjusts the calendar to begin from July 1 and wrap to June 30.
  chillframe = data.frame(tasmax,tasmin,daylen,DOY)
  if(nrow(chillframe)==365){
    chillframe$DOYadj = ifelse(DOY<=180,DOY+185,DOY-180)
  } else {
    chillframe$DOYadj = ifelse(DOY<=180,DOY+186,DOY-180)
  }
  chillframe = chillframe[order(chillframe$DOYadj),]
  
  ###
  # Calculates chill hours for each day and accumulated chill hours
  chillframe$CH = NA
  chillframe$ACH = NA
  for(i in 1:nrow(chillframe)){
    chillframe$CH[i] = chillh(chillframe$tasmax[i],chillframe$tasmin[i],chillframe$daylen[i],baset=baset)
    if(i==1) chillframe$ACH[i]=chillframe$CH[i]
    if(i>1) chillframe$ACH[i]=chillframe$ACH[(i-1)]+chillframe$CH[i]
  }
  
  ###
  # identifies the first day the where the chill requirements are met.
  if(plant==1) chillidx = which(chillframe$ACH >= 1375)
  if(plant==2 | plant==3) chillidx = which(chillframe$ACH >= 1250)
  
  chillDOY = chillframe$DOY[chillidx[1]] # chill DOY is returned in the normal calendar year DOY.
  chillDOY
}

###

# Growing degree hours calculator - supporting function for SI-x FLI and FBI calcs, not called directly.
growdh = function(tmax,tmin,daylen,baset){
  # INPUTS
  # tmax - scalar daily high temperature for a given day
  # tmin - scalar daily low temperature for a given day
  # daylen - scalar daylength for a given day (calculated from daylength function)
  # baset - base temperature for determining growing degree hours. Typically this is 0C, 31F, or 273.15K
  #
  # OUTPUTS
  # GDHOUR - total growing degree hours for a given day. Similar calc to growing degree days.
  
  if(tmin==0) tmin=0.01
  if(tmax==tmin) tmax=tmax+0.01
  
  DTR = tmax-tmin # diurnal temperature range
  IDL = floor(daylen) # round daylength to nearest integer
  
  ###
  # Estimates hourly temperature curves during daylight hours
  Thr = c()
  Thr[1] = tmin
  Thr[2:(IDL+1)]=DTR*sin(pi/(daylen+4)*(1:IDL))+tmin
  
  ###
  # Estimates hourly temperature curves during night hours and strings together with daylight hours.
  TS1 = DTR*sin(pi/(daylen+4)*daylen)+tmin
  if(TS1 <=0) TS1 = 0.01
  Thr[(IDL+2):24] = TS1-(TS1-tmin)/(log(24-daylen))*log((1:(23-IDL)))
  
  ###
  # growing degree hours calculation
  GHRS = Thr-baset # like growing degree days, this uses the actual temperature difference
  GHRS = ifelse(GHRS<0,0,GHRS) # values less than 0 are made zero
  GDHOUR = sum(GHRS)
  GDHOUR
}

###

# Synoptic Event Identifier - supporting function for SI-x FLI and FBI calcs, not called directly.
synval = function(GDH,lagGDH,type="leaf"){
  # INPUTS
  # GDH - scalar growing degree hours for a given day
  # lagGDH - vector growing degree hours for the 7 previous days
  # type - for FLI calcs, type="leaf"
  #        for FBI calcs, type="bloom"
  #
  # OUTPUTS - list with the following
  # synflag - synoptic event flag. Has a value of 1 if a synoptic event happened, value of 0 otherwise.
  # dde2 - sum of growing degree hours for current day and previous two days
  # dd57 - sum of growing degree hours 5-7 prior to current day.
  
  ##
  # Synoptic event definition depends on if one is interested in FLI or FBI
  SYNFLAG=0
  if(type == "leaf") LIMIT=637
  if(type == "bloom") LIMIT=2001
  if(type != "leaf" & type !="bloom") stop("Improper event type",.call=TRUE)
  
  ##
  # Synoptic event occurs if the growing degrees of most recent three days are >= LIMIT
  VALUE=GDH+lagGDH[1]+lagGDH[2]
  if(VALUE >= LIMIT){
    SYNFLAG=1
  } else {
    SNYFLAG=0
  }
  
  ##
  # gathering outputs
  DDE2=GDH+lagGDH[1]+lagGDH[2]
  DD57=sum(lagGDH[5:7])
  output = list(synflag = SYNFLAG,ddE2=DDE2,dd57=DD57)
  output
}

####

# FLI and FBI - primary workhorse for SI-x function, calls multiple supporting functions.
leafindex = function(tasmax,tasmin,daylen,baset,refdate,type="leaf",plant=1, verbose=FALSE){
  # INPUTS
  # tasmax - vector of daily high temperatures for a given year
  # tasmin - vector of daily low temperatures for a given year
  # daylen - vector of daylength for the given year
  # *** tasmax, tasmin, and daylen must be the same length.
  #
  # baset - scalar, base temperature for growing degree hours calculation (typically 0C)
  # refdate - scalar, reference day of year at which to begin calculations of FLI or FBI.
  #           for FLI, refdate = 1
  #           for FBI, refdate = FLI
  # type - for FLI calcs, type="leaf"
  #        for FBI calcs, type="bloom"
  # plant - scalar, signifies plant type 1, 2, or 3. plant = 1 for lilac, plant = 2 for arnold red, plant =3 for zabelli
  # verbose - FALSE by default, set this to TRUE to help with debugging.
  #
  # OUTPUTS
  # OUTDOY - The value of FLI or FBI (depending on the type input), for the specified plant
  
  ###
  # legacy code intialization inherited from Matlab
  ID2 = refdate
  daymax = 240
  Eflag = 0
  AGDH = 0
  SYNOP = 0
  CNT = refdate-1
  lagGDH = rep(0,7)
  
  ###
  # coefficients for estimating FLI and FBI by plant
  Alf=list(c(3.306,13.878,0.201,0.153),c(4.266,20.899,0.000,0.248),c(2.802,21.433,0.266,0.000))
  Alb=list(c(-23.934,0.116),c(-24.825,0.127),c(-11.368,0.096))
  
  ###
  # while loop goes through all days until the requirements for FLI or FBI are met.
  parametersout = matrix(NA,nrow=366,ncol=6)
  lagGDHvals = matrix(NA,nrow=366,ncol=7)
  
  while(CNT<daymax & Eflag==0){
    
    ###
    # growing degree hours calculation
    CNT = CNT+1
    #message("DOY ",CNT," - doing calcs for MDSUM1")
    GDH = growdh(tasmax[CNT],tasmin[CNT],daylen[CNT],baset)
    
    ###
    # initializing lagGDH with faux values for the first day
    if(CNT==1 & type=="leaf"){
      lagGDH[1] = GDH
      lagGDH[2] = GDH
    }
    
    lagGDHvals[CNT,] = lagGDH
    
    ###
    # Identify if synoptic event occurred
    synlist = synval(GDH,lagGDH,type=type)
    
    ###
    # start accumulating synoptic events and growing degree hours once past the reference date
    if(CNT>=refdate){
      AGDH = GDH+AGDH
      if(synlist[[1]]==1) SYNOP=SYNOP+1
    }
    
    ###
    # calculates the indicator (MDSUM1) for if first leaf or first bloom is achieved
    # functions inherited from Matlab code translated to R
    if(CNT>=refdate){
      MDS0 = CNT-refdate # number of days past the reference date
      parametersout[CNT,1] = MDS0
      parametersout[CNT,2] = SYNOP
      parametersout[CNT,3] = synlist[[2]]
      parametersout[CNT,4] = synlist[[3]]
      parametersout[CNT,5] = AGDH
      parametersout[CNT,6] = GDH
      if(type=="leaf"){
        if(plant==1) MDSUM1=(Alf[[1]][1]*MDS0)+(Alf[[1]][2]*SYNOP)+(Alf[[1]][3]*synlist[[2]])+(Alf[[1]][4]*synlist[[3]])
        if(plant==2) MDSUM1=(Alf[[2]][1]*MDS0)+(Alf[[2]][2]*SYNOP)+(Alf[[2]][3]*synlist[[2]])+(Alf[[2]][4]*synlist[[3]])
        if(plant==3) MDSUM1=(Alf[[3]][1]*MDS0)+(Alf[[3]][2]*SYNOP)+(Alf[[3]][3]*synlist[[2]])+(Alf[[3]][4]*synlist[[3]])
      }
      if(type=="bloom"){
        if(plant==1) MDSUM1=(Alb[[1]][1]*MDS0)+(Alb[[1]][2]*AGDH)
        if(plant==2) MDSUM1=(Alb[[2]][1]*MDS0)+(Alb[[2]][2]*AGDH)
        if(plant==3) MDSUM1=(Alb[[3]][1]*MDS0)+(Alb[[3]][2]*AGDH)
      }
    } else {
      MDSUM1 = 1
    }
    #message("MDSUM1 = ",MDSUM1)
    
    ###
    # First leaf or first bloom has been achieved when MDSUM1 >= 1000
    # the if statements below trigger the end of the while loop and determine OUTDOY.
    if(MDSUM1>=999.5 & Eflag==0){
      if(type=="leaf"){
        if(plant==1 | plant==3) OUTDOY = CNT
        if(plant==2) OUTDOY = CNT+1
      }
      if(type=="bloom"){
        OUTDOY = CNT
      }
      Eflag=1
    }
    
    ###
    # holding on to up to 7 previous day so of growing degree hours
    lagGDH[2:7]=lagGDH[1:6]
    lagGDH[1] = GDH
  }
  if(CNT>=daymax){
    OUTDOY=NA
  }
  if(verbose==TRUE){
    list(OUTDOY=round(OUTDOY),parameters=parametersout,lagGDHvals = lagGDHvals) # this is FLI or FBI depending on the type argument in the function call.
  } else {
    round(OUTDOY)
  }
}

####

# Freeze information calculator - supporting function, called directly and output used in other functions.
freezedates = function(tasmin,frzval,DOY){
  # INPUTS
  # tasmin - vector of daily low temperatures for a given year
  # DOY - vector of day of year values for the given year
  # *** tasmin and DOY must be the same length.
  #
  # frzval - scalar, freeze temperature of interest (hard freeze is tasmin < -2.2C)
  #
  # OUTPUTS - list including the following
  # firstfreeze - scalar, first day of year a hard freeze occurs in the fall
  # lastfreeze - scalar, last day of year a hard freeze occurs in the spring
  # freeseperiod - scalar, range of days between firstfreeze and lastfreeze
  # freezedata - data.frame, table with tasmin, DOY, and adjusted DOYs for all hard freeze days
  
  ###
  # initialize data table, adjust calendar, re-order to July 1 - June 30 calendar
  frzframe = data.frame(tasmin,DOY)
  if(nrow(frzframe)==365){
    frzframe$DOYadj = ifelse(DOY<=180,DOY+185,DOY-180)
  } else {
    frzframe$DOYadj = ifelse(DOY<=180,DOY+186,DOY-180)
  }
  frzframe = frzframe[order(frzframe$DOYadj),]
  
  ###
  # identify which rows of the table are hard freeze days
  idx = which(frzframe$tasmin<=frzval)
  
  ###
  # determine firstfreeze, lastfreeze, and freeze period, and write outputs list
  lastfreeze = frzframe[idx[length(idx)],2]
  firstfreeze = frzframe[idx[1],2]
  freezeperiod = frzframe[idx[length(idx)],3]- frzframe[idx[1],3]
  if(length(lastfreeze)==0) lastfreeze=NA
  freezelist = list(firstfreeze=firstfreeze,lastfreeze=lastfreeze,freezeperiod=freezeperiod,freezedata=frzframe[idx,])
  freezelist
}

####

# Damage Index calculator - SI-x index, called directly
damageindex = function(leafidx,lastfreeze){
  # INPUTS
  # leafidx - scalar, the FLI value
  # *** in the literature FLI is used, but FBI could be input here also.
  # lastfreeze - scalar, the day of year for the last freeze
  #
  # OUTPUTS
  # DMG - damage index, increasing values indicate greater damage potential
  
  idx = lastfreeze - leafidx # damage index calculation
  DMG = ifelse(idx<0 | idx>=180,0,idx) # setting values where last freeze happens before first leaf to zero, i.e. no damage potential
  # Also, if the last freeze occurs in the winter following first leaf, this is set to zero. i.e. no damage potential that spring.
  DMG
}

####

# False Spring Indicator - based on Allstadt et al. (2015)
falsesprings = function(SI=c(60,65),freezedata){
  # INPUTS
  # SI - vector with two arguments, the FLI and the FBI
  # *** FLI should be first.
  # freezedata - data.frame of freeze data from the freezedates function
  #
  # OUTPUTS - FSlist - list with the following
  # EFS - Early False Spring Indicator
  #       EFS occurs when a hard freeze happens 7 or more days after FLI and before FBI
  #       EFS = 1 when an early false spring occurs, EFS = 0 otherwise
  # LFS - Late False Spring Indicator
  #       LFS occurs when a hard freeze happens any day after FBI
  #       LFS = 1 when a late false spring occurs, LFS = 0 otherwise
  
  ###
  # determines the number of days where EFS or LFS conditions are met
  EFSidx = which(freezedata$DOY>=(SI[1]+7) & freezedata$DOY<=SI[2])
  LFSidx = which(freezedata$DOY>SI[2])
  
  ###
  # Checks length of the above and sets EFS or LFS to 1 if the length is greater than zero.
  EFS = ifelse(length(EFSidx)>0,1,0)
  LFS = ifelse(length(LFSidx)>0,1,0)
  
  ###
  # gathering output
  FSlist= list(EFS=EFS,LFS=LFS)
  FSlist
}

####

# False Spring Exposure Index - from Peterson and Abatzoglou, 2014
FSEI = function(falsespringvector){
  # INPUTS
  # falsespringvector - vector with N years of the EFS or LFS from false springs
  #
  # OUTPUTS
  # FSEI - scalar, False Spring Exposure Index value.
  
  FSEI = (sum(falsespringvector,na.rm=TRUE)/length(falsespringvector))*100
  FSEI
}
