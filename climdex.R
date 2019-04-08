#####
# rolling sums calculator

calcrollsum = function(inputdata,size=5){
  require(zoo)
  rollsum = rollapply(inputdata,size,sum,na.rm=TRUE)
  max(rollsum,na.rm=TRUE)
}

#####
# day of last freeze event

lastfreeze = function(inputdata,inputtimes,startdate,threshold=273.15){
  #threshold in K
  #inputdata should be one year of daily data
  
  #inputdata = tempdata[1,69,]
  #inputtimes= times[yearidx]
  #startdate = startdate
  #threshold = 273.15
  #inputdata=tempdata[r,c,];startdate=startdate;inputtimes=times[yearidx]
  if(length(inputdata)<365) stop("This is not a full year of data!",.call=TRUE)
  
  if(length(inputdata)==366){
    inputdata=inputdata[-60]
    inputtimes = inputtimes[-60]
    #message("Removed leap day")
  } 
  dates = as.character(substr(seq(as.Date("1999-01-01"),as.Date("1999-06-30"),by="day"),6,10))
  endpoint = length(dates)
  doy = 1:endpoint
  doyframe = data.frame(doy,dates)
  #message("Made doy reference table")
  
  useinput = inputdata[1:endpoint]
  #message("determined input to use")
  if(all(is.na(useinput)==FALSE)==TRUE){
   # message("data is available checking for last freeze")
  freezeidx = which(useinput<=threshold)
    if(length(freezeidx)>0){
      #message("there are freeze days - checking for last freeze")
      #lastfreezeday = inputtimes[freezeidx[length(freezeidx)]]
      #message("lastfreezeday is ",lastfreezeday)
      #date_of_last_freeze = as.Date(startdate)+lastfreezeday
      #message("date of last freeze is ",date_of_last_freeze)
      #doy_of_last_freeze = as.numeric(doyframe[which(doyframe[,2]==substr(date_of_last_freeze,6,10)),1])
      doy_of_last_freeze = freezeidx[length(freezeidx)]
      #message("doy for last freeze is ",doy_of_last_freeze)
    } else {
      #message("No freeze happened, printing 0")
      doy_of_last_freeze = 0
    }
  } else {
    #message("no data, printing NA")
    doy_of_last_freeze=NA
  }
  doy_of_last_freeze
  
}

########
# growing season length

GSLcalc = function(tmaxdata,tmindata,inputtimes,startdate){
  require(zoo)
  #threshold in K, temperatures in K
  #tmindata and tmaxdata should be one year of daily data
  if(length(tmindata)<365 | length(tmaxdata)<365) stop("This is not a full year of data! Check on tmindata or tmaxdata",.call=TRUE)
  
  #tmaxdata = tempdata1[1,69,]
  #tmindata = tempdata2[1,69,]
  #inputtimes = times[yearidx]
  #startdate = startdate
  
  inputdata = (tmaxdata+tmindata)/2
  
  if(all(is.na(inputdata)==FALSE)==TRUE){
  if(length(inputdata)==366){
    inputdata=inputdata[-60]
    inputtimes = inputtimes[-60]
  } 
  
    midpoint = 181
    dates = substr(seq(as.Date("1999-01-01"),as.Date("1999-12-31"),by="day"),6,10)
  
    doy = 1:365
    doyframe = data.frame(doy,dates)
    
    events_e = rep(0,365)
    events_l = rep(0,365)
    
    eventdays_early = which(inputdata[1:midpoint] >= 278.15)
    eventdays_late = which(inputdata[(midpoint+1):365] < 278.15)
  
    events_e[eventdays_early]=1
    events_l[eventdays_late+181]=1
    
    early6 = rollapply(events_e,6,sum)
    earidx = which(early6==6)
    if(length(earidx)>0){
    if(length(earidx)>1) earidx = earidx[1]
    timeear = inputtimes[earidx]+2.5
    } else {
      timeear = floor(inputtimes[1])
    }
    
    late6 = rollapply(events_l,6,sum)
    lateidx = which(late6==6)
    if(length(lateidx)>0){
    if(length(lateidx)>1) lateidx = lateidx[1]
    timelate = inputtimes[lateidx]+2.5
    } else {
      timelate = floor(inputtimes[length(inputtimes)])
    }
    GSL = timelate - timeear # Growing season length in days
  } else {
    GSL=NA
  }
  GSL
}

#############
# spell counter

spell_length_calc = function(data_vec,premasked,cond,spell_len=3,thold=0.254,outtype="count"){
  
  if(all(is.na(data_vec)==FALSE)==TRUE){
  if(premasked==TRUE & outtype=="ratio") stop("Ratio cannot be calculated on premasked data!",call.=TRUE)
  
  if(premasked==FALSE){
    if(cond=="LT") eventdays = which(data_vec < thold)
    if(cond=="LE") eventdays = which(data_vec <= thold)
    if(cond=="GE") eventdays = which(data_vec >= thold)
    if(cond=="GT") eventdays = which(data_vec > thold)
  } else {
    eventdays=which(data_vec== 1)
  }
  
  eventdays = c(eventdays,max(eventdays)+2)
  eventdiff = diff(eventdays)
  
  if(length(eventdays)!=length(data_vec) & any(eventdiff!=1)==TRUE){
    eventbreaks = which(eventdiff>1)
    event_lengths = c(eventbreaks[1],eventbreaks[2:length(eventbreaks)]-eventbreaks[1:(length(eventbreaks)-1)])
    if(outtype=="ratio") spellfreq = length(which(event_lengths>=spell_len))/length(event_lengths)
    if(outtype=="count") spellfreq = length(which(event_lengths>=spell_len))
    if(outtype=="max") spellfreq = max(event_lengths,na.rm=TRUE)
    if(outtype=="mean") spellfreq = mean(event_lengths,na.rm=TRUE)
  } else {
    if(outtype=="ratio") spellfreq = 1
    if(outtype=="count") spellfreq = round(length(eventdays)/spell_len)
    if(outtype=="max") spellfreq = length(data_vec)
    if(outtype=="mean") spellfreq = length(data_vec)/2
  }
  return(spellfreq)
  } else {
    spellfreq = NA
    return(spellfreq)
  }
}

###########################
# heatwaves

heatwave.calc = function(tmaxdata,tmindata,tmaxq95,tminq95){
  
  if(length(tmindata)<365 | length(tmaxdata)<365) stop("This is not a full year of data! Check on tmindata or tmaxdata",.call=TRUE)
  if(length(tmindata)!=length(tmaxdata)) stop("The length of tmaxdata isn't the same as the length of tmindata",.call=TRUE)
  tmaxmask = ifelse(tmaxdata>=tmaxq95,1,0)
  tminmask = ifelse(tmindata>=tminq95,1,0)
  inputdata = ifelse(tmaxmask==1 & tminmask==1,1,0)
  
  rm(tmaxmask)
  rm(tminmask)
  if(length(inputdata)==366){
    inputdata=inputdata[-60]
  } 
  
  if(all(is.na(inputdata)==FALSE)==TRUE){
   
    hwtotal = spell_length_calc(inputdata,premasked=FALSE,cond="GE",spell_len=3,thold=1,outtype="count")
  } else {
    hwtotal = NA
  }
  hwtotal
}
