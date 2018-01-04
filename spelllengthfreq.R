#y=1
#yearidx = which(years==yrs[y])
#inputdata= prens1[100,65,yearidx]
#maskfrominput = ifelse(inputdata>=0.254,1,0)
#maskdata= prens_mask1[100,65,yearidx]
#spell_length=3
#thres=0.254
#outtype="count"
#type="dry"
#inputdata=maskdata

spell_length_calc = function(inputdata,premasked=FALSE,type="lte",spell_length=3,thres=0.254,outtype="ratio"){

  #inputdata=prfut1ens3[44,6,yearidx]
  if(premasked==FALSE){
    if(type=="lt") eventdays = which(inputdata < thres)
    if(type=="lte") eventdays = which(inputdata <= thres)
    if(type=="gte") eventdays = which(inputdata>= thres)
    if(type=="gt") eventdays = which(inputdata> thres)
  } else {
    eventdays=which(inputdata== 1)
  }
  
  eventdays = c(eventdays,max(eventdays)+2)
  eventdiff = diff(eventdays)
  
  if(length(eventdays)!=length(inputdata) & any(eventdiff!=1)==TRUE){
    eventbreaks = which(eventdiff>1)
    event_lengths = c(eventbreaks[1],eventbreaks[2:length(eventbreaks)]-eventbreaks[1:(length(eventbreaks)-1)])
    if(outtype=="ratio") spellfreq = length(which(event_lengths>=spell_length))/length(event_lengths)
    if(outtype=="count") spellfreq = length(which(event_lengths>=spell_length))
  } else {
    if(outtype=="ratio") spellfreq = 1
    if(outtype=="count") spellfreq = round(length(eventdays)/spell_length)
  }
  spellfreq
}