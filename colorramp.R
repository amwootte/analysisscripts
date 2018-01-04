
colorramp = function(inputdata,colorchoice,Blimit,type = "difference",use_fixed_scale = FALSE, fixed_scale=c(-100,100)){
  
  if(use_fixed_scale==FALSE){
    message("Not using a fixed scale")
    datarange = range(inputdata,na.rm=TRUE)
    datarange[1]=floor(datarange[1])
    if(datarange[1] %% 2 != 0) datarange[1]=datarange[1]-1
    datarange[2]=ceiling(datarange[2])
    if(datarange[2] %% 2 != 0) datarange[2]=datarange[2]+1
  } else {
    message("Using a fixed scale")
    #tmp = strsplit(fixed_scale,",")
    #datarange = c(as.numeric(tmp[[1]][1]),as.numeric(tmp[[1]][2]))
    datarange = fixed_scale
  }
  
  if(datarange[1]>=0 & type=="difference"){centerpoint = 0; startpoint=0; datarange[1]=0; message("type=difference")}
  if(datarange[1]<0 & type=="difference"){centerpoint = 0; startpoint=datarange[1]; message("type=difference"); if(datarange[2]<0){datarange[2]=0}}
  
  if(type=="ratio"){centerpoint=1; startpoint=datarange[1]; message("type=ratio")}
  if(type=="raw"){centerpoint=datarange[1]; startpoint=datarange[1]; message("type=raw")}
  
  breakcheck = 1
  breaklist = c(0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,50,100,200,250,300,500,1000)
  
  actualbins = diff(datarange)/breaklist
  actidx = which(actualbins<Blimit)
  dataact = actualbins[actidx]-floor(actualbins[actidx])
  
  if(any(dataact==0)==TRUE){
    message("exact match for bins")
    dataidx = which(dataact==0)
    breakcheck=actidx[dataidx[1]]
  } else {
    message("no exact match going through while loop")
    checkpoint = any(dataact==0)
    counter=1
    while(checkpoint==FALSE){
      datarange[1] = floor(datarange[1]/(10^counter))*10^counter
      datarange[2] = ceiling(datarange[2]/(10^counter))*10^counter
      actualbins = diff(datarange)/breaklist
      actidx = which(actualbins<Blimit)
      dataact = actualbins[actidx]-floor(actualbins[actidx])
      dataidx = which(dataact==0)
      
      if(length(dataidx)>=1){
        breakcheck=actidx[dataidx[1]]
        checkpoint = any(dataact==0)
        break
      } else {
        counter=counter+1
        checkpoint = any(dataact==0)
      }
      
    }
  }
  
  if(datarange[2]==0 & colorchoice=="redtoblue") colorchoice="redtowhite"
  if(startpoint==0 & centerpoint==0 & colorchoice=="bluetored") colorchoice="whitetored"
  if(startpoint==0 & centerpoint==0 & colorchoice=="browntogreen") colorchoice="whitetogreen"
  
  zlimdiff = datarange
  breaksdiff = c(seq(datarange[1],datarange[2],by=breaklist[breakcheck]))
  
  if(any(breaksdiff==centerpoint)==FALSE & zlimdiff[1]<centerpoint & zlimdiff[2]>centerpoint){
    idx = which(abs(breaksdiff)==min(abs(breaksdiff)))
    if(length(idx)==1) breaksdiff[idx]=centerpoint
    if(length(idx)>1){
      breaksdiff = c(breaksdiff[1:(idx[1]-1)],centerpoint,breaksdiff[(idx[2]+1):length(breaksdiff)])
    }
  } 
  
  message("zlimdiff: ",zlimdiff)
  message("centerpoint: ",centerpoint)
  message("startpoint: ",startpoint)
  
  if(startpoint==centerpoint){
    message("startpoint matches centerpoint")
    if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
    if(colorchoice == "yellowtored") colorbardiff = colorRampPalette(c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))(length(breaksdiff)-1)
    
    if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
  } else {
    if(datarange[2]>centerpoint){
      message("datarange[2] > centerpoint")
      message("colorchoice = ",colorchoice)
      zeroidx = which(breaksdiff==centerpoint)
      if(colorchoice == "bluetored"){
        colorbardiff = c(colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "redtoblue"){
        colorbardiff = c(colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "browntogreen"){
        colorbardiff = c(colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "greentobrown"){
        colorbardiff = c(colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-zeroidx))
      } 
      
    } else {
      message("odd if")
      if(colorchoice == "redtowhite"){
        colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
      } 
      
    }
  }
  
  output = list(zlimdiff,breaksdiff,colorbardiff)
  output
  
}
