###################
#
# Revised HS09 Component Calculator Functions

########
# combine data

combinedata = function(arrayvals,idx){
arrayout = array(arrayvals,dim=c(dim(arrayvals)[1],dim(arrayvals)[2],dim(arrayvals)[3]*dim(arrayvals)[4]))
arrayout
}

########
# Internal Variability Calculator -- specialized version, time series version

intvar_spts = function(resvaltable,weights=list(weightsGCM,weightsDS),GCMs,DSs,files=files){
  
  #resvaltable = resvals
  #weights = list(GCMweights,DSweights)
  #files = memnames
  
  varsout = c()
  
  for(i in 1:length(GCMs)){
    
    varsoutDS = c()
    weightsDSused = c()
    for(j in 1:length(DSs)){
      #message("I'm on ",i)
      
      idx1 = grep(GCMs[i],files)
      idx2 = grep(DSs[j],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      IDX = idx+1
      
      if(length(idx)>0){
        
        if(length(idx)==1){
          
          #test = apply(resvalarray[,IDX],c(1,2),var,na.rm=TRUE)
          #varsoutDS[,,j] = weights[[2]][,,j] * test
          varsoutDS[j] = var(resvaltable[,IDX],na.rm=TRUE)
          
        } else {
          arrayout = c()
          for(k in 1:length(IDX)){
            arrayout = c(arrayout,resvaltable[,IDX[k]])
          }
          #test = apply(arrayout,c(1,2),var,na.rm=TRUE)
          #varsoutDS[,,j] = weights[[2]][,,j] * test
          varsoutDS[j] = var(arrayout,na.rm=TRUE)
        }
        
        weightsDSused[j]=weights[[2]][j]
        
      } 
      
    }
    
    for(c in 1:length(varsoutDS)) varsoutDS[c] = varsoutDS[c]*(weightsDSused[c]/sum(weightsDSused,c(1,2),na.rm=TRUE))
    
    varsout[i] = sum(varsoutDS,na.rm=TRUE) * weights[[1]][i]
    
  } 
  
  V = sum(varsout,na.rm=TRUE)
  V = ifelse(V==0,NA,V)  
  
  V
}

########
# Internal Variability Calculator -- specialized version

intvar_sp = function(resvalarray,weights=list(weightsGCM,weightsDS),GCMs,DSs,files=files){
  
  varsout = array(NA,dim=c(dim(resvalarray)[1],dim(resvalarray)[2],length(GCMs)))
  
  for(i in 1:length(GCMs)){
    
    varsoutDS = array(NA,dim=c(dim(resvalarray)[1],dim(resvalarray)[2],length(DSs)))
    weightsDSused = array(NA,dim=c(dim(resvalarray)[1],dim(resvalarray)[2],length(DSs)))
    for(j in 1:length(DSs)){
    #message("I'm on ",i)
    
    idx1 = grep(GCMs[i],files)
    idx2 = grep(DSs[j],files)
      
    idx = idx2[which(idx2 %in% idx1)]

    if(length(idx)>0){

      if(length(idx)==1){

	test = apply(resvalarray[,,,idx],c(1,2),var,na.rm=TRUE)
	#varsoutDS[,,j] = weights[[2]][,,j] * test
	varsoutDS[,,j] = test

      } else {
	      arrayout = combinedata(resvalarray[,,,idx],idx)        
	      test = apply(arrayout,c(1,2),var,na.rm=TRUE)
	      #varsoutDS[,,j] = weights[[2]][,,j] * test
		varsoutDS[,,j] = test
      }
      
	weightsDSused[,,j]=weights[[2]][,,j]

    } 
  
    }

  for(c in 1:dim(varsoutDS)[3]) varsoutDS[,,c] = varsoutDS[,,c]*(weightsDSused[,,c]/apply(weightsDSused,c(1,2),sum,na.rm=TRUE))

  varsout[,,i] = apply(varsoutDS,c(1,2),sum,na.rm=TRUE) * weights[[1]][,,i]
  
  } 

  V = apply(varsout,c(1,2),sum,na.rm=TRUE)
  V = ifelse(V==0,NA,V)  

V
}

########
# Internal Variability Calculator -- Generalized version

intvar = function(resvalarray,weights=list(weightsGCM,weightsDS),GCMs,DSs,files=files){
  
  varsout = array(NA,dim=c(dim(resvalarray)[1],dim(resvalarray)[2],length(GCMs)))
  
  for(i in 1:length(GCMs)){
    
    varsoutDS = array(NA,dim=c(dim(resvalarray)[1],dim(resvalarray)[2],length(DSs)))
    for(j in 1:length(DSs)){
    #message("I'm on ",i)
    
    idx1 = grep(GCMs[i],files)
    idx2 = grep(DSs[j],files)
      
    idx = idx2[which(idx2 %in% idx1)]

    if(length(idx)>0){

      if(length(idx)==1){

	test = apply(resvalarray[,,,idx],c(1,2),var,na.rm=TRUE)
	varsoutDS[,,j] = weights[[2]][,,j] * test

      } else {
	      arrayout = combinedata(resvalarray[,,,idx],idx)        
	      test = apply(arrayout,c(1,2),var,na.rm=TRUE)
	      varsoutDS[,,j] = weights[[2]][,,j] * test
      }
      
    } 
  
    }
  
  varsout[,,i] = apply(varsoutDS,c(1,2),sum,na.rm=TRUE) * weights[[1]][,,i]
  
  } 

  V = apply(varsout,c(1,2),sum,na.rm=TRUE)
  V = ifelse(V==0,NA,V)  

V
}

########
# GCM Uncertainty Calculator -- Specified version, time series

GCMuncts = function(fittedvaltable,weights=GCMweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  #fittedvaltable = fittedvals
  #weights = GCMweights
  #files = memnames
  
  wvartable = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(scens))
  
  for(s in 1:length(scens)){
    
    wvartemp = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(DSs))
    
    for(m in 1:length(DSs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(DSs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      IDX = idx+1
      
      wvar = c()
      
      if(length(idx)>1){
        namelist = do.call(rbind,strsplit(files[idx],"_"))
        
        #print(which(GCMs %in% namelist[,2]))
        
        W = weights[which(GCMs %in% namelist[,1])]
        checkthis = which(GCMs %in% namelist[,1])
        Wsum = sum(W,na.rm=TRUE)
        for(v in 1:length(checkthis)) W[v]=W[v]/Wsum
        
        #print(W[18,12,])
        
        for(i in 1:nrow(fittedvaltable)){
          if(all(is.na(fittedvaltable[i,IDX])==TRUE)==FALSE){
            wvar[i] = weighted.varts(fittedvaltable[i,IDX],w=W,na.rm=TRUE)
          }
          #message("I'm on point ",i)	
        }
      } 
      
      wvartemp[,m] = wvar
      #message("I'm on DS ",m)        
    }
    
    wvartemp2 = apply(wvartemp,1,sum,na.rm=TRUE)
    wvartable[,s] = wvartemp2
    #message("I'm on scenario ",s)
  }
  
  wvarsum = apply(wvartable,1,sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  
  
  M = wvarsum * (1/length(scens)) * (1/length(DSs))
  M
}


########
# GCM Uncertainty Calculator -- Specified version

GCMunc = function(fittedvalarray,weights=GCMweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  wvartable = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
    
    wvartemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(DSs)))

    for(m in 1:length(DSs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(DSs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      
      wvar = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))
      
      #message("m is ",s)
      #message("m is ",m)
      #message("length(idx)=",length(idx))

      if(length(idx)>1){
      namelist = do.call(rbind,strsplit(files[idx],"_"))

      #print(which(GCMs %in% namelist[,2]))

      W = weights[,,which(GCMs %in% namelist[,2])]
	checkthis = which(GCMs %in% namelist[,2])
	Wsum = apply(W,c(1,2),sum,na.rm=TRUE)
	for(v in 1:length(checkthis)) W[,,v]=W[,,v]/Wsum
   
     #print(W[18,12,])

      for(i in 1:dim(fittedvalarray)[3]){
 	 wvar[,,i] = weighted.var(fittedvalarray[,,i,idx],w=W,na.rm=TRUE)
	#message("I'm on point ",i)	
	}
      } 
      
     wvartemp[,,,m] = wvar
     #message("I'm on DS ",m)        
    }
    
    wvartemp2 = apply(wvartemp,c(1,2,3),sum,na.rm=TRUE)
    wvartable[,,,s] = wvartemp2
    #message("I'm on scenario ",s)
  }
  
  wvarsum = apply(wvartable,c(1,2,3),sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  

  M = wvarsum * (1/length(scens)) * (1/length(DSs))
  M
}

########
# GCM Uncertainty Calculator -- Generalized version

GCMuncold = function(fittedvalarray,weights=GCMweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  wvartable = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
    
    wvartemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(DSs)))

    for(m in 1:length(DSs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(DSs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      
      wvar = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))
      
      #message("m is ",s)
      #message("m is ",m)
      #message("length(idx)=",length(idx))

      if(length(idx)>1){
      namelist = do.call(rbind,strsplit(files[idx],"_"))

      #print(which(GCMs %in% namelist[,2]))

      W = weights[,,which(GCMs %in% namelist[,2])]
      
      #print(W[18,12,])

      for(i in 1:dim(fittedvalarray)[3]){
 	 wvar[,,i] = weighted.var(fittedvalarray[,,i,idx],w=W,na.rm=TRUE)
	#message("I'm on point ",i)	
	}
      } 
      
     wvartemp[,,,m] = wvar
     #message("I'm on DS ",m)        
    }
    
    wvartemp2 = apply(wvartemp,c(1,2,3),sum,na.rm=TRUE)
    wvartable[,,,s] = wvartemp2
    #message("I'm on scenario ",s)
  }
  
  wvarsum = apply(wvartable,c(1,2,3),sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  

  M = wvarsum * (1/length(scens)) * (1/length(DSs))
  M
}

########
# DS Uncertainty Calculator -- Specialized version, timeseries

DSuncts = function(fittedvaltable,weights=DSweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  #fittedvaltable = fittedvals
  #weights = DSweights
  #files = memnames
  
  wvartable = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(scens))
  
  for(s in 1:length(scens)){
    
    wvartemp = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(GCMs))
    
    for(m in 1:length(GCMs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(GCMs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      IDX = idx+1
      
      wvar = c()
      
      #message("m is ",s)
      #message("m is ",m)
      #message("length(idx)=",length(idx))
      
      if(length(idx)>1){
        namelist = do.call(rbind,strsplit(files[idx],"_"))
        #namelist = files[idx]
        #print(which(GCMs %in% namelist[,2]))
        
        W = weights[which(DSs %in% paste(namelist[,2],namelist[,3],sep="_"))]
        #W = weights[,,which(DSs %in% namelist)]
        checkthis = which(DSs %in% paste(namelist[,2],namelist[,3],sep="_"))
        #checkthis = which(DSs %in% namelist)
        Wsum = sum(W,na.rm=TRUE)
        for(v in 1:length(checkthis)) W[v]=W[v]/Wsum
        
        #print(W[18,12,])
        
        for(i in 1:nrow(fittedvaltable)){
          if(all(is.na(fittedvaltable[i,IDX])==TRUE)==FALSE){
            wvar[i] = weighted.varts(fittedvaltable[i,IDX],w=W,na.rm=TRUE)
          }
        }
      } 
      
      wvartemp[,m] = wvar
      
    }
    
    wvartemp2 = apply(wvartemp,1,sum,na.rm=TRUE)
    wvartable[,s] = wvartemp2
    
  }
  
  wvarsum = apply(wvartable,1,sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  
  
  D = wvarsum * (1/length(scens)) * (1/length(GCMs))
  D
}



########
# DS Uncertainty Calculator -- Specialized version

DSunc = function(fittedvalarray,weights=DSweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  wvartable = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
    
    wvartemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(GCMs)))

    for(m in 1:length(GCMs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(GCMs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      
      wvar = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))
      
      #message("m is ",s)
      #message("m is ",m)
      #message("length(idx)=",length(idx))

      if(length(idx)>1){
      namelist = do.call(rbind,strsplit(files[idx],"_"))
        #namelist = files[idx]
      #print(which(GCMs %in% namelist[,2]))

      W = weights[,,which(DSs %in% paste(namelist[,3],namelist[,4],sep="_"))]
        #W = weights[,,which(DSs %in% namelist)]
	  checkthis = which(DSs %in% paste(namelist[,3],namelist[,4],sep="_"))
        #checkthis = which(DSs %in% namelist)
	Wsum = apply(W,c(1,2),sum,na.rm=TRUE)
	for(v in 1:length(checkthis)) W[,,v]=W[,,v]/Wsum

      #print(W[18,12,])

      for(i in 1:dim(fittedvalarray)[3]){
 	 wvar[,,i] = weighted.var(fittedvalarray[,,i,idx],w=W,na.rm=TRUE)
	}
      } 
      
     wvartemp[,,,m] = wvar
            
    }
    
    wvartemp2 = apply(wvartemp,c(1,2,3),sum,na.rm=TRUE)
    wvartable[,,,s] = wvartemp2
    
  }
  
  wvarsum = apply(wvartable,c(1,2,3),sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  

  D = wvarsum * (1/length(scens)) * (1/length(GCMs))
  D
}

########
# DS Uncertainty Calculator -- Generalized version

DSuncold = function(fittedvalarray,weights=DSweights,scens=scens,DSs=DSs,GCMs=GCMs,files=files){
  
  wvartable = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
    
    wvartemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(GCMs)))

    for(m in 1:length(GCMs)){
      
      idx1 = grep(scens[s],files)
      idx2 = grep(GCMs[m],files)
      
      idx = idx2[which(idx2 %in% idx1)]
      
      wvar = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))
      
      #message("m is ",s)
      #message("m is ",m)
      #message("length(idx)=",length(idx))

      if(length(idx)>1){
      namelist = do.call(rbind,strsplit(files[idx],"_"))

      #print(which(GCMs %in% namelist[,2]))

      W = weights[,,which(DSs %in% namelist[,3])]
      
      #print(W[18,12,])

      for(i in 1:dim(fittedvalarray)[3]){
 	 wvar[,,i] = weighted.var(fittedvalarray[,,i,idx],w=W,na.rm=TRUE)
	}
      } 
      
     wvartemp[,,,m] = wvar
            
    }
    
    wvartemp2 = apply(wvartemp,c(1,2,3),sum,na.rm=TRUE)
    wvartable[,,,s] = wvartemp2
    
  }
  
  wvarsum = apply(wvartable,c(1,2,3),sum,na.rm=TRUE)
  wvarsum = ifelse(wvarsum==0,NA,wvarsum)  

  D = wvarsum * (1/length(scens)) * (1/length(GCMs))
  D
}

######
# weighted variance, time series
weighted.varts = function(x,w,na.rm=TRUE){
  
  x= as.vector(x)
  
  NAidx = which(is.na(x)==TRUE) # remove points with no data
  
  if(length(NAidx)>=1){
    x = x[-NAidx]
    w = w[-NAidx]
  }
  
  sum.w = sum(w,na.rm=TRUE)
  sum.w2 = sum(w^2,na.rm=TRUE)
  
  fit1 = x * w
  mean.w = sum(fit1,na.rm=TRUE)/sum.w
  
  part1 = (sum.w / (sum.w^2 - sum.w2))
  part2a = x
  for(i in 1:length(x)) part2a[i] = (x[i]-mean.w)^2
  part2b = w*part2a
  
  part2 = sum(part2b,na.rm=TRUE)
  
  output = part1*part2
  
  #for(r in 1:dim(x)[1]){ # check to make sure there's enough data for this variance calculation
  #for(c in 1:dim(x)[2]){
  #output[r,c] = ifelse(length(which(is.na(x[r,c,])==FALSE))<=1,NA,output[r,c])
  #}
  #}	
  
  output
}


########
# Weighted Variance

weighted.var = function(x,w,na.rm=TRUE){

      NAidx = which(is.na(apply(x,3,mean,na.rm=TRUE))==TRUE) # remove points with no data
      
      if(length(NAidx)>=1){
      x = x[,,-NAidx]
      w = w[,,-NAidx]
      }

	if(length(dim(x))>2 & dim(x)[3]!=0){

      sum.w = apply(w,c(1,2),sum,na.rm=TRUE)
      sum.w2 = apply(apply(w,c(1,2),"^",2),c(2,3),sum,na.rm=TRUE)
      
      fit1 = x * w
      mean.w = apply(fit1,c(1,2),sum,na.rm=TRUE)/sum.w
      
      part1 = (sum.w / (sum.w^2 - sum.w2))
      part2a = x
	for(i in 1:dim(x)[3]) part2a[,,i] = (x[,,i]-mean.w)^2
      part2b = w*part2a

      part2 = apply(part2b,c(1,2),sum,na.rm=TRUE)

      output = part1*part2
	
	#for(r in 1:dim(x)[1]){ # check to make sure there's enough data for this variance calculation
	#for(c in 1:dim(x)[2]){
	#output[r,c] = ifelse(length(which(is.na(x[r,c,])==FALSE))<=1,NA,output[r,c])
	#}
	#}	

	} else {
	output = matrix(NA,nrow=nrow(w),ncol=ncol(w))
	}

	output
}

######################
# Scenario Uncertainty Calculator - specified version

scenuncts = function(fittedvaltable,scens,DSs,GCMs,weightsGCM=GCMweights,weightsDS=DSweights,files=files){
  
  #fittedvaltable = fittedvals
  #weightsGCM=GCMweights
  #weightsDS = DSweights
  #files=memnames
  
  #scens = scens
  #GCMs = GCMs
  
  ENSout = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(scens))
  
  for(s in 1:length(scens)){
    
    ENStemp = matrix(NA,nrow=nrow(fittedvaltable),ncol=length(GCMs))
    
    for(g in 1:length(GCMs)){
      
      scenidx = grep(scens[s],files)
      idx1 = grep(GCMs[g],files)
      
      idx = idx1[which(idx1 %in% scenidx)] # finds what matches for GCM and scenario
      IDX = idx+1
      
      WSdtemp = matrix(NA,nrow=nrow(fittedvaltable),length(idx))
      WSgtemp = c()
      
      if(length(idx)!=0){
        if(length(idx)>1 & length(idx)!=0){
          
          filesplit = do.call(rbind,strsplit(files[idx],"_"))
          
          for(i in 1:length(idx)){
            #Wdnorm = (weightsDS[,,which(DSs %in% filesplit[i,3])]/apply(weightsDS[,,which(DSs %in% filesplit[,3])],c(1,2),sum))
            #
            Wdnorm = (weightsDS[which(DSs %in% paste(filesplit[i,2],filesplit[i,3],sep="_"))]/sum(weightsDS[which(DSs %in% paste(filesplit[,2],filesplit[,3],sep="_"))]))
            for(t in 1:nrow(fittedvaltable)) WSdtemp[t,i] = fittedvaltable[t,IDX[i]]*Wdnorm         
          }
          
          WSgtemp = apply(WSdtemp,1,sum,na.rm=TRUE)
          #for(i in 1:dim(fittedvalarray)[3]) WSgtemp[,,i] = WSgtemp[,,i]*weightsGCM[,,g]
          
        } else {
          
          for(i in 1:nrow(fittedvaltable)) WSdtemp[i,1] = fittedvaltable[i,IDX]
          for(i in 1:nrow(fittedvaltable)) WSgtemp[i] = WSdtemp[i,1]
          
        }
      }
      
      ENStemp[,g] = WSgtemp 
    }
    
    GCMcheck = c()
    for(c in 1:dim(ENStemp)[2]) GCMcheck[c] = ifelse(all(is.na(ENStemp[,c])==TRUE)==TRUE,0,1)
    
    GCMidx=which(GCMcheck!=0)
    ENStempuse = ENStemp[,GCMidx]
    
    for(n in 1:length(GCMidx)){
      Wgnorm = (weightsGCM[GCMidx[n]]/sum(weightsGCM[GCMidx],na.rm=TRUE))
      for(i in 1:nrow(fittedvaltable)) ENStempuse[i,n] = ENStempuse[i,n]*Wgnorm
    }
    
    ENSout[,s] = apply(ENStempuse,1,sum,na.rm=TRUE) 
    
  }
  
  S = apply(ENSout,1,var,na.rm=TRUE)
  S = ifelse(S==0,NA,S)
  S
}


######################
# Scenario Uncertainty Calculator - specified version

scenunc = function(fittedvalarray,scens,DSs,GCMs,weightsGCM=GCMweights,weightsDS=DSweights,files=files){
  
  scens = scens
  GCMs = GCMs
  
  ENSout = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
  
    ENStemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(GCMs)))

    for(g in 1:length(GCMs)){
    
        scenidx = grep(scens[s],files)
        idx1 = grep(GCMs[g],files)
        
        idx = idx1[which(idx1 %in% scenidx)] # finds what matches for GCM and scenario
      
        WSdtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(idx)))
        WSgtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))

        if(length(idx)!=0){
        if(length(idx)>1 & length(idx)!=0){
        
	filesplit = do.call(rbind,strsplit(files[idx],"_"))

        for(i in 1:length(idx)){
	   #Wdnorm = (weightsDS[,,which(DSs %in% filesplit[i,3])]/apply(weightsDS[,,which(DSs %in% filesplit[,3])],c(1,2),sum))
	   #
          Wdnorm = (weightsDS[,,which(DSs %in% paste(filesplit[i,3],filesplit[i,4],sep="_"))]/apply(weightsDS[,,which(DSs %in% paste(filesplit[,3],filesplit[,4],sep="_"))],c(1,2),sum))
           for(t in 1:dim(fittedvalarray)[3]) WSdtemp[,,t,i] = fittedvalarray[,,t,idx[i]]*Wdnorm         
        }
        
        WSgtemp = apply(WSdtemp,c(1,2,3),sum)
        #for(i in 1:dim(fittedvalarray)[3]) WSgtemp[,,i] = WSgtemp[,,i]*weightsGCM[,,g]
        
        } else {
          
          for(i in 1:dim(fittedvalarray)[3]) WSdtemp[,,i,1] = fittedvalarray[,,i,idx]
          for(i in 1:dim(fittedvalarray)[3]) WSgtemp[,,i] = WSdtemp[,,i,1]
          
        }
        }

        ENStemp[,,,g] = WSgtemp 
    }

	GCMcheck = c()
	for(c in 1:dim(ENStemp)[4]) GCMcheck[c] = ifelse(all(is.na(ENStemp[ceiling(dim(fittedvalarray)[1]/2),ceiling(dim(fittedvalarray)[2]/2),,c])==TRUE)==TRUE,0,1)

	GCMidx=which(GCMcheck!=0)
	ENStempuse = ENStemp[,,,GCMidx]

	for(n in 1:length(GCMidx)){
	Wgnorm = (weightsGCM[,,GCMidx[n]]/apply(weightsGCM[,,GCMidx],c(1,2),sum))
	for(i in 1:dim(fittedvalarray)[3]) ENStempuse[,,i,n] = ENStempuse[,,i,n]*Wgnorm
	}

    ENSout[,,,s] = apply(ENStempuse,c(1,2,3),sum,na.rm=TRUE) 
    
  }

  S = apply(ENSout,c(1,2,3),var,na.rm=TRUE)
  S = ifelse(S==0,NA,S)
  S
}

######################
# Scenario Uncertainty Calculator - generalized version

scenuncold = function(fittedvalarray,scens,DSs,GCMs,weightsGCM=GCMweights,weightsDS=DSweights,files=files){
  
  scens = scens
  GCMs = GCMs
  
  ENSout = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(scens)))

  for(s in 1:length(scens)){
  
    ENStemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(GCMs)))

    for(g in 1:length(GCMs)){
    
        scenidx = grep(scens[s],files)
        idx1 = grep(GCMs[g],files)
        
        idx = idx1[which(idx1 %in% scenidx)]
      
        WSdtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(idx)))
        WSgtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))

        if(length(idx)!=0){
        if(length(idx)>1 & length(idx)!=0){
        
        for(i in 1:length(idx)){
           for(t in 1:dim(fittedvalarray)[3]) WSdtemp[,,t,i] = fittedvalarray[,,t,idx[i]]*weightsDS[,,which(DSs==strsplit(files[idx[i]],"_")[[1]][3])]         
        }
        
        WSgtemp = apply(WSdtemp,c(1,2,3),sum)
        for(i in 1:dim(fittedvalarray)[3]) WSgtemp[,,i] = WSgtemp[,,i]*weightsGCM[,,g]
        
        } else {
          
          for(i in 1:dim(fittedvalarray)[3]) WSdtemp[,,i,1] = fittedvalarray[,,i,idx]*weightsDS[,,which(DSs==strsplit(files[idx],"_")[[1]][3])]
          for(i in 1:dim(fittedvalarray)[3]) WSgtemp[,,i] = WSdtemp[,,i,1]*weightsGCM[,,g]
          
        }
        }
        
        ENStemp[,,,g] = WSgtemp 
    }

    ENSout[,,,s] = apply(ENStemp,c(1,2,3),sum,na.rm=TRUE) 
    
  }

  S = apply(ENSout,c(1,2,3),var,na.rm=TRUE)
  S = ifelse(S==0,NA,S)
  S
}

#################
# Mean Change Calculator - Generalized,time series

meanchangets=function(fittable,scens,GCMs,DSs,weightsGCM,weightsDS,files=files){
  
  fittable = fittedvals
  files= memnames
  weightsGCM=GCMweights
  weightsDS = DSweights
  
  ENStemp = matrix(NA,nrow=nrow(fittable),length(GCMs))
  
  for(g in 1:length(GCMs)){
    
    idx = grep(GCMs[g],files)
    IDX=idx+1
    
    #if(g==2) idx=idx[c(1,2,5,6,7)]
    #if(g==17) idx=idx[c(3,4,8,9,10)]
    WSgtemp = c()
    
    if(length(idx)>0){
      
      WSdtemp = matrix(NA,nrow=nrow(fittable),length(idx))
      
      if(length(idx)>1){
        
        filesplit = do.call(rbind,strsplit(files[idx],"_"))
        
        for(i in 1:length(idx)){
          #Wdnorm = (weightsDS[,,which(DSs %in% filesplit[i,3])]/apply(weightsDS[,,which(DSs %in% filesplit[,3])],c(1,2),sum))
          Wdnorm = (weightsDS[which(DSs %in% paste(filesplit[i,2],filesplit[i,3],sep="_"))]/sum(weightsDS[which(DSs %in% paste(filesplit[,2],filesplit[,3],sep="_"))]))
          for(t in 1:nrow(fittable)) WSdtemp[t,i] = fittable[t,IDX[i]]*Wdnorm
        }
        WSgtemp = apply(WSdtemp,1,sum,na.rm=TRUE)
        
      } else {
        
        for(t in 1:nrow(fittable)) WSdtemp[t,1] = fittable[t,IDX[i]]
        
        WSgtemp = WSdtemp[,1]
        
      }
      
      #for(t in 1:dim(fittable)[3]) WSgtemp[,,t] = WSgtemp[,,t]*weightsGCM[,,g]
      
      ENStemp[,g] = WSgtemp
    }
  }
  
  GCMcheck = c()
  for(c in 1:ncol(ENStemp)) GCMcheck[c] = ifelse(all(is.na(ENStemp[,c])==TRUE)==TRUE,0,1)
  
  GCMidx=which(GCMcheck!=0)
  ENStempuse = ENStemp[,GCMidx]
  
  for(n in 1:length(GCMidx)){
    Wgnorm = (weightsGCM[GCMidx[n]]/sum(weightsGCM[GCMidx]))
    for(i in 1:nrow(fittable)) ENStempuse[i,n] = ENStempuse[i,n]*Wgnorm
  }
  
  ENStempuse = apply(ENStempuse,1,sum,na.rm=TRUE)
  G = ENStempuse/length(scens)
  for(t in 1:nrow(fittable)) G[t] = ifelse(is.na(weightsGCM[1])==TRUE,NA,G[t])
  G
  
}


#################
# Mean Change Calculator - Generalized

meanchange=function(fittable,scens,GCMs,DSs,weightsGCM,weightsDS,files=files){

ENStemp = array(NA,dim=c(dim(fittable)[1],dim(fittable)[2],dim(fittable)[3],length(GCMs)))

    for(g in 1:length(GCMs)){

      idx = grep(GCMs[g],files)

      #if(g==2) idx=idx[c(1,2,5,6,7)]
      #if(g==17) idx=idx[c(3,4,8,9,10)]
      WSgtemp = array(NA,dim=c(dim(fittable)[1],dim(fittable)[2],dim(fittable)[3]))

      if(length(idx)>0){

WSdtemp = array(NA,dim=c(dim(fittable)[1],dim(fittable)[2],dim(fittable)[3],length(idx)))

      if(length(idx)>1){

	filesplit = do.call(rbind,strsplit(files[idx],"_"))

      for(i in 1:length(idx)){
        #Wdnorm = (weightsDS[,,which(DSs %in% filesplit[i,3])]/apply(weightsDS[,,which(DSs %in% filesplit[,3])],c(1,2),sum))
        Wdnorm = (weightsDS[,,which(DSs %in% paste(filesplit[i,3],filesplit[i,4],sep="_"))]/apply(weightsDS[,,which(DSs %in% paste(filesplit[,3],filesplit[,4],sep="_"))],c(1,2),sum))
	for(t in 1:dim(fittable)[3]) WSdtemp[,,t,i] = fittable[,,t,i]*Wdnorm
      }
      WSgtemp = apply(WSdtemp,c(1,2,3),sum,na.rm=TRUE)

      } else {

       for(t in 1:dim(fittable)[3]) WSdtemp[,,t,1] = fittable[,,t,idx]

        WSgtemp = WSdtemp[,,,1]

      }

      #for(t in 1:dim(fittable)[3]) WSgtemp[,,t] = WSgtemp[,,t]*weightsGCM[,,g]

      ENStemp[,,,g] = WSgtemp
    }
    }

	GCMcheck = c()
	for(c in 1:dim(ENStemp)[4]) GCMcheck[c] = ifelse(all(is.na(ENStemp[ceiling(dim(fittable)[1]/2),ceiling(dim(fittable)[2]/2),,c])==TRUE)==TRUE,0,1)

	GCMidx=which(GCMcheck!=0)
	ENStempuse = ENStemp[,,,GCMidx]

	for(n in 1:length(GCMidx)){
	Wgnorm = (weightsGCM[,,GCMidx[n]]/apply(weightsGCM[,,GCMidx],c(1,2),sum))
	for(i in 1:dim(fittable)[3]) ENStempuse[,,i,n] = ENStempuse[,,i,n]*Wgnorm
	}

    ENStempuse = apply(ENStempuse,c(1,2,3),sum,na.rm=TRUE)
    G = ENStempuse/length(scens)
    for(t in 1:dim(fittable)[3]) G[,,t] = ifelse(is.na(weightsGCM[,,1])==TRUE,NA,G[,,t])
    G

  }



######################
# Mean Change Calculator - specialized

meanchangeold=function(fittable,scens,GCMs,DSs,weightsGCM,weightsDS,files=files){

ENStemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(GCMs)))

    for(g in 1:length(GCMs)){

      idx = grep(GCMs[g],files)

      #if(g==2) idx=idx[c(1,2,5,6,7)]
      #if(g==17) idx=idx[c(3,4,8,9,10)]
      WSgtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3]))

      if(length(idx)>0){

WSdtemp = array(NA,dim=c(dim(fittedvalarray)[1],dim(fittedvalarray)[2],dim(fittedvalarray)[3],length(idx)))

      if(length(idx)>1){
      for(i in 1:length(idx)){
        for(t in 1:dim(fittedvalarray)[3]) WSdtemp[,,t,i] = fittedvalarray[,,t,i]*weightsDS[,,which(DSs==strsplit(files[idx[i]],"_")[[1]][3])]
      }
      WSgtemp = apply(WSdtemp,c(1,2,3),sum,na.rm=TRUE)

      } else {

       for(t in 1:dim(fittedvalarray)[3]) WSdtemp[,,t,1] = fittedvalarray[,,t,idx]*weightsDS[which(DSs==strsplit(files[idx],"_")[[1]][3])]

        WSgtemp = WSdtemp[,,,1]

      }

      for(t in 1:dim(fittedvalarray)[3]) WSgtemp[,,t] = WSgtemp[,,t]*weightsGCM[,,g]

      ENStemp[,,,g] = WSgtemp
    }
    }

    ENStemp = apply(ENStemp,c(1,2,3),sum,na.rm=TRUE)
    G = ENStemp/length(scens)
    for(t in 1:dim(fittedvalarray)[3]) G[,,t] = ifelse(is.na(weightsGCM[,,1])==TRUE,NA,G[,,t])
    G

  }


