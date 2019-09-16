###########################
#
# SKC ANALYSIS FUNCTIONS
#
###########################

##################
#
# EOF Analysis - first test with the base function svd 

SVDanalysis = function(combodat,trunc=10,return.all=FALSE){
  
  svdvals = svd(combodat)
  
  U = svdvals$u # this is an m by m matrix, at this untruncated to the t modes
  V = svdvals$v # this is the n by m matrix - eventually n by t also untruncated
  
  lambda = svdvals$d # diagonal eigenvalues lambda
  
  ####
  # what's the order to eigenvalues in decreasing order?
  # if it's not already in decreasing order
  
  lorder = order(lambda,decreasing=TRUE)
  lambda = lambda[lorder]
  U = U[,lorder]
  V = V[,lorder]
  
  #####
  # truncate U, V, and lambda to top 10 modes
  # not really any guidance in SKC on how things were truncated so I'm unsure
  
  Ut = U[,1:trunc] # this will now be m rows by t cols as in SKC
  Vt = V[,1:trunc] # this will now be n rows by t cols as in SKC also
  lt = lambda[1:trunc] # keeping highest 10 modes
  
  if(return.all==FALSE){
    return(Ut)
  } else {
    return(list(Ut,Vt,lt))
  }
  
}

#########################
# Euclidean intermodel distances

####
# intermodel distance calc function

distfunc = function(Ut,i,j){
  diff1 = Ut[i,]-Ut[j,]
  diff2 = diff1^2
  innersum = sum(diff2)
  distval = sqrt(innersum)
  return(distval)
}

##########
# Distance matrix calculation loop

distmatcalc = function(Ut){
  
  distmat = matrix(NA,nrow=nrow(Ut),ncol=nrow(Ut))
  
  for(I in 1:nrow(Ut)){
    for(J in 1:nrow(Ut)){
      distmat[I,J] = distfunc(Ut,I,J) 
    }
  }
  
  return(distmat)
  
}

######################
# MDS calculator

MDScalc = function(delta,initlist=inits){
  
  results = list()
  delta = as.dist(delta)
  for(i in 1:length(initlist)){
    results[[i]]= mds(delta, ndim=2,init=initlist[[i]])$conf
  }
  
  return(results)
  
}

######################
# MDS calculator 2

MDScalc2 = function(delta,listlength=20){
  results=list()
  for(i in 1:listlength){
    results[[i]] = mds(delta, ndim=2,init=NULL,itmax=1000)$conf
  }
  return(results)  
}


########################
# Intermodel distance calculator

intermoddist = function(mdsvals){
  
  distmat = matrix(NA,nrow=nrow(mdsvals),ncol=nrow(mdsvals))
  
  for(i in 1:nrow(mdsvals)){
    for(j in 1:nrow(mdsvals)){
      
      diff1 = (mdsvals[i,1]-mdsvals[j,1])^2
      diff2 = (mdsvals[i,2]-mdsvals[j,2])^2
      
      distmat[i,j] = sqrt((diff1+diff2))
      
    }
  }
  
  return(distmat)  
}

#######################
# intermod distance - initializations to array

intermodarray = function(mdsvals,numsims=20){
  
  distvals = array(data=NA,dim=c(nrow(mdsvals[[1]]),nrow(mdsvals[[1]]),numsims))
  
  for(i in 1:numsims){
    distvals[,,i] = intermoddist(mdsvals[[i]])
  }
  avg = apply(distvals,c(1,2),mean,na.rm=TRUE)
  
  return(list(distvals,avg))
  
}

