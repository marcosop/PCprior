######################################################################
######## Auxiliar functions ##########################################
######################################################################
liklogitpcp=function(theta,matx,y){
  q=length(theta)
  p=q-1
  alpha=theta[p+1]
  plogit=(pplogis(as.vector(matx%*%theta[1:p])))^theta[p+1]
  plogitqlogit=log(dbinom(y,1,plogit))
  liklogit=sum(plogitqlogit)
  return(liklogit)
}

liklogit=function(theta,matx,y){
  plogit=(plogis(as.vector(matx%*%theta)))
  plogitqlogit=log(dbinom(y,1,plogit))
  liklogit=sum(plogitqlogit)
  return(liklogit)
}


termilogit=function(theta,matx,y){
  q=length(theta)
  p=q-1
  alpha=theta[p+1]
  plogit=(pplogis(as.vector(matx%*%theta[1:p])))^theta[p+1]
  plogitqlogit=log(dbinom(y,1,plogit))
  return(plogitqlogit)
}

termilogit1=function(theta,matx,y){
  plogit=plogis(as.vector(matx%*%theta))
  plogitqlogit=log(dbinom(y,1,plogit))
  return(plogitqlogit)
}


adjplogis=function(theta,matx){
  q=length(theta)
  p=q-1
  plogit=pplogis(as.vector(matx%*%theta[1:p]),lambda=theta[p+1])
  return(plogit)
}

adjlogis=function(theta,matx){
  plogit=plogis(as.vector(matx%*%theta))
  return(plogit)
}

quantres=function(resprob){
  yi=resprob[1]
  padi=resprob[2]
  a=pbeta(1-padi,2-yi,yi)
  b=pbeta(1-padi,1-yi,yi+1)
  ui=runif(1,a,b)
  res=qnorm(ui)
  return(res)
}


claslogis=function(pad,p,y){
  clasest=as.numeric(pad>p)
  N=length(y[y==0])
  P=length(y[y==1])
  tp=sum(y==1 & clasest==1)
  fn=sum(y==1 & clasest==0)
  fp=sum(y==0 & clasest==1)
  tn=sum(y==0 & clasest==0)
  tpr=tp/P
  tnr=tn/N
  csi = tp/(tp+fp+fn)
  acc=(tp+tn)/(P+N)
  ssi = tp/(tp + 2*(fp+fn))
  faith = (tp + 0.5*fn)/(P+N)
  pdif = (4*fp*fn) /((P+N)^2)
  table1=c(tp,fp)
  table2=c(fn,tn)
  table=rbind(table1,table2)
  return(list(table=table,tpr=tpr,tnr=tnr,acc=acc,csi=csi,ssi=ssi,faith=faith, pdif=pdif))
}

######################################################################
######## Auxiliar functions ##########################################
######################################################################

likprobitpcp=function(theta,matx,y){
  q=length(theta)
  p=q-1
  alpha=theta[p+1]
  pprobit=(ppnorm(as.vector(matx%*%theta[1:p])))^theta[p+1]
  pprobitqprobit=log(dbinom(y,1,pprobit))
  likprobit=sum(pprobitqprobit)
  return(likprobit)
}

likprobit=function(theta,matx,y){
  pprobit=(pnorm(as.vector(matx%*%theta)))
  pprobitqprobit=log(dbinom(y,1,pprobit))
  likprobit=sum(pprobitqprobit)
  return(likprobit)
}


termiprobit=function(theta,matx,y){
  q=length(theta)
  p=q-1
  alpha=theta[p+1]
  pprobit=(ppnorm(as.vector(matx%*%theta[1:p])))^theta[p+1]
  pprobitqprobit=log(dbinom(y,1,pprobit))
  return(pprobitqprobit)
}

termiprobit1=function(theta,matx,y){
  pprobit=pnorm(as.vector(matx%*%theta))
  pprobitqprobit=log(dbinom(y,1,pprobit))
  return(pprobitqprobit)
}



adjpnorm=function(theta,matx){
  q=length(theta)
  p=q-1
  pprobit=ppnorm(as.vector(matx%*%theta[1:p]),lambda=theta[p+1])
  return(pprobit)
}

adjprobit=function(theta,matx){
  pprobit=pnorm(as.vector(matx%*%theta))
  return(pprobit)
}


claslogis=function(pad,p,y){
  clasest=as.numeric(pad>p)
  N=length(y[y==0])
  P=length(y[y==1])
  tp=sum(y==1 & clasest==1)
  fn=sum(y==1 & clasest==0)
  fp=sum(y==0 & clasest==1)
  tn=sum(y==0 & clasest==0)
  tpr=tp/P
  tnr=tn/N
  csi = tp/(tp+fp+fn)
  acc=(tp+tn)/(P+N)
  ssi = tp/(tp + 2*(fp+fn))
  faith = (tp + 0.5*fn)/(P+N)
  pdif = (4*fp*fn) /((P+N)^2)
  table1=c(tp,fp)
  table2=c(fn,tn)
  table=rbind(table1,table2)
  return(list(table=table,tpr=tpr,tnr=tnr,acc=acc,csi=csi,ssi=ssi,faith=faith, pdif=pdif))
}

