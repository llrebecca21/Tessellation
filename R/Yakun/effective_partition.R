effective_partition=function(x,w,M,NumObs,Ntime,Tmin)
{
  stp=0
  ii=0
  while(stp==0)
  {
    ii=ii+1
    cat("ii is",ii,"\n")
    #S=sample(1:(NumObs*Ntime),M)
    S=sample(seq(Tmin,Ntime*NumObs-Tmin,by=2*Tmin),M)
    #S=sample(cutpoint,M)
    
    #prt=distance_partition(x,S,w)
    prt=distance_partitionC(as.matrix(x[,-c(1,2)]),S,w)
    can=time_interval_can(x,S,Tmin,prt,NumObs)
    if(can==1 & length(unique(prt))==M)
    {
      stp=1
    }
    
  }
  
  time_interval=time_interval(x,S,prt,NumObs)
  
  list(S=S,prt=prt,can=can,time_interval=time_interval)
  
}
