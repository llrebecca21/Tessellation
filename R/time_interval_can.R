### Decide if minimum time points > Tmin
time_interval_can=function(x,S,Tmin,prt,NumObs)
{
  ## decide whether the time interval > Tmin
  
  can=1 
  interval=matrix(0,length(S),NumObs)
  for(i in 1:length(S))
  {
    
    x_s=x[which(prt==i),]
    for(j in 1:NumObs)
    {
      
      idx=which(x_s[,"index"]==j)
      if(identical(idx, integer(0)))
      {
        interval[i,j]=-1
        
      }else{
        
        interval[i,j]=max(x_s[idx,"time"])-min(x_s[idx,"time"])+1
      }        
      
      if(interval[i,j]<Tmin & interval[i,j]!=-1)
      {
        can=0
        break
        
      }
      
    }
    
  }
  
  return(can)
  
}
