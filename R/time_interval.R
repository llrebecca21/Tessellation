time_interval=function(x,S,prt,NumObs)
{
  
  interval=matrix(0,length(S),NumObs)
  tmin=matrix(0,length(S),NumObs)
  tmax=matrix(0,length(S),NumObs)
  
  for(i in 1:length(S))
  {
    x_s=x[which(prt==i),]
    for(j in 1:NumObs)
    {
      idx=which(x_s[,"index"]==j)   # do not include jth obs in ith tessellation
      if(identical(idx, integer(0)))
      {
        tmin[i,j]=-1
        tmax[i,j]=-1
        interval[i,j]=-1
        
      }else{
        
        tmin[i,j]=min(x_s[idx,"time"])
        tmax[i,j]=max(x_s[idx,"time"])
        interval[i,j]=max(x_s[idx,"time"])-min(x_s[idx,"time"])+1
        
      }
      
      
    }
    
  }
  
  list(interval=interval,tmin=tmin,tmax=tmax)
  
}


# ### Decide if minimum time points > Tmin
# time_interval_can=function(x,S,Tmin,prt,NumObs)
# {
#   ## decide whether the time interval > Tmin
#   
#   can=1 
#   interval=matrix(0,length(S),NumObs)
#   for(i in 1:length(S))
#   {
#     
#     x_s=x[which(prt==i),]
#     for(j in 1:NumObs)
#     {
#       
#       idx=which(x_s[,"index"]==j)
#       if(identical(idx, integer(0)))
#       {
#         interval[i,j]=-1
#         
#       }else{
#         
#         interval[i,j]=max(x_s[idx,"time"])-min(x_s[idx,"time"])+1
#       }        
#       
#       if(interval[i,j]<Tmin & interval[i,j]!=-1)
#       {
#         can=0
#         break
#         
#       }
#       
#     }
#     
#   }
#   
#   return(can)
#   
# }
# 

#time_interval=time_interval(x,S,prt,NumObs)
