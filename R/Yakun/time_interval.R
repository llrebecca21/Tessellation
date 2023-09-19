time_interval=function(x,S,prt,NumObs)
{
  # inputs: x     : 20_000 x 5 matrix
  #         S     : 8 x 1 col vec
  #         prt   : 20_000 x 1 matrix
  #         NumObs: "numeric" 20
  
  # initialize variables
  M <- length(S)
  # Create matrix of 0 to store interval, tmin, tmax
  interval=matrix(0, nrow = M, ncol = NumObs) # 8 x 20
  tmin=matrix(0, nrow = M, ncol = NumObs)     # 8 x 20
  tmax=matrix(0, nrow = M, ncol = NumObs)     # 8 x 20
  
  # For each center element in S
  for(i in 1:length(S))
  {
    # subset x where the ith group belongs to the ith center
    x_s = x[which(prt == i),]
    
    # For 1 through 20:
    for(j in 1:NumObs)
    {
      # idx keeps track of 'true' groups
      idx=which(x_s[,"index"]==j) 
      # do not include jth obs in ith tessellation
      if(identical(idx, integer(0)))
      {
        tmin[i,j]= -1
        tmax[i,j]= -1
        interval[i,j]= -1
        
      }else{
        minmax <- range(x_s[idx,"time"])
        # Find minimum time of jth group within ith tessellation
        # tmin[i,j]=min(x_s[idx,"time"])
        tmin[i,j] = minmax[1]
        # Find maximum time of jth group within ith tessellation
        # tmax[i,j]=max(x_s[idx,"time"])
        tmax[i,j] = minmax[2]
        # Find inclusive range of the max and min times
        # interval[i,j]=max(x_s[idx,"time"])-min(x_s[idx,"time"])+1
        interval[i,j] = diff(minmax) + 1
        
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
