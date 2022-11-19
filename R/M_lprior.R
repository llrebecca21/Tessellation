M_lprior=function(n,p,M,Mmax)
{
  a=-log(p) - log(Mmax)
  
  b=1*(n/M)
  if(M==1)
  {
    b=b*1
    
  }else{
    
    for(i in 1:(M-1))
    {
      b=b*((n-i)/(M-i))
    }
    
  }
  
  a=a-log(b)
  
  return(a)
  
}