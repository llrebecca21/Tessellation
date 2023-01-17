rm(list = ls())
# setwd("C:/Users/yakun_wang/Desktop/tessellation/code_project2")
source('M_lprior.R')
source('tessellation_likelihood.R')
source('time_interval.R')
source('MHstep_BFGS_11-19.R')
source('effective_partition.R')
Rcpp::sourceCpp('distance_partitionC_Lee.cpp.')
Rcpp::sourceCpp('distance_partition.cpp')


# Load required libraries
library(pracma)
library(trust)
library(mvtnorm)
library(MASS)
library(NPflow)
library(plotly)
library(invgamma)
library(gtools)
library(forecast)



range01 <- function(x){(x-min(x))/(max(x)-min(x))}


###################################
# Abrupt-abrupt simulate data
#################################
set.seed(1080)
Ntime=1000
NumX=2
NumXT=NumX+1
NumObs=20

## create covariates matrix of x
xx=matrix(0,NumObs,NumX)


## True S is 2,8,13,18
xx[1,]=c(0.2,0.4)
xx[2,]=c(0.3,0.3)
xx[3,]=c(0.4,0.2)
xx[4,]=c(0.2,0.2)
xx[5,]=c(0.4,0.4)
xx[6,]=c(0.2,0.6)  
xx[7,]=c(0.2,0.8)
xx[8,]=c(0.3,0.7)  
xx[9,]=c(0.4,0.6)  
xx[10,]=c(0.4,0.8)  
xx[11,]=c(0.6,0.2)  
xx[12,]=c(0.6,0.4)  
xx[13,]=c(0.7,0.3)
xx[14,]=c(0.8,0.2)
xx[15,]=c(0.8,0.4)
xx[16,]=c(0.6,0.6)
xx[17,]=c(0.6,0.8)
xx[18,]=c(0.7,0.7)
xx[19,]=c(0.8,0.6)
xx[20,]=c(0.8,0.8)

# Redo code to create column names #
x = matrix(data = 0, nrow = Ntime*NumObs, ncol = NumXT+2)
# Create column names for x1 and x2
xname = paste("x", 1:(NumX), sep = "")
# Create name vector
names = append(c("index", "time", "scaled_time"),xname)
# change matrix x into a dataframe
x = data.frame(x)
# rename the columns of x
colnames(x) = names


#x=matrix(0,Ntime*NumObs,NumXT+2)
#names=c("index","time","scaled_time")
#for(i in 1:NumX)
#{
#  xname=paste("x",i,sep = "")
#  names=append(names,xname)
#}
#x=data.frame(x)
#colnames(x)=names


# Fill in the columns of x #
# replace index with 1:20 with each repeated 1_000 times
x$index=rep(1:NumObs,each=Ntime)
# replace time with 1:1000 repeated 20 times
x$time=rep(1:Ntime,NumObs)
# replace scaled_time with proportion of length completed; repeats each 1_000 index  
x$scaled_time=range01(x$time)


# Repeat each row of xx 1000 times 
x[,4:(NumXT+2)]=apply(X = xx, MARGIN = 2, FUN = function(i) rep(i,each=Ntime))


## create time series x_t, abrupt change at x2, and at half time points
#ts.sim=matrix(rep(0,NumObs*Ntime),NumObs,Ntime)
#ts.sim1=matrix(rep(0,NumObs*Ntime/2),NumObs,Ntime/2)
#ts.sim2=matrix(rep(0,NumObs*Ntime/2),NumObs,Ntime/2)

## Redo ts.sim, ts.sim1, ts.sim2 ##
# ts.sim = matrix(data = 0, nrow = NumObs, ncol = Ntime)
ts.sim1 = matrix(data = 0, nrow = Ntime/2, ncol = NumObs)
ts.sim2 = matrix(data = 0, nrow = Ntime/2, ncol = NumObs)


for(i in 1:NumObs)
{
  if(xx[i,1]<0.5)
  {
    
    if(xx[i,2]<0.5)
    {
      ts.sim1[i,]=arima.sim(list(order=c(1,0,0),ar=0.3),n=Ntime/2)
      ts.sim2[i,]=arima.sim(list(order=c(1,0,0),ar=-0.5),n=Ntime/2)
      
    }else{
      
      ts.sim1[i,]=arima.sim(list(order=c(1,0,0),ar=-0.9),n=Ntime/2)
      ts.sim2[i,]=arima.sim(list(order=c(1,0,0),ar=0.7),n=Ntime/2)
      
    }
    
    
  }else{
    
    if(xx[i,2]<0.5)
    {
      
      ts.sim1[i,]=arima.sim(list(order=c(1,0,0),ar=-0.3),n=Ntime/2)
      ts.sim2[i,]=arima.sim(list(order=c(1,0,0),ar=0.5),n=Ntime/2)
      
    }else{
      
      ts.sim1[i,]=arima.sim(list(order=c(1,0,0),ar=0.9),n=Ntime/2)
      ts.sim2[i,]=arima.sim(list(order=c(1,0,0),ar=-0.7),n=Ntime/2)
      
    }
    
  }
}


x_t = rbind(ts.sim1,ts.sim2)  




# standardized
#for (i in 1:NumObs)
#{
#  xmat=cbind(matrix(1,dim(x_t)[1],1), matrix(seq(1,dim(x_t)[1],1),dim(x_t)[1],1))
#  linfit=solve(t(xmat)%*%xmat)%*%t(xmat)%*%x_t[,i]
#  x_t[,i]=x_t[,i]-xmat%*%linfit
#}
#ts.plot(x_t[,20])


#plot(x[,c(4,5)])

# Redo Standardized #
xmat = cbind(1,1:Ntime)
linfit = solve(crossprod(xmat), crossprod(xmat,x_t))
x_t = x_t - xmat %*% linfit

ts.plot(x_t[,20])

plot(x[,c(4,5)])





###################################
# slowly-abrupt simulate data
#################################

set.seed(1080)
Ntime=1000
NumX=2
NumXT=NumX+1
NumObs=20

## create covariates matrix of x
xx=matrix(0,NumObs,NumX)


## True S is 2,8,13,18
xx[1,]=c(0.2,0.4)
xx[2,]=c(0.3,0.3)
xx[3,]=c(0.4,0.2)
xx[4,]=c(0.2,0.2)
xx[5,]=c(0.4,0.4)
xx[6,]=c(0.2,0.6)  
xx[7,]=c(0.2,0.8)
xx[8,]=c(0.3,0.7)  
xx[9,]=c(0.4,0.6)  
xx[10,]=c(0.4,0.8)  
xx[11,]=c(0.6,0.2)  
xx[12,]=c(0.6,0.4)  
xx[13,]=c(0.7,0.3)
xx[14,]=c(0.8,0.2)
xx[15,]=c(0.8,0.4)
xx[16,]=c(0.6,0.6)
xx[17,]=c(0.6,0.8)
xx[18,]=c(0.7,0.7)
xx[19,]=c(0.8,0.6)
xx[20,]=c(0.8,0.8)







#x = matrix(0,Ntime*NumObs,NumXT+2)
#names=c("index","time","scaled_time")
#for(i in 1:NumX)
#{
#  xname=paste("x",i,sep = "")
#  names=append(names,xname)
#}
#x=data.frame(x)
#colnames(x)=names

# Fill in the columns of x #
# replace index with 1:20 with each repeated 1_000 times
x$index=rep(1:NumObs,each=Ntime)
# replace time with 1:1000 repeated 20 times
x$time=rep(1:Ntime,NumObs)
# replace scaled_time with proportion of length completed; repeats each 1_000 index  
x$scaled_time=range01(x$time)


# Repeat each row of xx 1000 times 
x[,4:(NumXT+2)]=apply(X = xx, MARGIN = 2, FUN = function(i) rep(i,each=Ntime))



## create slowly change of time series x_t, abrupt change at X1, x2
sd_true = 1
e = matrix(data = rnorm(n = NumObs * Ntime, mean = 0, sd = sd_true), nrow = NumObs, ncol = Ntime)
ts.sim = matrix(data = 0, nrow = NumObs, ncol = Ntime)


for(i in 1:NumObs)
{
  if(xx[i,1] < 0.5)
  {
    
    if(xx[i,2] < 0.5)
    {
      # slowly varying -0.3 to 0.3
      phi_true = -0.3 + ((1:Ntime) / Ntime) * 0.6
      for (j in 1:Ntime)
      {
        if (j==1)
        {
          ts.sim[i,j] = e[i,j]
          
        }else
        {
          ts.sim[i,j] = phi_true[j] * ts.sim[i,j-1] + e[i,j]
        }
        
      }
      
      
    }else{
      
      # slowly varying from 0.9 to -0.9
      phi_true=0.9-((1:Ntime)/Ntime)*1.8
      for (j in 1:Ntime)
      {
        if (j==1)
        {
          ts.sim[i,j]=e[i,j]
          
        }else
        {
          ts.sim[i,j]=phi_true[j]*ts.sim[i,j-1]+e[i,j]
        }
      } 
      
    }
    
    
  }else{
    
    if(xx[i,2]<0.5)
    {
      
      # slowly varying from 0.7 to -0.7
      phi_true=0.7-((1:Ntime)/Ntime)*1.4
      for (j in 1:Ntime)
      {
        if (j==1)
        {
          ts.sim[i,j]=e[i,j]
          
        }else
        {
          ts.sim[i,j]=phi_true[j]*ts.sim[i,j-1]+e[i,j]
        }
      } 
      
    }else{
      
      # slowly varying from -0.5 to 0.5
      phi_true=-0.5+((1:Ntime)/Ntime)
      for (j in 1:Ntime)
      {
        if (j==1)
        {
          ts.sim[i,j]=e[i,j]
          
        }else
        {
          ts.sim[i,j]=phi_true[j]*ts.sim[i,j-1]+e[i,j]
        }
      } 
      
    }
    
  }
}


x_t=t(ts.sim)



# standardized
for (i in 1:NumObs)
{
  xmat=cbind(matrix(1,dim(x_t)[1],1), matrix(seq(1,dim(x_t)[1],1),dim(x_t)[1],1))
  linfit=solve(t(xmat)%*%xmat)%*%t(xmat)%*%x_t[,i]
  x_t[,i]=x_t[,i]-xmat%*%linfit
}
ts.plot(x_t[,1])


plot(x[,c(4,5)])






###############################################
#  Hyperparameters for tessellation partition
##################################################
set.seed(2333)


## Initialization of tessellation
M=8 # Number of centers
S=c(1250,1750,7250,7750,12250,12750,17250,17750)
w=c(0.04190068, 0.63559150, 0.32250782)


# s_num=seq(5,995,by=(1000-1)/(10-1))
# s2=s_num+1000
# s8=s_num+7000
# s13=s_num+12000
# s18=s_num+17000
# 
# 
# S=c(s2,s8,s13,s18)
# w=c(0.07359252, 0.46466762, 0.46173986)



prt=distance_partitionC(as.matrix(x[,-c(1,2)]),S,w)
interval_curr=time_interval(x,S,prt,NumObs)


# Max number of tessellation
Mmax=10

# Tmin
Tmin=50



## MCMC parameters
nbasis=7
nbeta=nbasis+1
sigmasqalpha=100



# Initialization of half-t distribution
tau=rep(10,M)
g=rep(1,M)


nus=10 #2
Gs=2 #10
param_random=0.2
nb_alpha=nbeta



## calculate log power spectrum
nfreq=50  # fixed the number of Fourier frequencies
freq <- (0:nfreq)/(2 * nfreq)
nu_mat <- lin_basis_func(freq, nbeta)
fhat_prop=matrix(0,(nfreq+1),M)



loglike=rep(0,M)
log_beta=rep(0,M)
beta_lprior=rep(0,M)
beta_prop=matrix(0,M,nbeta)


# update beta and loglike
for(i in 1:M)
{
  prop_obs=x[which(prt==i),]
  prop_obs_index=unique(prop_obs$index)
  
  # create log periodograms for replications in ith time seg and jth cov seg
  nseg_time_temp=interval_curr$interval[i,]
  tmin=interval_curr$tmin[i,]
  tmax=interval_curr$tmax[i,]
  
  y <- list()
  yy <- list()
  uu=0
  
  for(l in prop_obs_index)
  {
    uu=uu+1
    
    nfreq <- floor(nseg_time_temp[l] / 2)
    
    y[[uu]] <- log(abs(fft(x_t[tmin[l]:tmax[l], l])) ^ 2 / nseg_time_temp[l])
    yy[[uu]] <- y[[uu]][1:(nfreq + 1)]
    
  }
  
  
  # initial beta
  param=rep(0,nbeta)
  
  
  # BFGS
  postbeta_BFGS <- optim(param,fn, gr, method="BFGS", prop_obs_index, interval_curr$interval[i,], yy, tau[i], nbeta, nbasis, sigmasqalpha)
  Q=he(postbeta_BFGS$par,prop_obs_index, interval_curr$interval[i,], yy, tau[i], nbeta, nbasis, sigmasqalpha)
  Lt=chol(Q)
  
  
  ### sampling
  beta_prop[i,]=Chol_sampling(Lt,nbeta,postbeta_BFGS$par)
  
  ## Density
  log_beta[i]=Chol_density(beta_prop[i,],Lt,nbeta,postbeta_BFGS$par)
  
  
  
  ## prior density
  beta_lprior[i]=dmvnorm(beta_prop[i,],matrix(0,nbeta,1),diag(c(sigmasqalpha, tau[i]*matrix(1,nbasis,1))),log = T) # Prior Density of beta
  
  
  # estimated log of power spectrum and Whittlelikelihood
  fhat_prop[,i]=nu_mat%*%matrix(beta_prop[i,],nbeta,1)
  log_prop_spec_dens=whittle_like(y,prop_obs_index,interval_curr$interval[i,],beta_prop[i,],nbasis) # Loglikelihood  at proposed values
  loglike[i]=log_prop_spec_dens
  
  
}



## New prior of tessellation
p=NumXT
n=dim(x)[1]
log_M_prior=M_lprior(n,p,M,Mmax)

### current tessellation structure
TT=list("S"=S,"M"=M,"w"=w,"loglike"=loglike,"beta"=beta_prop,"tau"=tau,"g"=g,"log_beta"=log_beta,"M_lprior"=log_M_prior,"beta_lprior"=beta_lprior,"fhat_prop"=fhat_prop,"prt"=prt)





######################
#        MCMC 
######################
nwarmup=5000
nloop=10000
pre_beta=FALSE


result=matrix(list(),nloop,1)


for(p in 1:nloop)
{
  cat("p is ",p,"\n")
  result[[p,1]]=Whithin(x,x_t,TT,NumObs,nbasis,sigmasqalpha,pre_beta)
  
  
  ## update tau & g
  TT=result[[p,1]]$TT
  for(ss in 1:TT$M)
  {
    
    g_a=(nus + 1)/2
    g_b=(nus)/TT$tau[ss] + 1/(Gs^2)
    TT$g[ss] = 1/rgamma(1,g_a,scale=1/g_b)
    
    tau_a=(nb_alpha - 1 + nus)/2
    tau_b=sum(TT$beta[ss,2:nbeta]^2)/2+nus/TT$g[ss]
    TT$tau[ss] = 1/rgamma(1,tau_a,rate=tau_b)
    
  }
  
  
}




################################
#   Plot of Mean of posterior
###################################
spec_est=array(0,c(51,Ntime,NumObs))

for(p in (nwarmup+1):nloop)
{
  cat("p is",p,"\n")
  Tnew=result[[p,1]]$TT
  
  # now we have current partition
  interval_curr=time_interval(x,Tnew$S,Tnew$prt,NumObs)
  
  for(i in 1:Tnew$M)
  {
    prop_obs=x[which(Tnew$prt==i),]
    prop_obs_index=unique(prop_obs$index)
    
    for(j in prop_obs_index)
    {
      spec_est[,interval_curr$tmax[i,j]:interval_curr$tmin[i,j],j]=spec_est[,interval_curr$tmin[i,j]:interval_curr$tmax[i,j],j]+
        repmat(as.matrix(Tnew$fhat_prop[,i]),1,interval_curr$interval[i,j])/(nloop-nwarmup)
    }
    
  }
  
}





# plot estimation of spectrum of the first observation
nfreq=floor(100/2)
freq_hat=(0:nfreq)/100



for(i in 1:NumObs)
{
  plot=plot_ly(x=1:Ntime,y=freq_hat, z = spec_est[ , ,i ]) %>% add_surface() %>% layout(title = paste(
    "Estimated Log Spectrum for data with the", i , "th observation"), 
    scene=list(zaxis=list(title="Log Power Spectrum"),xaxis = list(title = 'Time'), yaxis = list(title = 'Frequency')))
  print(plot)  
  
}








#####################
#    diagnostics
######################
# accept rate
pdf(file= "acceptance_11-19_smooth(seed2333)_prebeta.pdf" )
par(mfrow=c(2,2))
for(i in 1:TT$M)
{
  accept_vec=rep(0,nloop)
  for(p in 1:nloop)
  {
    accept_vec[p]=result[[p,1]]$accepted[i]
    
  }
  
  accept_rate=round(length(which(accept_vec==1))/nloop,2)
  
  plot(accept_vec,ylab="accept",xlab="iteration",main=paste("acceptance for ",i,"th center with ratio ",accept_rate,sep=""))
  
}

dev.off() 



# alpha
pdf(file= "alpha_11-19_smooth(seed2333)_prebeta.pdf" )
par(mfrow=c(2,2))

for(i in 1:TT$M)
{
  
  alpha_vec=rep(0,nloop)
  for(p in 1:nloop)
  {
    alpha_vec[p]=result[[p,1]]$alpha[i]
    
  }
  
  plot(alpha_vec,ylab="alpha",xlab="iteration",main=paste("alpha for ",i,"th center",sep=""))
  
  
}

dev.off()



# plot for beta
pdf(file= "beta_hist_BFGS_11-19_smooth(seed2333))_prebeta.pdf")
for(i in 1:TT$M)
{
  par(mfrow=c(2,2))
  for(j in 1:nbeta)
  {
    
    ave_beta=rep(0,nloop)
    for(p in 1:nloop)
    {
      
      ave_beta[p]=result[[p,1]]$TT$beta[i,j]
      
    }
    
    # mean & sd
    mean=round(mean(ave_beta),2)
    deviance=round(sd(ave_beta),2)
    
    hist(ave_beta,main=paste("center",i,"base",j,"mean",mean,"sd",deviance))
    
    #plot(ave_beta, main = paste("Beta for ",i,"th center and ",j,"th basis",sep = ""))
    
  }
  
}

dev.off() 



# tau 
pdf(file= "tau_11-19_abrupt(seed2333).pdf" )
par(mfrow=c(2,2))
for(i in 1:TT$M)
{
  tau_vec=rep(0,nloop)
  for(p in 1:nloop)
  {
    tau_vec[p]=result[[p,1]]$TT$tau[i]
    
  }
  
  
  plot(tau_vec,ylab="tau",xlab="iteration",main=paste("tau for ",i,"th center ",sep=""))
  
}

dev.off() 



##############################
# sum of squared distances
##############################
# true
ARspec <- function(phi,sigma,freq)
{
  dim = dim(phi)
  len = dim[2]
  phispect = array(0,c(2,2,length(freq)))
  spect = array(0,c(2,2,length(freq)))
  
  for (k in 1:length(freq))
  {
    phispect[,,k] = diag(2)
    for(j in 1:(len/2))
    {
      if(j==1)
      {
        bigmat = phi[,1:2]*exp(-2*pi*sqrt(-1+0i)*freq[k])
      }else{
        bigmat = phi[,(2*j-1):(2*j)]*exp(-2*j*pi*sqrt(-1+0i)*freq[k])
      }
      
      phispect[,,k] = phispect[,,k] - bigmat
      
    }
    spect[,,k] = sigma%*%solve(phispect[,,k])%*%Conj(solve(phispect[,,k]))
    
  }
  
  return(spect)
  
}




spec_true=array(0,c(51,Ntime,NumObs))
ar1=diag(2)
ar2=matrix(0,2,2)
sig=diag(2)
nfreq=floor(100/2)
freq_hat=(0:nfreq)/100


for(l in 1:NumObs)
{
  
  if(xx[l,1]<0.5)
  {
    if(xx[l,2]<0.5)
    {
      
      a=c(repmat(0.3,Ntime/2,1),repmat(-0.5,Ntime/2,1))
      
    }else{
      
      a=c(repmat(-0.9,Ntime/2,1),repmat(0.7,Ntime/2,1))
      
    }
    
    
  }else{
    
    if(xx[l,2]<0.5)
    {
      a=c(repmat(-0.3,Ntime/2,1),repmat(0.5,Ntime/2,1))
      
    }else{
      
      a=c(repmat(0.9,Ntime/2,1),repmat(-0.7,Ntime/2,1))
    }
    
    
  }
  
  # for(t in 1:Ntime)
  # {
  #   spec=ARspec(cbind(a[t]*ar1,ar2),sig,freq_hat)
  #   spec_true[,t,l]=log(Re(spec[1,1,]))
  #   
  # }
}




# periodogram
period=list()
spec_true=list()

for(i in 1:M)
{
  cat("i is",i,"\n")
  
  prop_obs=x[which(prt==i),]
  prop_obs_index=unique(prop_obs$index)
  
  # create log periodograms for replications in ith time seg and jth cov seg
  nseg_time_temp=interval_curr$interval[i,]
  tmin=interval_curr$tmin[i,]
  tmax=interval_curr$tmax[i,]
  
  y <- list()
  yy <- list()
  temp<- list()
  
  uu=0
  
  for(l in prop_obs_index)
  {
    uu=uu+1
    
    nfreq <- floor(nseg_time_temp[l] / 2)
    freq_hat=(0:nfreq)/nseg_time_temp[l]
    
    
    y[[uu]] <- log(abs(fft(x_t[tmin[l]:tmax[l], l])) ^ 2 / nseg_time_temp[l])
    yy[[uu]] <- y[[uu]][1:(nfreq + 1)]
    
    
    for(t in 1:Ntime)
    {
      cat("t is",t,"\n")
      spec=ARspec(cbind(a[t]*ar1,ar2),sig,freq_hat)
      
      
    }
    
    temp[[uu]]=log(Re(spec[1,1,]))
    
    
  }
  
  
  period[[i]]=yy
  spec_true[[i]]=temp
  
  
}


SSE=0
for(i in 1:M)
{
  for(j in 1:5)
  {
    SSE=SSE+sum((period[[i]][[j]]-spec_true[[i]][[j]])^2)
  }
}

SSE

