#' Simulaton of the epidemic process to find the speed of the disease
#'
#' Epidemic simulation and maximum distance the wave can travel using gillespsie algorithm.
#'
#'\code{simul} provide the simulation of the epidemic process and the
#'       maximum distance the wave can travel at each particular obaservation date.
#'
#' @param alpha, beta, gama indicating the dispersal kernel paremeter,
#'       the contact parameter and the infectious period.
#' @param tim  A vector of observation times
#' @param Tmax  Final observation time
#' @param l Takes values 1, 2, 3.
#'        \enumerate{
#'        \item indicates the rayleigh kernerl
#'        \item indicates the exponential kernel
#'        \item indicates the cauchy kernel
#'        }
#' @return A list with components:
#'       \describe{
#'         \item{epidem}{A five-dimentional matrix giving the  dynamic of the process .
#'        Each column indicates respectively the times, the x cooedinate, the y-coordinate
#'        indicator (0 if infection and 1 if removal) and the index of individual
#'        removed (0 if the event is infection).}
#'         \item{maxdist}{A vector of maximum distances travelled by the wave on each time in tim.}
#'        }
#' @import VGAM
#' @examples
#'# Simulation with rayleigh kernel
#' alpha=0.00012
#' beta=0.012
#' tim=10:325
#' gama=110
#' Tmax=325
#' l=1
#' res=simul(alpha,beta,tim,gama,Tmax,l)
#'

#' @export
simul=function(alpha,beta,tim,gama,Tmax,l){
  t=0
  inf=0
  i=1
  S=0
  size=0
  dat=c(0,0,0,1,0)
  init=c(0,0,0)
  inf_lis=1
  ni=1
  indx=0
  ii=0
  remt=rem_tim=gama
  rp=0
  vec=vecx=vecy=NULL
  while(t<=Tmax){
    if(l==1){
      r=rrayleigh(1,1/sqrt(2*alpha))
    }
    if(l==2){
      r=rexp(1,1/alpha)
    }
    if(l==3){
      r=rcauchy(1,alpha)
    }
    theta=runif(1,0,2*pi)
    rate=ni*beta
    dt=log(1/runif(1))/rate
    if(length(inf_lis)==1){

      while(min(rem_tim)<dt){ #only one infection left thus has to be an unfection
        dt=log(1/runif(1))/rate
      }
      j=inf_lis[1]
      size=c(size,ni)

      if(ii==0){
        S=c(S,j)
        x=dat[2]+ r*cos(theta)
        y=dat[3]+ r*sin(theta)
        ii=1
      }
      else{
        S=c(S,j)
        x=dat[j,2]+ r*cos(theta)
        y=dat[j,3]+ r*sin(theta)
      }
      ni=ni+1
      dat=rbind(dat,c(t+dt,x,y,ni,0))
      inf_lis=c(inf_lis,nrow(dat))
      indx=c(indx,nrow(dat)-1)
      rem_tim=rem_tim-dt
      rem_tim=c(rem_tim,gama)
      remt=c(remt,gama)
      dd=sqrt((x-init[1])^2+(y-init[2])^2)
      if(dd>rp){
        rp=dd
      }
    }
    else{
      if(min(rem_tim)<dt){  #removal
       dt=min(rem_tim)
        ni=ni-1
        j=which(rem_tim==min(rem_tim))
        dat=rbind(dat,c(t+dt,0,0,ni,inf_lis[j]))
        inf_lis=inf_lis[-j]
        rem_tim=rem_tim[-j]-rem_tim[j]
      }
      else{  # infection
        ni=ni+1
        j=sample(inf_lis,1)
        x=dat[j,2]+ r*cos(theta)
        y=dat[j,3]+ r*sin(theta)
        dat=rbind(dat,c(t+dt,x,y,ni,0))
        inf_lis=c(inf_lis,nrow(dat))
        indx=c(indx,nrow(dat)-1)
        rem_tim=rem_tim-dt
        rem_tim=c(rem_tim,gama)
        remt=c(remt,gama)
        size=c(size,ni-1)
        S=c(S,j)
        dd=sqrt((x-init[1])^2+(y-init[2])^2)
        if(dd>rp){
          rp=dd
        }
      }
    }
    if(length(tim)!=0){
      if(t+dt>Tmax){
        vec=c(vec,rep(rp,length(tim)))
      }
      else{
        while(tim[1]>t&&tim[1]<=(t+dt)){
          if(t==0){
            vec=c(vec,0)
          }
          else{
            vec=c(vec,rp)

          }
          tim=tim[-1]
        }
      }

    }
    if(t<Tmax && t+dt>Tmax){
      dat[nrow(dat),]=c(t+dt,0,0,ni,0)
      break
    }
    t=t+dt


  }
  return(list(maxdist=vec, epidem=dat))
}

#' Posterior distribution of Wave speed
#'
#'\code{speed} provide sample from posterior wave speed using the posterior
#'             distribution of the model parameters.
#'
#' @param alpha, beta, gama indicating sample from the posterior distribution
#'        of the dispersal kernel paremeter, the contact parameter and the
#'        infectious period.
#' @param tim  A vector of observation times i.e. the equence of times at which the system progress is observed.
#' @param Tmax  Final observation time
#' @param samp  A vector of sample to draw from the posteriors distributions.
#' @param l Takes values 1, 2, 3.
#'        \enumerate{
#'        \item indicates the rayleigh kernerl
#'        \item indicates the exponential kernel
#'        \item indicates the cauchy kernel
#'        }
#'
#' @details The parameterisation of  distributions used are:
#'           \enumerate{
#'        \item Rayleigh kernerl: \eqn{f(r;\alpha)=2\alpha*r*exp(\alpha*r^2)}
#'        \item Exponential kernel: \eqn{f(r;\alpha)=\alpha*exp(\alpha*r)}
#'        \item Cauchy kernel:  \eqn{f(r;\alpha)=1/(\pi\alpha(1 + r^2/\alpha^2))}
#'        }
#'
#' @return A vector indicating the speed of progation.
#'@examples
#'# Simulation with Rayleigh kernel
#' data(postray)  # Posterior distribution of the model parameters obtained from the MCMC
#' samp=sample(10000:100000,1000)
#' alpha=postray[,1][samp]
#' beta=postray[,2][samp]
#' Tmax=325  # Increase Tmax to 600 for ex for more sample and better estimate of the speed.
#' tim=10:Tmax
#' gama=postray[,3][samp]
#' l=1
#' speed=speed(alpha,beta,tim,gama,Tmax,samp,l)
#' # plot speed with 95\% credible interval along with the median
#' hist(speed,breaks=50,col='blue',main=' ',xlab='Wave speed in any direction',xaxt='n')
#' abline(v=round(median(speed),2),col='red',lwd=2)
#' abline(v=round(quantile(speed,.025),2),col='red',lty=2)
#' abline(v=round(quantile(speed,.975),2),col='red',lty=2)
#' axis(side=1,at=c(round(quantile(speed,.025),2),round(median(speed),2),round(quantile(speed,.975),2)),las=2)
#' box()
#' @export
speed=function(alpha,beta,tim,gama,Tmax,samp,l){
  Mat=array(0,c(length(samp),length(tim)))
  for(i in 1:length(alpha)){
    vec=simul(alpha[i],beta[i],tim,gama[i],Tmax,l)
    Mat[i,]=vec[[1]]
  }
  speed=apply(Mat,2,max)/tim
  return(speed)
}

