% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{raymcmc}
\alias{caumcmc}
\alias{expmcmc}
\alias{raymcmc}
\title{Parameter estimation with MCMC.}
\usage{
raymcmc(Coo, Sourc, Inftim, nits = 1000L, lmin = 0L, lmax = 120L,
  Tmax = 325L, a = 1e-04, b = 1e-05, c = 0.001, d = 1e-04,
  theta = as.numeric(c(1, 2, 110)))

expmcmc(Coo, Sourc, Inftim, nits = 1000L, lmin = 0L, lmax = 120L,
  Tmax = 325L, a = 1e-04, b = 1e-05, c = 0.001, d = 1e-04,
  theta = as.numeric(c(1, 2, 110)))

caumcmc(Coo, Sourc, Inftim, nits = 1000L, lmin = 0L, lmax = 120L,
  Tmax = 325L, a = 1e-04, b = 1e-05, c = 0.001, d = 1e-04,
  theta = as.numeric(c(1, 2, 110)))
}
\arguments{
\item{Coo}{A NX2 matrix indicating the location (longitude latitude) of infected hosts with N the size of the
infected population.}

\item{Sourc, Inftim}{Vector of length N indicating the initial source of infection of every host
with -1 for the first infection and the corresponding infection times respectively.}

\item{nits, lmin, lmax, Tmax}{Indicate the number of iterations (10000 by default), the lower and upper bound of the prior
infectious period (0 and 110 by defaut) and the final observation time (325 by default)  respectively.}

\item{a, b, c, d}{The gamma conjugate priors for the parametes. Vague priors are considere as default.}

\item{theta}{A vector of size 3 indicating the initial values for the parameters alpha, beta and gamma.}
}
\value{
A nitsX3 matrix indicating sample from the posterior distributions of alpha, beta and gamma.
}
\description{
Sample from the posterior distribution of the contact rate, the kernel parameter and
the infectious period using rayleigh (\code{raymcmc}), exponential (\code{expmcmc}) and cauchy (\code{caumcmc}) kernels.
}
\details{
This is a Markov Chain Monte Carlo (MCMC) technique to sample from the posterior distribution
of the model parameters. It uses Gibbs sampling \url{https://en.wikipedia.org/wiki/Gibbs_sampling}
to sample from the model parametes and use data augmentation technique to change in the
transmission network.
}
\examples{
data(ppr)
attach(ppr)
#Premisses locations
Coo=matrix(c(ppr$Longitude,ppr$Latitude),nrow=26,ncol=2,byrow=F)

#Infection times
diff=as.Date(Date.of.start.of.the.outbreak)-as.Date(Date.of.start.of.the.outbreak)[1]
inf=as.numeric(diff)
inf[2:26]=inf[2:26]+runif(25)
Inftim=sort(inf)

#Initial sources
infper=110
rt=Inftim + infper
Sourc=numeric(26)
for(i in 2:26){
a=which(Inftim[i]>Inftim)
b=which(Inftim[i]<rt)
l=a[a\%in\%b]
Sourc[i]=sample(l,1)
}
Sourc=Sourc-1

#Rayleigh kernel
a=raymcmc(Coo, Sourc,Inftim)
a=raymcmc(Coo, Sourc,Inftim,nits=100000)
#a=caumcmc(Coo, Sourc,Inftim,nits=10000,theta=c(.1,.3,infper))
b=data.frame(a[10000:100000,])
colnames(b)=c("V1","V2","V3")
require("ggplot2")
qplot(1/sqrt(2*V1),data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(50,150))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")

#Exponential kernel
a=expmcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
colnames(b)=c("V1","V2","V3")
qplot(1/V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(25,250))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")

#Cauchy kernel
a=caumcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
b=data.frame(a[10000:100000,])
colnames(b)=c("V1","V2","V3")
qplot(V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(40,350))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0.0025,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
}

