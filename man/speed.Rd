% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speeddist.R
\name{speed}
\alias{speed}
\title{Posterior distribution of Wave speed}
\usage{
speed(alpha, beta, tim, gama, Tmax, samp, l)
}
\arguments{
\item{alpha, }{beta, gama indicating sample from the posterior distribution
of the dispersal kernel paremeter, the contact parameter and the
infectious period.}

\item{tim}{A vector of observation times i.e. the equence of times at which the system progress is observed.}

\item{Tmax}{Final observation time}

\item{samp}{A vector of sample to draw from the posteriors distributions.}

\item{l}{Takes values 1, 2, 3.
\enumerate{
\item indicates the rayleigh kernerl
\item indicates the exponential kernel
\item indicates the cauchy kernel
}}
}
\value{
A vector indicating the speed of progation.
}
\description{
\code{speed} provide sample from posterior wave speed using the posterior
            distribution of the model parameters.
}
\details{
The parameterisation of  distributions used are:
          \enumerate{
       \item Rayleigh kernerl: \eqn{f(r;\alpha)=2\alpha*r*exp(\alpha*r^2)}
       \item Exponential kernel: \eqn{f(r;\alpha)=\alpha*exp(\alpha*r)}
       \item Cauchy kernel:  \eqn{f(r;\alpha)=1/(\pi\alpha(1 + r^2/\alpha^2))}
       }
}
\examples{
# Simulation with Rayleigh kernel
data(postray)  # Posterior distribution of the model parameters obtained from the MCMC
samp=sample(10000:100000,1000)
alpha=postray[,1][samp]
beta=postray[,2][samp]
Tmax=325  # Increase Tmax to 600 for ex for more sample and better estimate of the speed.
tim=10:Tmax
gama=postray[,3][samp]
l=1
speed=speed(alpha,beta,tim,gama,Tmax,samp,l)
# plot speed with 95\\\% credible interval along with the median
hist(speed,breaks=50,col='blue',main=' ',xlab='Wave speed in any direction',xaxt='n')
abline(v=round(median(speed),2),col='red',lwd=2)
abline(v=round(quantile(speed,.025),2),col='red',lty=2)
abline(v=round(quantile(speed,.975),2),col='red',lty=2)
axis(side=1,at=c(round(quantile(speed,.025),2),round(median(speed),2),round(quantile(speed,.975),2)),las=2)
box()
}

