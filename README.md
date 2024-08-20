# 1 Instuctions

## Installations
The package is installable on any machine that supports R, as long as the following requirements are met:

1. Installation of R or Rstudio 

2. Make sure you have Rtools installed. This could be obtained from CRAN by following the link
https://cran.r-project.org/bin/windows/Rtools/

3. Make sure the R and Rtools paths are added in the environment variable.
                                                                                                               
4. Checks if Rtools is installed by:
      - Installing the R package **devtools**
      - loading the package into R (**library(devtools)**)
      - If the command **find_rtools()** results in **TRUE**, then it should be possible to proceed.
      
5. Installation of the dependent R packages Rcpp, RcppArmadillo, ggplot2 and VGAM (all available on CRAN)

6. The two files to be installed are **pprmcmc_0.1.0.tar.gz** and. The command line **install.packages("path_to_file/pprmcmc_0.1.0.tar.gz", repos = NULL, type = "source",dependencies=TRUE)** where path_to_file would represent the full path will install the package.
   
After the installation the pacakges can be loaded into R/Rstudio using *library(pprmcmc)*. The command line **??pprmcmc** and **vignette("pprmcmc_vignettes")** would provide the package help pages and its documentation respectively. The most important functions are listed in the description page and could be accessed. The list of all functions and some headers could then be accessed by clicking on **Index** seen at the bottom line of the description page.

The source codes (C++ code) are included in the folder **pprmcmc** at /src, the data at /data and any R code is included in /R. Note that the only routine written in R is the speed of propagation **speed.R**.
All the functions (C++) are contained in one file, so any modification of a function regarding the MCMC will be made in **postsamp.cpp** respectively. Note that the examples included are those used to obtain the graphs (regarding PPR) in the main report.

## Modification
The simplest approach for modifying this package would follow these steps:

1. Install the most recent Rstudio

2. Install the packages **devtools, roxygen2, testthat, knitr** (all available on CRAN)

3. Open the file **pprmcmc.Rproj** contained in the folder **pprmcmc**. 

4. Go to **Tools>Project Options>Build Tools>Configure** and check **Build & Reload**

5. After making the requiered modification click **Build & Reload** in the build panel to completely rebuild the package, including updating all the documentation, installs it in your own library.

@WIC15 provides more details of how to build and update such a package. 

# 2 A User’s Guide


## Tunisian Data and initials considerations
We now explain how to use **pprmcmc** on the Tunisia data in @EFSA15 which consist of the times the locations of infected sites. We reflect our vague prior knowledge on the parameters $\alpha$ and $\beta$ by setting $a=c=10^{-4}$ and $ b=d=10^{-6}$ but we choose an informative prior belief for the infectious period by setting $e=0$ and $f=120$. This choice allows premises to remain infectious for a period up to four months as opposed to three months period considered in @EFSA15 since there is no evidence showing an external source of infection for the $4^{th}$ infection which occurred $108$ days after the third one. 

We begin by loading the package and then attaching the PPR data from Tunisia.  Then we provide some initial values for the MCMC routine.
```{r, message=F, warning=F,echo=TRUE}
library("pprmcmc")
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
l=a[a%in%b]
Sourc[i]=sample(l,1)
}
Sourc=Sourc-1
```
## Sample from posterior distributions

**raymcmc**, **expmcmc** and **caumcmc** provide respectively samples from the posteriors distributions using the Rayleigh, the Exponential and the Cauchy kernls. 

*   **Rayleigh kernel**

Considering the Rayleigh kernel and using the default values for the priors and the initial values for the model parameters, (see the package help page or simply **help(raymcmc)**), 100000 iterations using the followng command produced the results in the main report. Note that the algorithm runs using the original parameterisation but presents results for the transformed parameter (see report)

```{r, message=F, warning=F,echo=TRUE, fig.align='center'}
a=raymcmc(Coo, Sourc,Inftim,nits=100000)
```
The result corresponds to a matrix where each column corresponds to the sample from a model parameter. So the posterior distribution of the parameter $\alpha$, $\beta$ and $\gamma$ correspond to the first, second and third column respectively. We discarde the first 10000 iterations and produce the histogram representing the distribution of the each parameter using the popular **ggplot2** package.

```{r, message=F, warning=F,echo=TRUE, fig.align='center'}
b=data.frame(a[10000:100000,])
colnames(b)=c("V1","V2","V3")
library("ggplot2")
q=qplot(1/sqrt(2*V1),data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(50,150))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
```
* **Exponential kernel**

Here, we run the MCMC for 100000 iterations and specify the initials values for the models parameters.
```{r, message=F, warning=F,echo=TRUE, fig.align='center'}
a=expmcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
b=data.frame(a[10000:100000,])
colnames(b)=c("V1","V2","V3")
q=qplot(1/V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(25,250))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
```

* **Cauchy kernel**

Similar to the previous case, we run the MCMC for 100000 iterations and specify the initials values for the models parameters.
```{r, message=F, warning=F,echo=TRUE, fig.align='center'}
a=caumcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
b=data.frame(a[10000:100000,])
colnames(b)=c("V1","V2","V3")
q=qplot(V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(25,350))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
q=qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')
q+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
```

## Posterior distribution of wave speed

The speed of wave speed is estimated with the function **speed**. See **help(speed)** for the function description
The following is the code to poduce the wave speed of propagation using the posterior distribution of the model parameters obtained from the above MCMC. We use the distribution of the parameters included in the package but it is straightforward to directly use the distributions generated by the MCMC. 

We begin by loading the joint posterior parameter distribution into R, randomy drawing 1000 samples from such posteriors, choosing the final observaton time and fixing the frequency at which the system will be observed.

```{r, message=F, warning=F,echo=TRUE}
data(postray)  # Posterior distribution of the model parameters obtained from the MCMC
samp=sample(10000:100000,1000)
alpha=postray[,1][samp]
beta=postray[,2][samp]
Tmax=600  
tim=10:Tmax
gama=postray[,3][samp]
```

The following lines produce the wave speed of propagation using for the Rayleigh kernel.
```{r, message=F, warning=F,echo=TRUE, fig.cap='95%-credible interval for wave speed from posterior predictive simulations of the process', fig.align='center', fig.width=6, fig.height=6}
l=1
speed=speed(alpha,beta,tim,gama,Tmax,samp,l)
# plot speed with 95\% credible interval along with the median  
hist(speed,breaks=50,col='blue',main=' ',xlab='Wave speed in any direction',xaxt='n')
abline(v=round(median(speed),2),col='red',lwd=2)
abline(v=round(quantile(speed,.025),2),col='red',lty=2)
abline(v=round(quantile(speed,.975),2),col='red',lty=2)
axis(side=1,at=c(round(quantile(speed,.025),2),round(median(speed),2),round(quantile(speed,.975),2)),las=2)
box()
```
