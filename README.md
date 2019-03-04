# pprmcmc
Markov Chain Monte Carlo (MCMC) implementation for the Bayesian network model
# 1 Introduction
The problem of analysis of a network diffusion in general is challenging, given that the transmission network is often hidden. This complexity is exacerbated for livestock diseases such as Pest des Petits Ruminants (PPR) since the locations of the nodes are not known in advance.  The location may only become known when a node is infected [@EFSA15]. Techniques relevant to modelling PPR spread within a country should therefore be capable of handling the problem of unobserved susceptible populations when estimating system parameters.  There is a large body of literature that uses machine-learning methods for estimating the transmission networks e.g [@EPL09; @KSAX10;@LOSI10;@MYLE10; @GRLEKR10; @GRBASC11; @DSWZ13]. These approaches to estimation of transmission networks consider only the temporal dynamic and do not take account of spatial aspects. Therefore, to capture within-country spread we have opted to adapt explicitly spatio-temporal models exemplified by Paper citep{CMBG07}, and [@CGC08; @SSF09 @Par09; @DCMGB11; @CLNDG14; @PGPG14; @CSDGG15] utilising modern Bayesian computational techniques in developing bespoke fitting algorithms. These are illustrated using the data from Tunisia. 

**pprmcmc** performs the popular Markov Chain Monte Carlo (MCMC) coupled with the Bayesian data augmentation techinques by considering the missing transmission network as addional parameters. In addition to generating the joint posterior distribution of the model parameters and the transmission network, it provides users with a simple way to generate simulate replicates of the model process and to estimate the wave speed of propagation.

# 2 What **pprmcmc** does
MCMC and Bayesian data augmentation in this context involves imputing, at each iteration of the MCMC, a transimission network so that the likelihood becomes tractable. The model parameters are then drawn considering their full conditional distributions. **pprmcmc** then provides samples of the model parameters and the transmission network which are treated as their posterior distributions. Under normal circumstances, you only need to provide some initials parameters including the models parameters and the transimission network and can then use the posterior distribution for estimating the wave speed of propagation. The posterior distribution of the wave speed is the provided by **pprmcmc** by simpling providing the distribution of the model parameters.

The advantage of **pprmcmc** is that it combines the comparative speed and ease-of-use of our algorithm with the power of transmission network imputation.

## Assumptions 
We represent the spread of PPR through a landscape using a model that can be fitted to data where only infections are recorded.

The main assumptions are as follows:

1. An infection premises remain infectious for fixed period, $\gamma$

2. Premises cause new infections according to Poisson process with fixed rate $\beta$ and the new infection is located at random distance $r\sim f(r;\alpha)$ from the source.

3. The unknown susceptible population is implicitly modelled in this approach as a spatial Poisson process.

We note that the representation of the susceptible population means that the effect of depletion of the susceptible population is not explicitly represented in the population.  However, the approach is valid for the linear phase of the epidemic and, arguably, is appropriate for the estimation of speeds of the disease ‘wave-front'.

For a given choice of distribution f, there are 3 key parameters in this model – the contact rate $\beta$, the infectious period $\gamma$ and the parameter $\alpha$ which controls the spatial dispersal process.  All of these play a role in the estimation of rates of spatial spread.   We have proposed a number of candidate models: 
\begin{align}\label{kernel}
\textbf{Rayleigh($\alpha$), i.e} \qquad f(r;\alpha)=&2r\alpha\exp\left(-\alpha r^2\right)\nonumber\\
\textbf{Exponential($\alpha$), i.e }\qquad f(r;\alpha)=&\alpha\exp\left(-\alpha r\right)\\
\textbf{Cauchy($\alpha$), i.e                   }\qquad  f(r;\alpha)=&\frac{2}{\pi \alpha\left(1+\frac{r^2}{\alpha^2}\right)}
\end{align}

in which $r$ is the Euclidean distance between a given pair of premises and measured in $km$. Note that these one dimensional kernels are isotropic defined on the positive real axis. These models represent scenarios in which the propensity for long range interactions increases – with the Cauchy model exhibiting ‘jumps’ in the infection process that are very long-range. The kernels defined above ensure that different pattern of the epidemic could be compared directly. While thin-tailed such as Rayleigh and exponential kernels give rise to slow spreading wave of new infections, the thick-tailed kernels such as Cauchy result in a rapid and long dispersal ahead of the source [@MOL77; @CSKM99]. Techniques developed here can readily be adapted to other models in addition to the three variants considered above.

Analyses under **pprmcmc** were carried out assuming that the infectious period $\gamma$ was known, but this assumption is relaxed in the present pacakage where $\gamma$ is inferred.


## Inference methods

### Bayesian imputation of the network
