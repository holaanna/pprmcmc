
#include <RcppArmadilloExtensions/sample.h>


// #include <Rcpp.h>
 using namespace std;
 //using namespace arma;
 using namespace Rcpp;
//// [[Rcpp::interfaces(r, cpp)]]

/*-------------------------- Interagate beta ---------------------------------------------*/
// [[Rcpp::export]]
double integ(NumericMatrix dyn, double Tmax,int N){
  int k=1;
  double sum=0;
  while(k<=N){
    if(k==N-1 && dyn(N-1,0)<Tmax){
      sum+=dyn(k-1,3)*(dyn(k,0)-dyn(k-1,0)+(Tmax-dyn(k-1,0)));
      break;

    }
    if(k==N-1 && dyn(N-1,0)>Tmax){
      sum+=dyn(k-1,3)*(Tmax-dyn(k-1,0));
      break;

    }
    sum+=dyn(k-1,3)*(dyn(k,0)-dyn(k-1,0));

    k++;
  }
  return(sum);
}

/*----------------------------------Distance between two locations----------------------------*/
//' Distance between two locations.
//'
//'\code{distance} returns the square of the distance between two locations.
//'
//' @param phi1,phi2  Longititudes of location 1 and 2 respectively.
//' @param lamda1,lamda2  Latitudes of location 1 and 2 respectively.
//'
//' @return The distance between locations \code{phi1,lamda1} and \code{phi2,lamda2}.
//'
//' @examples
//' distance(10,20,5,10)
//' @export
// [[Rcpp::export]]
double distance(double phi1,double phi2,double lamda1,double lamda2){
  const double pi = 3.1415926535897;
  double R=6371;
  double x1=2*R*tan(pi/4-phi1*pi/360)*sin(lamda1*pi/180);
  double y1=-2*R*tan(pi/4-phi1*pi/360)*cos(lamda1*pi/180);
  double x2=2*R*tan(pi/4-phi2*pi/360)*sin(lamda2*pi/180);
  double y2=-2*R*tan(pi/4-phi2*pi/360)*cos(lamda2*pi/180);
  return((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

/*------------------------Sample  a source of infection--------------------*/

//' Sample  a source of infection to the current premisse.
//'
//'\code{Soursamp} returns the index of the source that infect a premisse.
//'
//' @param rtime,time Vector of removals and infections times respectively.
//' @param ind,size,ele Integers indicating the index of the premisse in the vector time,
//'        the size of the potential sources and the premisse selected respectively.
//'
//' @return The index of the source of the premisse at the position \code{ele}.
//' @export
// [[Rcpp::export]]
int Soursamp(NumericVector rtime, NumericVector time, int ind, int siz,int ele){
  IntegerVector S(siz);
  int k=0;
  for(int j=0;j<ind;j++){

    if((time[ind]>time[j])&& (rtime[j]>time[ind])){
      S[k]=j;
      k++;

    }

  }
  return(S[ele]);
}

/*-----------------------Index of an element in a vector-------------*/
// [[Rcpp::export]]
int Which(SEXP rhs, double key){
  NumericVector vec(rhs);
  int i=0;
  for (i=0;i<vec.size();i++) {
    if(vec[i]==key){
      break;
    }
  }
  return i;
}
//------------------------------------------------------------------------------------


//' Check for existence of an element in a vector.
//'
//'\code{contains} returns a boolean either true or false.
//'
//' @param x A vector.
//' @param z An element to be checked.
//'
//' @return True if \code{z} is in vector \code{X} and False otherwise.
//'
//' @examples
//' distance(10,20,5,10)
//' @export
// [[Rcpp::export]]
bool contains(SEXP X, double z) {
  NumericVector vec(X);
  return std::find(vec.begin(), vec.end(), z)!=vec.end();
}


/*-----------------------Construction of the process-----------------------------------*/

//' Construction of the epidemic process.
//'
//'\code{construct} reconvers the dynamic of the epidemic process. It reorders the events
//'                 as they happened'
//'
//' @param Coo Matrix indicating the locations of all premisses (Longitude,latitude).
//' @param Remtime,Inftime Vectors indicating their corresponding removal and infection
//'        times respectively .
//' @param infper The infectious period.
//' @param NR,Tmax The number of hosts removed and the final observation time.
//' @param M1 2NX5 Matrix to store the dynamic of the process where colomn:
//'       \enumerate{
//'       \item  Contains the events times
//'       \item  Contains the latitude of the infected premisse.
//'       \item  Contains the longititu of the infected premisse.
//'       \item  Indicates the type of event, particulary we o for infection else the
//'              index of host removed.
//'       \item  The trajectory of the epidemic ie the size of the infections.
//'       }
//' @param Indx,size Empty vectors to store the position of premisses in the matrix M1 that
//'        store the process.
//' @param nrows Integer to keep the number of rows used in M1 since there might be more
//'        rows than events.
//'
//' @return The dynamic of the epidemic process, the position of each infection, the size
//'         of the source of infection of each infected premisse and the number of events.
//'
//' @export
// [[Rcpp::export]]
void construct(NumericMatrix Coo, NumericVector Inftime,double infper,int NR,int Tmax,IntegerVector &Indx,IntegerVector &size,  NumericMatrix &M1, int &nrows ){

  int NI=Inftime.size();
  NumericVector Remtime(Inftime + infper);
  int i,k,k1; int ni=0;
  NumericVector tsort= union_(Inftime,Remtime);
  std::sort(tsort.begin(),tsort.end()); // Sort the events
  //
  for(i=0;i<(NI+NR);i++){
    if(tsort[i]>Tmax){
      --i;
      break;
    }
    if(std::find(Inftime.begin(), Inftime.end(), tsort[i]) != Inftime.end() ) {
      k=Which(Inftime,tsort[i]);
      k1=0;
      ni++;
      Indx[k]=i;
      size[k]=ni-1;
      M1(i,0)=tsort[i];
      M1(i,1)=Coo(k,0);
      M1(i,2)=Coo(k,1);
      M1(i,4)=0;
      M1(i,3)=ni;
    } else {  // removal
      k=Which(Remtime,tsort[i]);
      k1=k;
      --ni;
      M1(i,0)=tsort[i];
      M1(i,1)=0;
      M1(i,2)=0;
      M1(i,4)=k1+1;
      M1(i,3)=ni;
    }


  }
  nrows=i+1;
}

//' @rdname construct
//' @export
// [[Rcpp::export]]
void construct1(NumericMatrix Coo, NumericVector Inftime,NumericVector Remtime,int NR,int Tmax,IntegerVector &Indx,IntegerVector &size,  NumericMatrix &M1, int &nrows ){

  int NI=Inftime.size();
  int i,k,k1; int ni=0;
  NumericVector tsort= union_(Inftime,Remtime);
  std::sort(tsort.begin(),tsort.end()); // Sort the events
  //
  for(i=0;i<(NI+NR);i++){

    if(tsort[i]>Tmax){
      i--;
      break;
    }
    if(std::find(Inftime.begin(), Inftime.end(), tsort[i]) != Inftime.end() ) {
      k=Which(Inftime,tsort[i]);
      k1=0;
      ni++;
      Indx[k]=i;
      size[k]=ni-1;
      M1(i,0)=tsort[i];
      M1(i,1)=Coo(k,0);
      M1(i,2)=Coo(k,1);
      M1(i,4)=0;
      M1(i,3)=ni;
    } else {  // removal
      k=Which(Remtime,tsort[i]);
      k1=k;
      --ni;
      M1(i,0)=tsort[i];
      M1(i,1)=0;
      M1(i,2)=0;
      M1(i,4)=k1+1;
      M1(i,3)=ni;
    }


  }
  //Rcpp::Rcout<<i<<std::endl;
  nrows=i+1;
}

/*-------------------------------------Likelihood---------------------------------*/

//' Likelihood computation using rayleigh distribution as kernel.
//'
//'\code{likelihood} compute the log-likelihood of the parameters when the kernl is
//' the rayleigh distribution.
//'
//' @param Coo,Inftime Vectors indicating the locations of all premisses and their corresponding
//'        infection times respectively .
//' @param dyn 2NX5 Matrix indicating dynamic of the process where colomn:
//'       \enumerate{
//'       \item  Contains the events times
//'       \item  Contains the latitude of the infected premisse.
//'       \item  Contains the longititu of the infected premisse.
//'       \item  Indicates the type of event, particulary we o for infection else the
//'              index of host removed.
//'       \item  The trajectory of the epidemic ie the size of the infections.
//'       }
//' @param Sou_siz,S Vectors indicating the potential source size, the position of premisses
//'        in the matrix M1 that store the process.
//' @param beta,alpha The contact and the dispersal parameters.
//' @param nrows Integer indicating the number of rows used in dyn since there might be more
//'        rows than events.
//' @return The log-likelihood of the parameters.
//' @export
// [[Rcpp::export]]
double likelihood(NumericMatrix Coo,IntegerVector Sou_siz,IntegerVector S, NumericMatrix dyn,double beta,double alpha,int nrows){
  double sum=0;
  // int  N=Coo.nrow();
  for(int i=1;i<nrows;i++){
    sum+=-beta*dyn(i-1,3)*(dyn(i,0)-dyn(i-1,0));
  }

  sum-=Rcpp::sum(Rcpp::log(Sou_siz));
  // Rcpp::Rcout<<sum<<' '<<Rcpp::sum(Rcpp::log(Sou_siz))<<std::endl;
  // for(int i=1;i<N;i++){
  //  // double r=distance(dyn(S(i),1),dyn(i,1),dyn(S(i),0),dyn(i,0));
  //   sum+=-log(Sou_siz[i]);//-alpha*r;
  // }
  //return(sum+N*log(beta)+N*log(2*alpha));
  return(sum);
}
//---------------------------------------------------------------------------------------
// [[Rcpp::export]]
double lik(double alpha,IntegerVector S,NumericMatrix mat){
  double sum=0;
  for(int i=0;i<(S.size()-1);i++){
    double dist2=mat(std::min(S[i+1],i+1),std::max(S[i+1],i+1));
    sum+=log(alpha*(1+dist2/(alpha*alpha)));

  }
  return sum;
}

/*-----------------------------------MCMC with rayleigh kernel---------------------------------------------*/

//' Parameter estimation with MCMC.
//'
//' Sample from the posterior distribution of the contact rate, the kernel parameter and
//' the infectious period using rayleigh (\code{raymcmc}), exponential (\code{expmcmc}) and cauchy (\code{caumcmc}) kernels.
//'
//' This is a Markov Chain Monte Carlo (MCMC) technique to sample from the posterior distribution
//' of the model parameters. It uses Gibbs sampling \url{https://en.wikipedia.org/wiki/Gibbs_sampling}
//' to sample from the model parametes and use data augmentation technique to change in the
//' transmission network.
//'
//' @param Coo A NX2 matrix indicating the location (longitude latitude) of infected hosts with N the size of the
//'        infected population.
//' @param Sourc,Inftim  Vector of length N indicating the initial source of infection of every host
//'        with -1 for the first infection and the corresponding infection times respectively.
//' @param nits,lmin,lmax,Tmax Indicate the number of iterations (10000 by default), the lower and upper bound of the prior
//'        infectious period (0 and 110 by defaut) and the final observation time (325 by default)  respectively.
//' @param a,b,c,d The gamma conjugate priors for the parametes. Vague priors are considere as default.
//' @param theta   A vector of size 3 indicating the initial values for the parameters alpha, beta and gamma.
//' @return A nitsX3 matrix indicating sample from the posterior distributions of alpha, beta and gamma.
//'
//' @import ggplot2
//'
//' @examples
//' data(ppr)
//' attach(ppr)
//' #Premisses locations
//' Coo=matrix(c(ppr$Longitude,ppr$Latitude),nrow=26,ncol=2,byrow=F)
//'
//' #Infection times
//' diff=as.Date(Date.of.start.of.the.outbreak)-as.Date(Date.of.start.of.the.outbreak)[1]
//' inf=as.numeric(diff)
//' inf[2:26]=inf[2:26]+runif(25)
//' Inftim=sort(inf)
//'
//' #Initial sources
//' infper=110
//' rt=Inftim + infper
//' Sourc=numeric(26)
//' for(i in 2:26){
//' a=which(Inftim[i]>Inftim)
//' b=which(Inftim[i]<rt)
//' l=a[a%in%b]
//' Sourc[i]=sample(l,1)
//' }
//' Sourc=Sourc-1
//'
//'#Rayleigh kernel
//' a=raymcmc(Coo, Sourc,Inftim)
//' a=raymcmc(Coo, Sourc,Inftim,nits=100000)
//' #a=caumcmc(Coo, Sourc,Inftim,nits=10000,theta=c(.1,.3,infper))
//' b=data.frame(a[10000:100000,])
//' colnames(b)=c("V1","V2","V3")
//' require("ggplot2")
//' qplot(1/sqrt(2*V1),data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(50,150))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//'
//'#Exponential kernel
//' a=expmcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
//' colnames(b)=c("V1","V2","V3")
//' qplot(1/V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(25,250))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//'
//'#Cauchy kernel
//' a=caumcmc(Coo, Sourc,Inftim,nits=100000,theta=c(.1,.3,infper))
//' b=data.frame(a[10000:100000,])
//' colnames(b)=c("V1","V2","V3")
//' qplot(V1,data=b,geom='histogram',xlab=bquote(alpha~(km)),ylab='',xlim=c(40,350))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V2,data=b,geom='histogram',xlab=bquote(beta~(days^-1~km)),ylab='',xlim=c(0.0025,.02))+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' qplot(V3,data=b,geom='histogram',xlab=bquote('Infectious period'~(days^-1)),ylab='')+ geom_histogram(aes(fill = ..count..)) +theme(legend.position="none")
//' @export
// [[Rcpp::export]]
NumericMatrix raymcmc(NumericMatrix Coo, IntegerVector Sourc,NumericVector Inftim, int nits=1000,int lmin=0,int lmax=120,int Tmax=325,double a=1e-4, double b=10e-6,
                      double c=1e-3, double d=1e-4, NumericVector theta=NumericVector::create(1.0,2.0,110)){

  int N=Inftim.size(); int nrows;
  // Define a large matrix to store all events
  //Initial data construction
  NumericMatrix dat(2*N,5);
  IntegerVector Sou_siz(N);
  IntegerVector Indx(N);
  NumericVector rt=Inftim + theta(2);
  construct(Coo, Inftim, theta(2),N,Tmax,Indx,Sou_siz,dat,nrows);
  Sou_siz[0]=1;

  double pacc;
  int acc=1;
  //Define the outputs
  NumericMatrix Par(nits,3);   // Paremters including alpha,beta,lamda(infectious period)
  NumericMatrix Sour(nits,N);  // Source of each infection

  double rtnew;
  double sum=0;
  double rtim=theta[2];
  NumericMatrix mat(Indx.size(),Indx.size()); // Interdistances
  for(int i=0;i<(Indx.size());i++){   // computation of the distance between an infection site and its source
    for(int j=0;j<Indx.size();j++){
      if(i<j){
        mat(i,j)=distance(Coo(i,1),Coo(j,1),Coo(i,0),Coo(j,0)) ;
      }

    }


  }
  for(int i=0;i<(Indx.size()-1);i++){
    sum+=mat(std::min(Sourc(i+1),i+1),std::max(Sourc(i+1),i+1));
  }


  Par(0,_)=theta;

  for(int j=1;j<N;j++){
    Sour(0,j)=Sourc[j]+1;
  }

  double sum1=integ(dat,Tmax,nrows); // initial integration in beta
  int i; int nrows1;

  IntegerVector Sou_siz1= Rcpp::clone(Sou_siz); IntegerVector Indx1= Rcpp::clone(Indx);

  //--------------------
  for( i=1;i<nits;i++){
    Par(i,1)= Rcpp::rgamma(2,N-1+c,1/(sum1+d))[0];
    //update Alpha

    Par(i,0)= Rcpp::rgamma(2,N-1+a,1/(sum+b))[0];

    //Update the infection times
    int k=1;

    while(k<N){

      //Symptomatic individual
      int ll=k;  // choose site seqentially but with larger data can scan choise ie randomly chosen site
      IntegerVector samp=seq_len(Sou_siz[ll]);
      RNGScope scope;
      int j= Rcpp::RcppArmadillo::sample(samp,1,false)[0]-1;
      int ii=Soursamp(rt,Inftim,ll,Sou_siz[ll],j);
      double rnew=mat(std::min(ii,ll),std::max(ii,ll));
      double rold=mat(std::min(Sourc[k],ll),std::max(Sourc[k],ll));
      //Rcpp::Rcout<<rnew<<" "<<rold<<std::endl;
      pacc=exp(-Par(i,0)*(rnew-rold));

      if(pacc> R::runif(0,1)){
        Sourc[k]=ii;
        sum=sum-rold+rnew;
        if(k==10){
          acc+=1;
        }

      }

      k++;
    }
    //Update to the infectious period
    rtnew=rtim + R::runif(-2,2);
    int j;
    for(j=1;j<N;j++){
      if(Inftim[j]>(Inftim[Sourc[j]]+rtnew)){
        break;
      }
    }
    if(rtnew>0 && rtnew<lmax && j==(N)){
      //   Rcpp::Rcout<<rtnew<<std::endl;
      construct(Coo, Inftim, rtnew,N,Tmax,Indx1,Sou_siz1,dat,nrows1);
      Sou_siz1[0]=1;

     // if(is_false(any(Sou_siz1==0))){  // no disjoint chain
        double sum_new=integ(dat,Tmax,nrows1);
       // double product = std::accumulate(Sou_siz.begin(), Sou_siz.end(), 1.0, std::multiplies<double>());
       // double product1 = std::accumulate(Sou_siz1.begin(), Sou_siz1.end(), 1.0, std::multiplies<double>());

        pacc=exp(-Par(i,1)*(sum_new-sum1));
        if(pacc>R::runif(0,1)){
          rtim=rtnew;
          sum1=sum_new;
          nrows=nrows1;
          for(int ll=0;ll<N;ll++){
            Sou_siz[ll]=Sou_siz1[ll];
            Indx[ll]=Indx1[ll];
            rt[ll]=Inftim[ll]+rtim;
          }
        }
     // }
    }

    Par(i,2)=rtim;
    for(int j=1;j<N;j++){
      Sour(i,j)=Sourc[j]+1;
    }
  }

  return Par;


}

//' @rdname raymcmc
//' @export
// [[Rcpp::export]]
NumericMatrix expmcmc(NumericMatrix Coo, IntegerVector Sourc,NumericVector Inftim, int nits=1000,int lmin=0,int lmax=120,int Tmax=325,double a=1e-4, double b=10e-6,
                      double c=1e-3, double d=1e-4, NumericVector theta=NumericVector::create(1.0,2.0,110)){

  int N=Inftim.size(); int nrows;
  // Define a large matrix to store all events
  //Initial data construction
  NumericMatrix dat(2*N,5);
  IntegerVector Sou_siz(N);
  IntegerVector Indx(N);
  NumericVector rt=Inftim + theta(2);
  construct(Coo, Inftim, theta(2),N,Tmax,Indx,Sou_siz,dat,nrows);
  Sou_siz[0]=1;

  double pacc;
  int acc=1;
  //Define the outputs
  NumericMatrix Par(nits,3);   // Paremters including alpha,beta,lamda(infectious period)
  NumericMatrix Sour(nits,N);  // Source of each infection

  double rtnew;
  double sum=0;
  double rtim=theta[2];
  NumericMatrix mat(Indx.size(),Indx.size()); // Interdistances
  for(int i=0;i<(Indx.size());i++){   // computation of the distance between an infection site and its source
    for(int j=0;j<Indx.size();j++){
      if(i<j){
        mat(i,j)=sqrt(distance(Coo(i,1),Coo(j,1),Coo(i,0),Coo(j,0)));
      }

    }


  }
  for(int i=0;i<(Indx.size()-1);i++){
    sum+=mat(std::min(Sourc(i+1),i+1),std::max(Sourc(i+1),i+1));
  }


  Par(0,_)=theta;

  for(int j=1;j<N;j++){
    Sour(0,j)=Sourc[j]+1;
  }

  double sum1=integ(dat,Tmax,nrows); // initial integration in beta
  int i; int nrows1;

  IntegerVector Sou_siz1= Rcpp::clone(Sou_siz); IntegerVector Indx1= Rcpp::clone(Indx);

  //--------------------
  for( i=1;i<nits;i++){
    Par(i,1)= Rcpp::rgamma(2,N-1+c,1/(sum1+d))[0];
    //update Alpha

    Par(i,0)= Rcpp::rgamma(2,N-1+a,1/(sum+b))[0];

    //Update the infection times
    int k=1;

    while(k<N){

      //Symptomatic individual
      int ll=k;  // choose site seqentially but with larger data can scan choise ie randomly chosen site
      IntegerVector samp=seq_len(Sou_siz[ll]);
      RNGScope scope;
      int j= Rcpp::RcppArmadillo::sample(samp,1,false)[0]-1;
      int ii=Soursamp(rt,Inftim,ll,Sou_siz[ll],j);
      double rnew=mat(std::min(ii,ll),std::max(ii,ll));
      double rold=mat(std::min(Sourc[k],ll),std::max(Sourc[k],ll));
      pacc=rold/rnew*exp(-Par(i,0)*(rnew-rold));

      if(pacc> R::runif(0,1)){
        Sourc[k]=ii;
        sum=sum-rold+rnew;
        if(k==10){
          acc+=1;
        }

      }

      k++;
    }
    //Update to the infectious period
    rtnew=rtim + R::runif(-2,2);
    int j;
    for(j=1;j<N;j++){
      if(Inftim[j]>(Inftim[Sourc[j]]+rtnew)){
        break;
      }
    }
    if(rtnew>0 && rtnew<lmax && j==(N)){
      construct(Coo, Inftim, rtnew,N,Tmax,Indx1,Sou_siz1,dat,nrows1);
      Sou_siz1[0]=1;

      if(is_false(any(Sou_siz1==0))){  // no disjoint chain
        double sum_new=integ(dat,Tmax,nrows1);
        //double product = std::accumulate(Sou_siz.begin(), Sou_siz.end(), 1.0, std::multiplies<double>());
        //double product1 = std::accumulate(Sou_siz1.begin(), Sou_siz1.end(), 1.0, std::multiplies<double>());
        pacc=exp(-Par(i,1)*(sum_new-sum1));
        if(pacc>R::runif(0,1)){
          rtim=rtnew;
          sum1=sum_new;
          nrows=nrows1;
          for(int ll=0;ll<N;ll++){
            Sou_siz[ll]=Sou_siz1[ll];
            Indx[ll]=Indx1[ll];
            rt[ll]=Inftim[ll]+rtim;
          }
        }
      }
    }

    Par(i,2)=rtim;
    for(int j=1;j<N;j++){
      Sour(i,j)=Sourc[j]+1;
    }
  }

  return Par;

}



//' @rdname raymcmc
//' @export
// [[Rcpp::export]]
NumericMatrix caumcmc(NumericMatrix Coo, IntegerVector Sourc,NumericVector Inftim, int nits=1000,int lmin=0,int lmax=120,int Tmax=325,double a=1e-4, double b=10e-6,
                     double c=1e-3, double d=1e-4, NumericVector theta=NumericVector::create(1.0,2.0,110)){
  int N=Inftim.size(); int nrows;
  double lnew,lold;

  // Define a large matrix to store all events
  //Initial data construction
  NumericMatrix dat(2*N,5);
  IntegerVector Sou_siz(N);
  IntegerVector Indx(N);
  NumericVector rt=Inftim + theta(2);
  construct(Coo, Inftim, theta(2),N,Tmax,Indx,Sou_siz,dat,nrows);
  Sou_siz[0]=1;

  double pacc;
  //int acc=1;
  //Define the outputs
  NumericMatrix Par(nits,3);   // Paremters including alpha,beta,lamda(infectious period)
  NumericMatrix Sour(nits,N);  // Source of each infection

  double rtnew;
  //double sum=0;
  double rtim=theta[2];
  //Interdistances between all teh premisses
  NumericMatrix mat(N,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      mat(i,j)=distance(Coo(i,0),Coo(j,0),Coo(i,1),Coo(j,1));    // distance from the source to the neighborgh country
    }

  }

  Par(0,_)=theta;
  lold=lik(Par(0,0),Sourc,mat);

  for(int j=1;j<N;j++){
    Sour(0,j)=Sourc[j]+1;
  }

  double sum1=integ(dat,Tmax,nrows); // initial integration in beta
  int i; int nrows1;
  IntegerVector Sou_siz1= Rcpp::clone(Sou_siz); IntegerVector Indx1= Rcpp::clone(Indx);

  //--------------------
  for( i=1;i<nits;i++){
    Par(i,1)= Rcpp::rgamma(2,N-1+c,1/(sum1+d))[0];
    //update Alpha

    double alphanew=Rcpp::rnorm(1,Par(i-1,0),5)[0];
    if(alphanew>0){
      lnew=lik(alphanew,Sourc,mat);
      pacc=exp(lold-lnew);
      if(pacc>R::runif(0,1)){
        Par(i,0)=alphanew;
        lold=lnew;
      }
      else{
        Par(i,0)=Par(i-1,0);
      }
    }
    else{
      Par(i,0)=Par(i-1,0);
    }


    //Update the infection times
    int k=1;
    while(k<N){
      //Symptomatic individual
      int ll=k;  // choose site seqentially but with larger data can scan choise ie randomly chosen site

      IntegerVector samp=seq_len(Sou_siz[ll]);
      int j= Rcpp::RcppArmadillo::sample(samp,1,false)[0]-1;
      int ii=Soursamp(rt,Inftim,ll,Sou_siz[ll],j);
      double rnew=mat(std::min(ii,ll),std::max(ii,ll));
      double rold=mat(std::min(Sourc[k],ll),std::max(Sourc[k],ll));
      pacc=(1+rold/(Par(i,0)*Par(i,0)))/(1+rnew/(Par(i,0)*Par(i,0)));
      if(pacc> R::runif(0,1)){
        Sourc[k]=ii;
        lold=lik(Par(i,0),Sourc,mat);
      }


      k++;
    }
    //Update to the infectious period

    rtnew=rtim + R::runif(-2,2);
    int j;
    for(j=1;j<N;j++){
      if(Inftim[j]>(Inftim[Sourc[j]]+rtnew)){
        break;
      }
    }
    if(rtnew>0 && rtnew<lmax && j==(N)){
      construct(Coo, Inftim, rtnew,N,Tmax,Indx1,Sou_siz1,dat,nrows1);
      Sou_siz1[0]=1;
     // if(is_false(any(Sou_siz1==0))){  // no disjoint chain
        double sum_new=integ(dat,Tmax,nrows1);
        //double product = std::accumulate(Sou_siz.begin(), Sou_siz.end(), 1.0, std::multiplies<double>());
        //double product1 = std::accumulate(Sou_siz1.begin(), Sou_siz1.end(), 1.0, std::multiplies<double>());
        pacc=exp(-Par(i,1)*(sum_new-sum1));
        if(pacc>R::runif(0,1)){
          rtim=rtnew;
          sum1=sum_new;
          nrows=nrows1;
          for(int ll=0;ll<N;ll++){
            Sou_siz[ll]=Sou_siz1[ll];
            Indx[ll]=Indx1[ll];
            rt[ll]=Inftim[ll]+rtim;
          }
        }
      //}
    }
    Par(i,2)=rtim;
    for(int j=1;j<N;j++){
      Sour(i,j)=Sourc[j]+1;
    }

  }

  return Par;

}

