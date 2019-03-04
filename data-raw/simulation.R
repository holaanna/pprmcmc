
data(ppr)
attach(ppr)
diff=as.Date(Date.of.start.of.the.outbreak)-as.Date(Date.of.start.of.the.outbreak)[1]
inf=as.numeric(diff)
inf[2:26]=inf[2:26]+runif(25)
inf=sort(inf)
dat=matrix(c(ppr$Longitude,ppr$Latitude),nrow=26,ncol=2,byrow=F)
re=110
rem_tim=re
Dat=array(0,c(26,5))
ni=1
t=inf[1]
Tmax=inf[26]
Dat=c(inf[1],dat[1,],ni,0)
inf_tim=inf[2:26]
k=2
i=2
size=0
indx=0
inf_lis=1
S=0
rt=re
while(t<Tmax){

  Min=min(inf_tim[1],min(rem_tim))
  t=t+Min
  if(Min==inf_tim[1]){ # Infection
    ni=ni+1
    Dat=rbind(Dat,c(t,dat[i,],ni,0))
    i=i+1
    j=sample(inf_lis,1)
    S=c(S,j)
    rem_tim=rem_tim-Min
    inf_tim=inf_tim[-1] - Min
    indx=c(indx,k-1)
    size=c(size,ni-1)
    rem_tim=c(rem_tim,re)
    inf_lis=c(inf_lis,k)
    rt=c(rt,re)
  }
  else{  # removal
    ni=ni-1
    j=which(rem_tim==Min)
    Dat=rbind(Dat,c(t,0,0,ni,inf_lis[j]))
    rem_tim=rem_tim[-j]-Min
    inf_lis=inf_lis[-j]
    inf_tim=inf_tim - Min
  }
  k=k+1
}

rt=inf+110
Sourc=numeric(26)
for(i in 2:26){
  a=which(inf[i]>inf)
  b=which(inf[i]<rt)
  l=a[a%in%b]
  Sourc[i]=sample(l,1)
}
S1=S-1
write.table(Dat,file="Coo.txt",sep="\t",col.names=F,row.names=F)
write.table(t(S1),file="Sour.txt",sep=" ",col.names=F,row.names=F)
write.table(t(size),file="Sour_siz.txt",sep=" ",col.names=F,row.names=F)
write.table(t(inf),file="inf_tim.txt",sep=" ",col.names=F,row.names=F)
write.table(t(indx),file="indx.txt",sep=" ",col.names=F,row.names=F)
write.table(t(rt),file="rt.txt",sep=" ",col.names=F,row.names=F)

