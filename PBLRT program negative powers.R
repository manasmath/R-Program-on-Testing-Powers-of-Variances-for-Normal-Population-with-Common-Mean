#Hypothesis Testing about Powers of Scale Parameters of Two Normal Populations With Common Mean#
#Pravash Jena, 5th December 2024#
#Hypothesis testing-Negative Powers#
#mu=0,sigma1=1# LRT+PBLRT+CAT+MCAT#
##########################################################################
library(MASS)
library(nleqslv)
M=10000
B=10000
n1=5
n2=5
alpha=0.05
c=1
mu=0
# sd2=1
# sigma2=sd2^2
lambda1=1
sigma1=(1/lambda1)^(1/c)
sd1=sqrt(sigma1)
lambda10=1
sigma10=(1/lambda10)^(1/c)
sd10=sqrt(sigma10)
# for(lambda1 in c(1.1,1.3,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0)){
#   sigma1=(1/lambda1)^(1/c)
#   sd1=sqrt(sigma1)
for(sd2 in c(0.25,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0)){
  sigma2=sd2^2
  #.....................................
  x1=matrix(0,n1,M)
  x2=matrix(0,n2,M)
  x1bar=array(0,M)
  x2bar=array(0,M)
  s1=array(0,M)
  s2=array(0,M)
  muml=array(0,M)
  sigma1ml=array(0,M)
  sd1ml=array(0,M)
  sigma2ml=array(0,M)
  sd2ml=array(0,M)
  lambda1ml=array(0,M)
  lambda1u=array(0,M)
  dml=array(0,M)
  murml=array(0,M)
  sigma2rml=array(0,M)
  sd2rml=array(0,M)
  L=array(0,M)
  L0=array(0,M)
  lL=array(0,M)
  Lml=array(0,M)
  Uml=array(0,M)
  z3=array(0,M)
  Lpml=array(0,M)
  Upml=array(0,M)
  v1=0
  w1=0
  v2=0
  w2=0
  u=0
  v=0
  #......Calling functions for deteriming MLEs based on original samples using Newton method.........
  f1newton=function(x){ 
    y=numeric(3)
    y[1]=x[1]-((((n1*x1bar[j])/x[2])+((n2*x2bar[j])/x[3]))/((n1/x[2])+(n2/x[3])))
    y[2]=x[2]-(s1[j]/n1)-((((n2*(x1bar[j]-x2bar[j]))/x[3])/((n1/x[2])+(n2/x[3])))^2)
    y[3]=x[3]-(s2[j]/n2)-((((n1*(x2bar[j]-x1bar[j]))/x[2])/((n1/x[2])+(n2/x[3])))^2)
    y
  }
  f2newton=function(x){ 
    y=numeric(2)
    y[1]=x[1]-((((n1*x1bar[j])/sigma10)+((n2*x2bar[j])/x[2]))/((n1/sigma10)+(n2/x[2])))
    y[2]=x[2]-(s2[j]/n2)-((((n1*(x2bar[j]-x1bar[j]))/sigma10)/((n1/sigma10)+(n2/x[2])))^2)
    y
  }
  #......Calling functions for deteriming MLEs based on bootstrap samples using Newton method.........
  f1newtonb=function(x){
    yy=numeric(3)
    yy[1]=x[1]-((((n1*x11bar[i])/x[2])+((n2*x22bar[i])/x[3]))/((n1/x[2])+(n2/x[3])))
    yy[2]=x[2]-(s11[i]/n1)-((((n2*(x11bar[i]-x22bar[i]))/x[3])/((n1/x[2])+(n2/x[3])))^2)
    yy[3]=x[3]-(s22[i]/n2)-((((n1*(x22bar[i]-x11bar[i]))/x[2])/((n1/x[2])+(n2/x[3])))^2)
    yy
  }
  f2newtonb=function(x){ 
    yy=numeric(2)
    yy[1]=x[1]-((((n1*x11bar[i])/sigma10)+((n2*x22bar[i])/x[2]))/((n1/sigma10)+(n2/x[2])))
    yy[2]=x[2]-(s22[i]/n2)-((((n1*(x22bar[i]-x11bar[i]))/sigma10)/((n1/sigma10)+(n2/x[2])))^2)
    yy
  }
  for(j in 1:M){
    x1[,j]=rnorm(n1,mu,sd1)
    x2[,j]=rnorm(n2,mu,sd2)
    x1bar[j]=mean(x1[,j])
    x2bar[j]=mean(x2[,j])
    s1[j]=sum((x1[,j]-x1bar[j])^2)
    s2[j]=sum((x2[,j]-x2bar[j])^2)
    #..................Newton Method for determining MLEs using original samples.........................................
    xstart=c((n2*x1bar[j]+n1*x2bar[j])/(n1+n2),s1[j]/(n1-1),s2[j]/(n2-1))
    muml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[1]
    sigma1ml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[2]
    sigma2ml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[3]
    lambda1ml[j]=sigma1ml[j]^(-c)
    sd1ml[j]=sqrt(sigma1ml[j])
    sd2ml[j]=sqrt(sigma2ml[j])
    dml[j]=(lambda1ml[j]-lambda10)^2
    #..................Newton Method for computing rectricted MLE under the null hypothesis........................
    xstart=c((n2*sigma10+n1*x2bar[j])/(n1+n2),s2[j]/(n2-1))
    murml[j]=nleqslv(xstart,f2newton,control=list(btol=0.001),method="Newton")$x[1]
    sigma2rml[j]=nleqslv(xstart,f2newton,control=list(btol=0.001),method="Newton")$x[2]
    sd2rml[j]=sqrt(sigma2rml[j])
    #...........................Likelihood ratio statistic using original samples...................
    L[j]=prod(dnorm(x1[,j], mean = muml[j], sd = sd1ml[j]))*prod(dnorm(x2[,j], mean = muml[j], sd = sd2ml[j]))
    L0[j]=prod(dnorm(x1[,j], mean = murml[j], sd= sd10))*prod(dnorm(x2[,j], mean = murml[j], sd = sd2rml[j]))
    lL[j]=-2*log(L0[j]/L[j])
    chi=qchisq(0.95,1,FALSE)
    if(lL[j]>chi){v=v+1}
    #............Bootstrap replication in innerloop.........................#
    x11=matrix(0,n1,B)
    x22=matrix(0,n2,B)
    x11bar=array(0,B)
    x22bar=array(0,B)
    s11=array(0,B)
    s22=array(0,B)
    mubml=array(0,B)
    sigma1bml=array(0,B)
    sigma2bml=array(0,B)
    lambda1bml=array(0,B)
    sd1bml=array(0,B)
    sd2bml=array(0,B)
    d0=array(0,B)
    murbml=array(0,B)
    sigma2rbml=array(0,B)
    sd2rbml=array(0,B)
    Lb=array(0,B)
    L0b=array(0,B)
    lLb=array(0,B)
    for(i in 1:B){
      x11[,i]=rnorm(n1,murml[j],sd10)
      x22[,i]=rnorm(n2,murml[j],sd2rml[j])
      x11bar[i]=mean(x11[,i])
      x22bar[i]=mean(x22[,i])
      s11[i]=sum((x11[,i]-x11bar[i])^2)
      s22[i]=sum((x22[,i]-x22bar[i])^2)
      #..................Newton Method for determining MLEs using the bootstrap replicates.........................................
      xstart=c((n2*x11bar[i]+n1*x22bar[i])/(n1+n2),s11[i]/(n1-1),s22[i]/(n2-1))
      mubml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[1]
      sigma1bml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[2]
      sigma2bml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[3]
      lambda1bml[i]=sigma1bml[i]^(-c)
      sd1bml[i]=sqrt(sigma1bml[i])
      sd2bml[i]=sqrt(sigma2bml[i])
      d0[i]=(lambda1bml[i]-lambda10)^2
      #..................Newton Method for determining bootstrap rectricted MLE under the null hypothesis...............
      xstart=c((n2*sigma10+n1*x22bar[i])/(n1+n2),s22[i]/(n2-1))
      murbml[i]=nleqslv(xstart,f2newtonb,control=list(btol=0.001),method="Newton")$x[1]
      sigma2rbml[i]=nleqslv(xstart,f2newtonb,control=list(btol=0.001),method="Newton")$x[2]
      sd2rbml[i]=sqrt(sigma2rbml[i])
      #...........................Likelihood ratio statistic using bootstrap replications...................
      Lb[i]=prod(dnorm(x11[,i], mean = mubml[i], sd = sd1bml[i]))*prod(dnorm(x22[,i], mean = mubml[i], sd = sd2bml[i]))
      L0b[i]=prod(dnorm(x11[,i], mean = murbml[i], sd= sd10))*prod(dnorm(x22[,i], mean = murbml[i], sd = sd2rbml[i]))
      lLb[i]=-2*log(L0b[i]/Lb[i])
    }
    #.......PBLRT..................
    Lpsml=sort(lLb)
    Upml[j]=Lpsml[(1-(alpha/2))*B]
    Lpml[j]=Lpsml[(alpha/2)*B]
    if(lL[j]>Upml[j]){v2=v2+1}
    if(lL[j]<Lpml[j]){w2=w2+1}
    #.......CAT..................
    Lsml=sort(lambda1bml)
    Uml[j]=Lsml[(1-(alpha/2))*B]
    Lml[j]=Lsml[(alpha/2)*B]
    if(lambda1ml[j]>Uml[j]){v1=v1+1}
    if(lambda1ml[j]<Lml[j]){w1=w1+1}
    #.......MCAT..................
    zs3=sort(d0)
    z3[j]=zs3[(1-alpha)*B]
    if(dml[j]>z3[j]){u=u+1}
  }
  rho=sd2/sd1
  catm=(v1+w1)/M
  pblrt=(v2+w2)/M
  mcat=u/M
  lrt=v/M
  cat("\n",rho,"\t",lrt,"\t",pblrt,"\t",catm,"\t",mcat)
}
