#Hypothesis Testing about Powers of Scale Parameters of Two Normal Populations With Common Mean#
#Pravash Jena, 15th January 2025#
#Computing p-values, Real Data program#
#LRT+PBLRT+CAT+MCAT#Positive Powers
##########################################################################
library(MASS)
library(nleqslv)
M=1
B=10000
n1=9
n2=13
# n1=12
# n2=14
alpha=0.05
c=1
lambda10=0.0032
#lambda10=0.04244007
sigma10=(1/lambda10)^(1/c)
sd10=sqrt(sigma10)
#....................................
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
#..............P-value..........
pvml1=array(0,M)
pvml2=array(0,M)
pvml11=array(0,M)
pvml22=array(0,M)
pvml=array(0,M)
pvmml=array(0,M)
pvmml1=array(0,M)
pvpb1=array(0,M)
pvpb2=array(0,M)
pvpb11=array(0,M)
pvpb22=array(0,M)
pvpb=array(0,M)
#................Newton method......Original Samples........#
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
#................Newton method......Bootstrap Samples........#
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
  x1[,j]=c(105, 83, 76, 75, 51, 76, 93, 75, 62)
  x2[,j]=c(84, 86, 85, 82, 77, 76, 77, 80, 83, 81, 78, 78, 78)
  x1bar[j]=mean(x1[,j])
  x2bar[j]=mean(x2[,j])
  s1[j]=sum((x1[,j]-x1bar[j])^2)
  s2[j]=sum((x2[,j]-x2bar[j])^2)
  
  # x1bar=109.75
  # x2bar=109.5
  # s1=11*20.74802
  # s2=13*2.729104
  #......................Newton Method (Computing MLEs using original samples)...............
  xstart=c((n2*x1bar[j]+n1*x2bar[j])/(n1+n2),s1[j]/(n1-1),s2[j]/(n2-1))
  muml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[1]
  sigma1ml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[2]
  sigma2ml[j]=nleqslv(xstart,f1newton,control=list(btol=0.001),method="Newton")$x[3]
  lambda1ml[j]=sigma1ml[j]^(-c)
  sd1ml[j]=sqrt(sigma1ml[j])
  sd2ml[j]=sqrt(sigma2ml[j])
  dml[j]=(lambda1ml[j]-lambda10)^2
  #................Resrtricted MLEs under Null hypothesis using original samples..................
  xstart=c((n2*sigma10+n1*x2bar[j])/(n1+n2),s2[j]/(n2-1))
  murml[j]=nleqslv(xstart,f2newton,control=list(btol=0.001),method="Newton")$x[1]
  sigma2rml[j]=nleqslv(xstart,f2newton,control=list(btol=0.001),method="Newton")$x[2]
  sd2rml[j]=sqrt(sigma2rml[j])
  
  L[j]=prod(dnorm(x1[,j], mean = muml[j], sd = sd1ml[j]))*prod(dnorm(x2[,j], mean = muml[j], sd = sd2ml[j]))
  L0[j]=prod(dnorm(x1[,j], mean = murml[j], sd= sd10))*prod(dnorm(x2[,j], mean = murml[j], sd = sd2rml[j]))
  lL[j]=-2*log(L0[j]/L[j])
  
  #................................For summary data...........................................
  # L[j]=((1/(2*(3.141)))^((n1+n2)/2))*(sigma1ml[j])^(-n1/2)*(sigma2ml[j])^(-n2/2)*exp((-1/(2*sigma1ml[j]))*(s1[j]+n1*(x1bar[j]-muml[j])^2))*exp((-1/(2*sigma2ml[j]))*(s2[j]+n2*(x2bar[j]-muml[j])^2))
  # L0[j]=((1/(2*(3.141)))^((n1+n2)/2))*(sigma10)^(-n1/2)*(sigma2rml[j])^(-n2/2)*exp((-1/(2*sigma10))*(s1[j]+n1*(x1bar[j]-murml[j])^2))*exp((-1/(2*sigma2rml[j]))*(s2[j]+n2*(x2bar[j]-murml[j])^2))
  # lL[j]=-2*log(L0[j]/L[j])
  
  chi=qchisq(0.95,1,FALSE)
  if(lL[j]>chi){v=v+1}
  #............Bootstrap loop.........................#
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
    #......................Newton Method (Computing MLEs using bootstrap samples)...............
    xstart=c((n2*x11bar[i]+n1*x22bar[i])/(n1+n2),s11[i]/(n1-1),s22[i]/(n2-1))
    mubml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[1]
    sigma1bml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[2]
    sigma2bml[i]=nleqslv(xstart,f1newtonb,control=list(btol=0.001),method="Newton")$x[3]
    lambda1bml[i]=sigma1bml[i]^(-c)
    sd1bml[i]=sqrt(sigma1bml[i])
    sd2bml[i]=sqrt(sigma2bml[i])
    d0[i]=(lambda1bml[i]-lambda10)^2
    #................Resrtricted MLEs under Null hypothesis using bootstrap samples..................
    xstart=c((n2*sigma10+n1*x22bar[i])/(n1+n2),s22[i]/(n2-1))
    murbml[i]=nleqslv(xstart,f2newtonb,control=list(btol=0.001),method="Newton")$x[1]
    sigma2rbml[i]=nleqslv(xstart,f2newtonb,control=list(btol=0.001),method="Newton")$x[2]
    sd2rbml[i]=sqrt(sigma2rbml[i])
    
    Lb[i]=prod(dnorm(x11[,i], mean = mubml[i], sd = sd1bml[i]))*prod(dnorm(x22[,i], mean = mubml[i], sd = sd2bml[i]))
    L0b[i]=prod(dnorm(x11[,i], mean = murbml[i], sd= sd10))*prod(dnorm(x22[,i], mean = murbml[i], sd = sd2rbml[i]))
    lLb[i]=-2*log(L0b[i]/Lb[i])
  }
  Lpsml=sort(lLb)
  pvpb1[j]=length(which(Lpsml<lL[j]))
  pvpb2[j]=length(which(Lpsml>lL[j]))
  pvpb11[j]=pvpb1[j]/B
  pvpb22[j]=pvpb2[j]/B
  pvpb[j]=2*min(pvpb11[j],pvpb22[j])
  
  Lsml=sort(lambda1bml)
  pvml1[j]=length(which(Lsml<lambda1ml[j]))
  pvml2[j]=length(which(Lsml>lambda1ml[j]))
  pvml11[j]=pvml1[j]/B
  pvml22[j]=pvml2[j]/B
  pvml[j]=2*min(pvml11[j],pvml22[j])
  
  zs3=sort(d0)
  pvmml1[j]=length(which(zs3>dml[j]))
  pvmml[j]=pvmml1[j]/B
}
cat("\n",lambda10,"\t",1-pchisq(lL,1), "\t",pvpb,"\t",pvml,"\t",pvmml)



