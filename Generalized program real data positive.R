#Hypothesis Testing about Powers of Scale Parameters of Two Normal Populations With Common Mean#
#Pravash Jena, 15th January 2025#
#Generalized Variable Approach#Real Data Program#Positive Powers
##########################################################################
library(MASS)
M=1
N=10000
n1=9
n2=13
# n1=12
# n2=14
alpha=0.05
c=1
lambda10=209.4191
#lambda10=14.06939
A=((2^c)*gamma(((n1-1)/2)+c))/gamma((n1-1)/2)
B=gamma((n2-1)/2)/(sqrt(2)*gamma(n2/2))
C=gamma((n1-1)/2)/(sqrt(2)*gamma(n1/2))
D=1.1531/2
E=1.0349/2


bu1=array(0,M)
bu2=array(0,M)
bu11=array(0,M)
bu22=array(0,M)
du=array(0,M)
cu=array(0,M)

bgm1=array(0,M)
bgm2=array(0,M)
bgm11=array(0,M)
bgm22=array(0,M)
dgm=array(0,M)
cgm=array(0,M)

bgd1=array(0,M)
bgd2=array(0,M)
bgd11=array(0,M)
bgd22=array(0,M)
dgd=array(0,M)
cgd=array(0,M)

bks1=array(0,M)
bks2=array(0,M)
bks11=array(0,M)
bks22=array(0,M)
dks=array(0,M)
cks=array(0,M)

bmk1=array(0,M)
bmk2=array(0,M)
bmk11=array(0,M)
bmk22=array(0,M)
dmk=array(0,M)
cmk=array(0,M)

btk1=array(0,M)
btk2=array(0,M)
btk11=array(0,M)
btk22=array(0,M)
dtk=array(0,M)
ctk=array(0,M)

bbc11=array(0,M)
bbc12=array(0,M)
bbc111=array(0,M)
bbc122=array(0,M)
dbc1=array(0,M)
cbc1=array(0,M)

bbc21=array(0,M)
bbc22=array(0,M)
bbc211=array(0,M)
bbc222=array(0,M)
dbc2=array(0,M)
cbc2=array(0,M)

for(j in 1:M){
  x1bar=77.33333
  x2bar=80.38462
  s1=1986
  s2=135.0769
  
  # x1bar=109.75
  # x2bar=109.5
  # s1=11*20.74802
  # s2=13*2.729104
  
  gm=(n1*x1bar+n2*x2bar)/(n1+n2)
  gd=((n1*(n1-1)*s2*x1bar)+(n2*(n2-1)*s1*x2bar))/((n1*(n1-1)*s2)+(n2*(n2-1)*s1))
  ks=((n1*(n1-3)*s2*x1bar)+(n2*(n2-3)*s1*x2bar))/((n1*(n1-3)*s2)+(n2*(n2-3)*s1))
  mk=((sqrt(n1*(n1-1)*s2)*x1bar)+(sqrt(n2*(n2-1)*s1)*x2bar))/((sqrt(n1*(n1-1)*s2))+(sqrt(n2*(n2-1)*s1)))
  tk=(((sqrt(n1*s2))*B*x1bar)+((sqrt(n2*s1))*C*x2bar))/((sqrt(n1*s2)*B)+(sqrt(n2*s1)*C))
  bc1=x1bar+(x2bar-x1bar)*((D*s1)/(n1*(n1-1)))/((s1/(n1*(n1-1)))+(s2/(n2*(n2+2)))+(((x2bar-x1bar)^2)/(n2+2)))
  bc2=x1bar+(((x2bar-x1bar)*E*(n2)*(n2-1)*s1)/((n2*(n2-1)*s1)+(n1*(n1-1)*s2)))
  
  #...........unbiased estimator.......
  lambda1u=(s1^c)/A
  
  #............grand mean..............
  lambda1gm=((s1+(n1*(x1bar-gm)^2))^c)/A
  
  #............Grybill and Deal.........
  lambda1gd=((s1+(n1*(x1bar-gd)^2))^c)/A
  
  #............Khatri and Shah..........
  lambda1ks=((s1+(n1*(x1bar-ks)^2))^c)/A
  
  #............Moore and Kishnamoorthy............
  lambda1mk=((s1+(n1*(x1bar-mk)^2))^c)/A
  
  #............Tripathy and Kumar............
  lambda1tk=((s1+(n1*(x1bar-tk)^2))^c)/A
  
  #............BC1............
  lambda1bc1=((s1+(n1*(x1bar-bc1)^2))^c)/A
  
  #............BC2............
  lambda1bc2=((s1+(n1*(x1bar-bc2)^2))^c)/A
  
  
  Tu=array(0,N)
  Tgm=array(0,N)
  Tgd=array(0,N)
  Tks=array(0,N)
  Tmk=array(0,N)
  Ttk=array(0,N)
  Tbc1=array(0,N)
  Tbc2=array(0,N)
  for(i in 1:N){
    u=rchisq(1,n1-1,0)
    v=rchisq(1,n2-1,0)
    w=rchisq(1,1,0)
    
    #...............unbiased.................
    Tu[i]=(A*lambda1u)/(u)^c
    
    #...............Grandmean.................
    Tgm[i]=(((s1/u)^c)*lambda1gm)/((s1+(n1*(n2/(n1+n2))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
    #...............Grybill and Deal.................
    Tgd[i]=(((s1/u)^c)*lambda1gd)/((s1+(n1*((n2*(n2-1)*s1)/((n1*(n1-1)*s2)+(n2*(n2-1)*s1)))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
    #...............Khatri and Shah.................
    Tks[i]=(((s1/u)^c)*lambda1ks)/((s1+(n1*((n2*(n2-3)*s1)/((n1*(n1-3)*s2)+(n2*(n2-3)*s1)))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
    #...............Moore and Kishnamoorthy.................
    Tmk[i]=(((s1/u)^c)*lambda1mk)/((s1+(n1*((sqrt(n2*(n2-1)*s1))/((sqrt(n1*(n1-1)*s2))+(sqrt(n2*(n2-1)*s1))))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
    #...............Tripathy and Kumar.................
    Ttk[i]=(((s1/u)^c)*lambda1tk)/((s1+(n1*(((sqrt(n2*s1))*C)/(((sqrt(n1*s2))*B)+((sqrt(n2*s1))*C)))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
    #...............BC1.................
    Tbc1[i]=(((s1/u)^c)*lambda1bc1)/(((s1)+n1*w*((s1/(n1*u))+(s2/(n2*v)))*(((D*s1)/(n1*(n1-1)))/((s1/(n1*(n1-1)))+(s2/(n2*(n2+2)))+((w*((s1/(n1*u))+(s2/(n2*v))))/(n2+2))))^2)^c/A)
    
    #...............BC2.................
    Tbc2[i]=(((s1/u)^c)*lambda1bc2)/((s1+(n1*((n2*(n2-1)*E*s1)/((n1*(n1-1)*s2)+(n2*(n2-1)*s1)))^2*w*((s1/(n1*u))+(s2/(n2*v)))))^c/A)
    
  }
  #...............unbiased.................
  bu1[j]=length(which(Tu<lambda10))
  bu2[j]=length(which(Tu>lambda10))
  bu11[j]=bu1[j]/N
  bu22[j]=bu2[j]/N
  du[j]=2*min(bu11[j],bu22[j])
  #...............Grandmean.................
  bgm1[j]=length(which(Tgm<lambda10))
  bgm2[j]=length(which(Tgm>lambda10))
  bgm11[j]=bgm1[j]/N
  bgm22[j]=bgm2[j]/N
  dgm[j]=2*min(bgm11[j],bgm22[j])
  #...............Grybill and Deal.................
  bgd1[j]=length(which(Tgd<lambda10))
  bgd2[j]=length(which(Tgd>lambda10))
  bgd11[j]=bgd1[j]/N
  bgd22[j]=bgd2[j]/N
  dgd[j]=2*min(bgd11[j],bgd22[j])
  #...............Khatri and Shah.................
  bks1[j]=length(which(Tks<lambda10))
  bks2[j]=length(which(Tks>lambda10))
  bks11[j]=bks1[j]/N
  bks22[j]=bks2[j]/N
  dks[j]=2*min(bks11[j],bks22[j])
  #...............Moore and Kishnamoorthy.................
  bmk1[j]=length(which(Tmk<lambda10))
  bmk2[j]=length(which(Tmk>lambda10))
  bmk11[j]=bmk1[j]/N
  bmk22[j]=bmk2[j]/N
  dmk[j]=2*min(bmk11[j],bmk22[j])
  #...............Tripathy and Kumar.................
  btk1[j]=length(which(Ttk<lambda10))
  btk2[j]=length(which(Ttk>lambda10))
  btk11[j]=btk1[j]/N
  btk22[j]=btk2[j]/N
  dtk[j]=2*min(btk11[j],btk22[j])
  #...............BC1.................
  bbc11[j]=length(which(Tbc1<lambda10))
  bbc12[j]=length(which(Tbc1>lambda10))
  bbc111[j]=bbc11[j]/N
  bbc122[j]=bbc12[j]/N
  dbc1[j]=2*min(bbc111[j],bbc122[j])
  #...............BC2.................
  bbc21[j]=length(which(Tbc2<lambda10))
  bbc22[j]=length(which(Tbc2>lambda10))
  bbc211[j]=bbc21[j]/N
  bbc222[j]=bbc22[j]/N
  dbc2[j]=2*min(bbc211[j],bbc222[j])
}
cat("\n",lambda10,"\t",du,"\t",dgd,"\t",dks,"\t",dmk,"\t",dtk,"\t",dbc1,"\t",dbc2,"\t",dgm)

