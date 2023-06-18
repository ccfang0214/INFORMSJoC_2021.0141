################ for version 1 #################
gc()
#set g and k
g=function(p,b) (p+b)/2
k=function(p,b,q) sqrt(g(p,b)^2-p*b*(1-q))
#set hyperbolic functions
HM=function(p,b,q,t) (g(p,b)-p*(1-q))*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)
HR=function(p,b,q,t) g(p,b)*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)
HM2=function(VA1,VM1,VR1,p,b,q,t) ((b-g(p,b))*VM1+p*q*VR1)*sinh(k(p,b,q)*t)+VM1*k(p,b,q)*cosh(k(p,b,q)*t)
HR2=function(VA1,VM1,VR1,p,b,q,t) (b*VM1+(p-g(p,b))*VR1)*sinh(k(p,b,q)*t)+VR1*k(p,b,q)*cosh(k(p,b,q)*t)

#set M(t),R(t),A(t),B(t)
M=function(a,p,b,q,t) a/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HM(p,b,q,t))
R=function(a,p,b,q,t) a/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HR(p,b,q,t))
A=function(a,p,b,q,t) a+q*R(a,p,b,q,t)
B=function(a,p,b,q,t) A(a,p,b,q,t)-R(a,p,b,q,t)

Mv=function(VA1,VM1,VR1,a,p,b,q,t) (VA1+a-q*VR1)/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HM(p,b,q,t))+exp(-g(p,b)*t)/k(p,b,q)*HM2(VA1,VM1,VR1,p,b,q,t)
Rv=function(VA1,VM1,VR1,a,p,b,q,t) (VA1+a-q*VR1)/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HR(p,b,q,t))+exp(-g(p,b)*t)/k(p,b,q)*HR2(VA1,VM1,VR1,p,b,q,t)
Av=function(VA1,VM1,VR1,a,p,b,q,t) VA1+a-q*VR1+q*Rv(VA1,VM1,VR1,a,p,b,q,t)
Bv=function(VA1,VM1,VR1,a,p,b,q,t) Av(VA1,VM1,VR1,a,p,b,q,t)-Rv(VA1,VM1,VR1,a,p,b,q,t)

c11=90
c12=100
c13=160
c14=1530
a1=29
p1=0.4568
b1=0.0507
q1=0.6861

c21=85
c22=90
c23=126
c24=180
a2=43
p2=0.123
b2=0.0632
q2=0.301

c31=95
c32=105
c33=170
c34=1180
a3=13
p3=0.1066
b3=0.0341
q3=0.4791

sCv1=50000
sCv2=38000
sCv3=20000

Cv1=function(t1)
{
  ifelse(t1>0,c11*t1+c12*M(a1,p1,b1,q1,t1)+c13*R(a1,p1,b1,q1,t1)+c14*B(a1,p1,b1,q1,t1),sCv1)
}
Cv2=function(VA1,VM1,VR1,tv)
{
  ifelse(tv>0,c21*tv+c22*(Mv(VA1,VM1,VR1,a2,p2,b2,q2,tv)-VM1)+c23*(Rv(VA1,VM1,VR1,a2,p2,b2,q2,tv)-VR1) +c24*Bv(VA1,VM1,VR1,a2,p2,b2,q2,tv),sCv2)
}
Cv3=function(VA1,VM1,VR1,tv)
{
  ifelse(tv>0,c31*tv+c32*(Mv(VA1,VM1,VR1,a3,p3,b3,q3,tv)-VM1)+c33*(Rv(VA1,VM1,VR1,a3,p3,b3,q3,tv)-VR1)+c34*Bv(VA1,VM1,VR1,a3,p3,b3,q3,tv),sCv3)
}

CostV3=function(pars)
{
  t1=pars[1]
  t2=pars[2]
  t3=pars[3]  
  
  Cv1(t1)+Cv2(A(a1,p1,b1,q1,t1),M(a1,p1,b1,q1,t1),R(a1,p1,b1,q1,t1),t2)+Cv3(Av(A(a1,p1,b1,q1,t1),M(a1,p1,b1,q1,t1),R(a1,p1,b1,q1,t1),a2,p2,b2,q2,t2),Mv(A(a1,p1,b1,q1,t1),M(a1,p1,b1,q1,t1),R(a1,p1,b1,q1,t1),a2,p2,b2,q2,t2),Rv(A(a1,p1,b1,q1,t1),M(a1,p1,b1,q1,t1),R(a1,p1,b1,q1,t1),a2,p2,b2,q2,t2),t3)
}

optim(par=c(100,90,50),CostV3,lower=c(0,0,0),method="L-BFGS-B")

