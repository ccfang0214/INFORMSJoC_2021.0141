# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
#################### Case III ####################

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

#set the parameters and the cost function for version 1
c11=90
c12=100
c13=160
c14=1530
a1=29
p1=0.4568
b1=0.0507
q1=0.6861


C1=function(t1)
{
  c11*t1+c12*M(a1,p1,b1,q1,t1)+c13*R(a1,p1,b1,q1,t1)+c14*B(a1,p1,b1,q1,t1)
}

print("### Version 1 ###")
#for sequence solution
print("=== for sequence solution: Testing Time, Total Cost, A, M, R ===")
min.t1=optimize(C1,c(0,500),tol=0.0001)
seqsol.t1=min.t1$minimum
seqsol.c1=C1(seqsol.t1)
m1A1=A(a1,p1,b1,q1,seqsol.t1)
m1M1=M(a1,p1,b1,q1,seqsol.t1)
m1R1=R(a1,p1,b1,q1,seqsol.t1)
m1RC1=c14*B(a1,p1,b1,q1,seqsol.t1)  #risk cost
sprintf("A1= %f ; M1= %f ; R1= %f " , m1A1, m1M1, m1R1)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , seqsol.t1, seqsol.c1, m1RC1)
curve(C1(x),xlim=c(0,seqsol.t1*1.5))  #draw the curve of the cost for sequence solution, version 1

#for DP
print("=== for DP solution: Testing Time, Total Cost, A, M, R ")
dp.t1=111.4287
dp.c1=C1(dp.t1)
m3A1=A(a1,p1,b1,q1,dp.t1)
m3M1=M(a1,p1,b1,q1,dp.t1)
m3R1=R(a1,p1,b1,q1,dp.t1)
m3RC1=c14*B(a1,p1,b1,q1,dp.t1)  #risk cost
sprintf("A1= %f ; M1= %f ; R1= %f " , m3A1, m3M1, m3R1)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , dp.t1, dp.c1, m3RC1)

#draw M, R, A for version 1
curve(M(a1,p1,b1,q1,x),xlim=c(0,dp.t1*1.2), col = 2)  #red
curve(R(a1,p1,b1,q1,x),xlim=c(0,dp.t1*1.2),add = TRUE, col = 4)  #blue
curve(A(a1,p1,b1,q1,x),xlim=c(0,dp.t1*1.2),add = TRUE, col = 3)   #green #blue

#################### for version 2 ####################
#set the parameters and the cost function for version 2
c21=85
c22=90
c23=126
c24=180
a2=43
p2=0.123
b2=0.0632
q2=0.301

Cv2=function(VA1,VM1,VR1,tv)
{
  c21*tv+c22*(Mv(VA1,VM1,VR1,a2,p2,b2,q2,tv)-VM1)  +c23*(Rv(VA1,VM1,VR1,a2,p2,b2,q2,tv)-VR1) +c24*Bv(VA1,VM1,VR1,a2,p2,b2,q2,tv)
}
Cv2T1=function(tv)
{ Cv2(m1A1,m1M1,m1R1,tv) }
Cv2T3=function(tv)
{ Cv2(m3A1,m3M1,m3R1,tv) }

print("### Version 2 ###")
#for sequence solution
print("=== for sequence solution: Testing Time, Total Cost, A, M, R ===")
min.tv=optimize(Cv2T1,c(0,500),tol=0.0001)
seqsol.t2=min.tv$minimum
seqsol.c2=Cv2T1(seqsol.t2)
m1A2=Av(m1A1,m1M1,m1R1,a2,p2,b2,q2,seqsol.t2)
m1M2=Mv(m1A1,m1M1,m1R1,a2,p2,b2,q2,seqsol.t2)
m1R2=Rv(m1A1,m1M1,m1R1,a2,p2,b2,q2,seqsol.t2)
m1RC2=c24*Bv(m1A1,m1M1,m1R1,a2,p2,b2,q2,seqsol.t2)    #risk cost
sprintf("A2= %f ; M2= %f ; R2= %f " , m1A2, m1M2, m1R2)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , seqsol.t2, seqsol.c2, m1RC2)
curve(Cv2T1(x),xlim=c(0,seqsol.t2*1.5))  # draw the curve of the cost for sequence solution, version 2


#for DP
print("=== for DP solution: Testing Time, Total Cost, A, M, R ======")
dp.t2=82.2313  # the value got from DP program in advance
dp.c2=Cv2T3(dp.t2)
m3A2=Av(m3A1,m3M1,m3R1,a2,p2,b2,q2,dp.t2)
m3M2=Mv(m3A1,m3M1,m3R1,a2,p2,b2,q2,dp.t2)
m3R2=Rv(m3A1,m3M1,m3R1,a2,p2,b2,q2,dp.t2)
m3RC2=c24*Bv(m3A1,m3M1,m3R1,a2,p2,b2,q2,dp.t2) #risk cost
sprintf("A2= %f ; M2= %f ; R2= %f " , m3A2, m3M2, m3R2)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , dp.t2, dp.c2, m3RC2)

#draw M, R, A for version 2
curve(Mv(m3A1,m3M1,m3R1,a2,p2,b2,q2,x),xlim=c(0,dp.t2*1.2), col = 2)  #red
curve(Rv(m3A1,m3M1,m3R1,a2,p2,b2,q2,x),xlim=c(0,dp.t2*1.2),add = TRUE, col = 4)  #blue
curve(Av(m3A1,m3M1,m3R1,a2,p2,b2,q2,x),xlim=c(0,dp.t2*1.2),add = TRUE, col = 3)   #green #blue
curve(Cv2T3(x),xlim=c(0,dp.t2*1.5))  # draw the curve of the cost for DP solution, version 2

#################### for version 3 ####################
#set the parameters and the cost function for version 3
c31=95
c32=105
c33=170
c34=1180
a3=13
p3=0.1066
b3=0.0341
q3=0.4791

Cv3=function(VA1,VM1,VR1,tv)
{
  c31*tv+c32*(Mv(VA1,VM1,VR1,a3,p3,b3,q3,tv)-VM1)+c33*(Rv(VA1,VM1,VR1,a3,p3,b3,q3,tv)-VR1)+c34*Bv(VA1,VM1,VR1,a3,p3,b3,q3,tv)
}
Cv3T1=function(tv)
{ Cv3(m1A2,m1M2,m1R2,tv) }
Cv3T3=function(tv)
{ Cv3(m3A2,m3M2,m3R2,tv) }

print("### Version 3 ###")
#for sequence solution
print("=== for sequence solution: Testing Time, Total Cost, A, M, R ===")
min.tv=optimize(Cv3T1,c(0,500),tol=0.0001)
seqsol.t3=min.tv$minimum
seqsol.c3=Cv3T1(seqsol.t3)
m1A3=Av(m1A2,m1M2,m1R2,a3,p3,b3,q3,seqsol.t3)
m1M3=Mv(m1A2,m1M2,m1R2,a3,p3,b3,q3,seqsol.t3)
m1R3=Rv(m1A2,m1M2,m1R2,a3,p3,b3,q3,seqsol.t3)
m1RC3=c34*Bv(m1A2,m1M2,m1R2,a3,p3,b3,q3,seqsol.t3)     #risk cost
sprintf("A3= %f ; M3= %f ; R3= %f " , m1A3, m1M3, m1R3)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , seqsol.t3, seqsol.c3, m1RC3)
curve(Cv3T1(x),xlim=c(0,seqsol.t3*1.5)) # draw the curve of the cost for sequence solution, version 3

#for DP
print("=== for dynamic programming: Testing Time, Total Cost, A, M, R ===")
dp.t3=49.1055     # the value got from DP program in advance
dp.c3=Cv3T3(dp.t3)
m3A3=Av(m3A2,m3M2,m3R2,a3,p3,b3,q3,dp.t3)
m3M3=Mv(m3A2,m3M2,m3R2,a3,p3,b3,q3,dp.t3)
m3R3=Rv(m3A2,m3M2,m3R2,a3,p3,b3,q3,dp.t3)
m3RC3=c34*Bv(m3A2,m3M2,m3R2,a3,p3,b3,q3,dp.t3)     #risk cost
sprintf("A3= %f ; M3= %f ; R3= %f " , m1A3, m1M3, m1R3)
sprintf("Optimal testing time = %f ; Optimal testing cost = %f ; Risk Cost = %f" , dp.t3, dp.c3, m3RC3)

#draw M, R, A for version 3
curve(Mv(m3A2,m3M2,m3R2,a3,p3,b3,q3,x),xlim=c(0,dp.t3*1.2), col = 2)  #red
curve(Rv(m3R2,m3M2,m3R2,a3,p3,b3,q3,x),xlim=c(0,dp.t3*1.2),add = TRUE, col = 4)  #blue
curve(Av(m3A2,m3M2,m3R2,a3,p3,b3,q3,x),xlim=c(0,dp.t3*1.2),add = TRUE, col = 3)   #green #blue
curve(Cv3T3(x),xlim=c(0,dp.t3*1.5))  # draw the curve of the cost for DP solution, version 3
