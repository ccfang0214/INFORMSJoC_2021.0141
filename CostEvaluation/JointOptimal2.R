
M1=function(a,p,b,q,t)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))
	a/(1-q)*( 1- exp(-g*t)*((g-p*(1-q))/k*sinh(k*t)+cosh(k*t)) )
}

R1=function(a,p,b,q,t)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))
	a/(1-q)*( 1- exp(-g*t)*(g/k*sinh(k*t)+cosh(k*t)) )
}

A1=function(a,p,b,q,t) a+q*R1(a,p,b,q,t)
B1=function(a,p,b,q,t) A1(a,p,b,q,t)-R1(a,p,b,q,t)

M2=function(a,p,b,q,t,aa,pp,bb,qq,tt)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))

	A1=A1(aa,pp,bb,qq,tt)
	M1=M1(aa,pp,bb,qq,tt)
	R1=R1(aa,pp,bb,qq,tt)
	alpha=A1-q*R1+a

	alpha/(1-q)*( 1- exp(-g*t)*((g-p*(1-q))/k*sinh(k*t)+cosh(k*t)) ) + exp(-g*t)*( ((b-g)*M1+p*q*R1)/k*sinh(k*t)+M1*cosh(k*t) )
}

R2=function(a,p,b,q,t,aa,pp,bb,qq,tt)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))

	A1=A1(aa,pp,bb,qq,tt)
	M1=M1(aa,pp,bb,qq,tt)
	R1=R1(aa,pp,bb,qq,tt)
	alpha=A1-q*R1+a

	alpha/(1-q)*( 1- exp(-g*t)*(g/k*sinh(k*t)+cosh(k*t)) ) + exp(-g*t)*( (b*M1+(p-g)*R1)/k*sinh(k*t)+R1*cosh(k*t) )
}

A2=function(a,p,b,q,t,aa,pp,bb,qq,tt)
{
	A1=A1(aa,pp,bb,qq,tt)
	R1=R1(aa,pp,bb,qq,tt)
	alpha=A1-q*R1+a

	alpha+q*R2(a,p,b,q,t,aa,pp,bb,qq,tt)
}

B2=function(a,p,b,q,t,aa,pp,bb,qq,tt)
{
	A2(a,p,b,q,t,aa,pp,bb,qq,tt)-R2(a,p,b,q,t,aa,pp,bb,qq,tt)
}

M3=function(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))

	A2=A2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	M2=M2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	R2=R2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	alpha=A2-q*R2+a

	alpha/(1-q)*( 1- exp(-g*t)*((g-p*(1-q))/k*sinh(k*t)+cosh(k*t)) ) + exp(-g*t)*( ((b-g)*M2+p*q*R2)/k*sinh(k*t)+M2*cosh(k*t) )
}

R3=function(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
{
	g=(p+b)/2
	k=sqrt(g^2-p*b*(1-q))
	
	A2=A2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	M2=M2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	R2=R2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	alpha=A2-q*R2+a

	alpha/(1-q)*( 1- exp(-g*t)*(g/k*sinh(k*t)+cosh(k*t)) ) + exp(-g*t)*( (b*M2+(p-g)*R2)/k*sinh(k*t)+R2*cosh(k*t) )
}

A3=function(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
{
	A2=A2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	M2=M2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	R2=R2(aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
	alpha=A2-q*R2+a

	alpha+q*R3(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
}

B3=function(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
{
	A3(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)-R3(a,p,b,q,t,aa,pp,bb,qq,tt,aaa,ppp,bbb,qqq,ttt)
}

CostV3=function(pars)
{
	t1=pars[1]
	t2=pars[2]
	t3=pars[3]

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
	c23=125
	c24=420
	a2=43
	p2=0.123
	b2=0.0632
	q2=0.301

	c31=95
	c32=105
	c33=170
	c34=915
	a3=13
	p3=0.1066
	b3=0.0341
	q3=0.4791

	c11*t1 + c12*M1(a1,p1,b1,q1,t1) + c13*R1(a1,p1,b1,q1,t1) + c14*B1(a1,p1,b1,q1,t1) + c21*t2 + c22*(M2(a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)-M1(a1,p1,b1,q1,t1)) + c23*(R2(a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)-R1(a1,p1,b1,q1,t1)) + c24*B2(a2,p2,b2,q2,t2,a1,p1,b1,q1,t1) + c31*t3 + c32*(M3(a3,p3,b3,q3,t3,a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)-M2(a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)) + c33*(R3(a3,p3,b3,q3,t3,a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)-R2(a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)) + c34*B3(a3,p3,b3,q3,t3,a2,p2,b2,q2,t2,a1,p1,b1,q1,t1)

}

optim(par=c(100,90,50),CostV3,lower=c(0,0,0),method="L-BFGS-B")

