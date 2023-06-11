# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
# Fitting Results for the dataset (OpenProj data)

# The datasets of the testing time, the number of detected errors, and the number of corrected errors for OpenProj
CPU=c(1:49)
detected=c(9,15,19,24,28,29,32,36,36,40,41,41,45,47,49,52,52,55,57,58,61,
	61,64,66,67,69,70,73,74,75,75,78,79,80,82,83,83,84,84,85,85,87,87,87,
	89,89,91,91,94)
corrected=c(0,2,11,12,12,13,14,15,20,21,22,22,24,25,30,30,31,31,31,31,31,
	31,33,33,33,33,33,33,33,36,38,38,40,41,41,41,44,44,44,44,45,45,46,46,
	47,62,62,65,65)
data1=data.frame(CPU,detected,corrected)

#Set the functions of g and k
g=function(p,b) (p+b)/2
k=function(p,b,q) sqrt(g(p,b)^2-p*b*(1-q))

#set hyperbolic functions
HR=function(p,b,q,t) g(p,b)*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)
HM=function(p,b,q,t) (g(p,b)-p*(1-q))*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)

#Set the functions of M(t), R(t) and A(t)
R=function(a,p,b,q,t) a/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HR(p,b,q,t))
M=function(a,p,b,q,t) a/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HM(p,b,q,t))
A=function(a,p,b,q,t) a+q*R(a,p,b,q,t)

#Set the function for calculating the MSE value
MSE=function(pars)
{
	a=pars[1]
	p=pars[2]
	b=pars[3]
	q=pars[4]
	mi=detected
	ri=corrected
	ti=CPU
	n=length(CPU)

	sum((M(a,p,b,q,ti)-mi)^2+(R(a,p,b,q,ti)-ri)^2)/(2*n)
}

#Minimize the MSE value and Obtain the values of the model parameters
(optPar=optim(par=c(40,0.27,0.04,0.92),fn=MSE,lower=c(0,0,0,0),method="L-BFGS-B")$par)
MSE(optPar)

#Generate the figure of the fitting result for OpenProj
png(
  filename = "OpenProjN.png",
  type = "cairo", 
  res = 300, # 300dpi 
  width = 2400, height = 2400,  #size
  bg = "transparent" # background color
)

par(cex=1.5)
plot(CPU,detected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),pch=17,xlab=" ",ylab=" ",main="OpenProj")
curve(M(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="black",lwd=2,add=T)
par(new=T)
par(cex=1.5)
plot(CPU,corrected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",pch=16,xlab="cumulative testing time",ylab="cumulative faults")
curve(R(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",lwd=2,add=T)
par(cex=1.2)
legend(0.56*max(CPU),0.3*max(detected),c("Actual Cumu. Detec.","Actual Cumu. Correc.","Fit. Detec. Func.","Fit. Correc. Func."),pch=c(17,16,NaN,NaN),lwd=c(NaN,NaN,2,2),col=c("1","gray","1","gray"))

dev.off()
