# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
# Fitting Results for the dataset (Firefox 3.5, Version 2)
# The datasets of the testing time, the number of detected errors, and the number of corrected errors

#Parameters' value for the previous version
A1=63.45664
M1=62.44345
R1=49.56172

# The datasets of the testing time, the number of detected errors, and the number of corrected errors (Firefox 3.5, Version 2)
CPU=c(1:28)
detected=c(66,73,76,81,83,87,88,92,94,94,94,99,102,104,105,105,106,106,107,
	108,108,109,112,113,113,115,115,116)
corrected=c(48,51,51,54,55,57,59,63,64,65,66,70,73,73,76,77,80,80,81,83,83,
	86,88,91,92,93,93,93)
data1=data.frame(CPU,detected,corrected)

#Set the functions of g and k
g=function(p,b) (p+b)/2
k=function(p,b,q) sqrt(g(p,b)^2-p*b*(1-q))

#Set hyperbolic functions
HR=function(p,b,q,t) g(p,b)*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)
HR2=function(p,b,q,t) (b*M1+(p-g(p,b))*R1)*sinh(k(p,b,q)*t)+R1*k(p,b,q)*cosh(k(p,b,q)*t)
HM=function(p,b,q,t) (g(p,b)-p*(1-q))*sinh(k(p,b,q)*t)+k(p,b,q)*cosh(k(p,b,q)*t)
HM2=function(p,b,q,t) ((b-g(p,b))*M1+p*q*R1)*sinh(k(p,b,q)*t)+M1*k(p,b,q)*cosh(k(p,b,q)*t)

#Set the functions of M(t), R(t) and A(t)
#set R(t)
R=function(a,p,b,q,t) (A1+a-q*R1)/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HR(p,b,q,t))+exp(-g(p,b)*t)/k(p,b,q)*HR2(p,b,q,t)
#set M(t)=R(t)+(1/b)*dR(t)/dt
M=function(a,p,b,q,t) (A1+a-q*R1)/(1-q)*(1-exp(-g(p,b)*t)/k(p,b,q)*HM(p,b,q,t))+exp(-g(p,b)*t)/k(p,b,q)*HM2(p,b,q,t)
#set A(t)
A=function(a,p,b,q,t) A1+a-q*R1+q*R(a,p,b,q,t)

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
(optPar=optim(par=c(43,0.12,0.06,0.30),fn=MSE,lower=c(0,0,0,0),method="L-BFGS-B")$par)
MSE(optPar)

#Generate the series data for the plot
(A2=A(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))
(M2=M(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))
(R2=R(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))

#Generate the figure of the fitting result for Firefox dataset 3.5
png(
  filename = "firefox2N.png",
  type = "cairo", 
  res = 300, # 300dpi 
  width = 2400, height = 2400,  #size
  bg = "transparent" # background color
)
par(cex=1.5)
plot(CPU,detected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),pch=17,xlab=" ",ylab=" ",main="Firefox 3.5")
curve(M(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="black",lwd=2,add=T)
par(new=T)
par(cex=1.5)
plot(CPU,corrected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",pch=16,xlab="cumulative testing time",ylab="cumulative faults")
curve(R(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",lwd=2,add=T)
par(cex=1.2)
legend(0.56*max(CPU),0.3*max(detected),c("Actual Cumu. Detec.","Actual Cumu. Correc.","Fit. Detec. Func.","Fit. Correc. Func."),pch=c(17,16,NaN,NaN),lwd=c(NaN,NaN,2,2),col=c("1","gray","1","gray"))

dev.off()

