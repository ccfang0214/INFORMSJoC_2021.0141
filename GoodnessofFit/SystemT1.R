# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
# Fitting Results for the dataset (SystemT1 data)

# The datasets of the testing time, the number of detected errors, and the number of corrected errors for SystemT1
CPU=c(4,8.3,10.3,10.9,13.2,14.8,16.6,31.3,56.4,60.9,70.4,78.9,108.4,130.4,
  	169.9,195.9,220.9,252.3,282.3,295.1,300.1)
detected=c(2,2,2,3,4,6,7,16,29,31,42,44,55,69,87,99,111,126,132,135,136)
corrected=c(1,2,2,3,4,4,5,7,13,17,18,32,37,56,75,85,97,117,129,131,136)
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
(optPar=optim(par=c(380,0.00154,0.0389,0),fn=MSE,lower=c(0,0,0,0),method="L-BFGS-B")$par)
MSE(optPar)

#Generate the figure of the fitting result for SystemT1
png(
  filename = "SystemT1N.png",
  type = "cairo", 
  res = 300, # 300dpi 
  width = 2400, height = 2400,  #size
  bg = "transparent" # background color
)

par(cex=1.5)
plot(CPU,detected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),pch=17,xlab=" ",ylab=" ",main="System T1")
curve(M(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="black",lwd=2,add=T)
par(new=T)
par(cex=1.5)
plot(CPU,corrected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",pch=16,xlab="cumulative testing time",ylab="cumulative faults")
curve(R(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",lwd=2,add=T)
par(cex=1.2)
legend(0.56*max(CPU),0.3*max(detected),c("Actual Cumu. Detec.","Actual Cumu. Correc.","Fit. Detec. Func.","Fit. Correc. Func."),pch=c(17,16,NaN,NaN),lwd=c(NaN,NaN,2,2),col=c("1","gray","1","gray"))

dev.off()