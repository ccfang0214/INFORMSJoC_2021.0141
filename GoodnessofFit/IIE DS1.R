# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
# Fitting Results for the dataset (IIE_2016_DS1 data)

# The datasets of the testing time, the number of detected errors, and the number of corrected errors for IIE_2016_DS1 data
CPU=c(1:35)
detected=c(20,40,56,70,85,95,103,109,120,127,131,139,146,148,156,160,165,
	170,173,176,177,179,183,184,186,188,188,192,193,194,194,194,195,195,196)
corrected=c(8,21,37,58,69,82,90,99,110,119,126,133,140,146,153,158,163,
	166,167,172,174,177,180,183,185,186,187,189,193,193,194,194,194,194,195)
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
(optPar=optim(par=c(57,0.47,0.96,0.72),fn=MSE,lower=c(0,0,0,0),method="L-BFGS-B")$par)
MSE(optPar)

#Generate the figure of the fitting result for IIE DS1
png(
  filename = "DS1N.png",
  type = "cairo", 
  res = 300, # 300dpi 
  width = 2400, height = 2400,  #size
  bg = "transparent" # background color
)

par(cex=1.5)
plot(CPU,detected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),pch=17,xlab=" ",ylab=" ",main="DS1")
curve(M(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="black",lwd=2,add=T)
par(new=T)
par(cex=1.5)
plot(CPU,corrected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",pch=16,xlab="cumulative testing time",ylab="cumulative faults")
curve(R(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",lwd=2,add=T)
par(cex=1.2)
legend(0.56*max(CPU),0.3*max(detected),c("Actual Cumu. Detec.","Actual Cumu. Correc.","Fit. Detec. Func.","Fit. Correc. Func."),pch=c(17,16,NaN,NaN),lwd=c(NaN,NaN,2,2),col=c("1","gray","1","gray"))

dev.off()