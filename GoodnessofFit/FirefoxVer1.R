# This is the code file for paper "A Study on Optimal Release Schedule for Multi-Version Software"
# Fitting Results for the dataset (Firefox 3.0, Version 1)

# The datasets of the testing time, the number of detected errors, and the number of corrected errors (Firefox 3.0, Version 1)
CPU=c(1:53) 
detected=c(9,12,16,25,27,29,29,32,34,35,36,36,39,39,40,40,40,41,42,43,43,44,
	45,45,46,47,47,49,50,50,50,50,51,52,53,54,55,55,55,55,56,59,60,60,
      60,61,62,62,62,62,62,64,65)
corrected=c(3,3,4,7,9,12,12,13,15,17,18,19,21,22,22,22,23,24,26,26,26,26,26,
      26,26,27,27,28,28,29,33,33,34,34,35,37,37,37,37,37,37,37,38,40,42,
	42,43,45,45,46,46,47,48)
data1=data.frame(CPU,detected,corrected)

#Set the functions of g and k
g=function(p,b) (p+b)/2
k=function(p,b,q) sqrt(g(p,b)^2-p*b*(1-q))

#Set hyperbolic functions
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
(optPar=optim(par=c(29,0.45,0.05,0.68),fn=MSE,lower=c(0,0,0,0),method="L-BFGS-B")$par)
MSE(optPar)

#Generate the series data for the plot
(A1=A(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))
(M1=M(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))
(R1=R(optPar[1],optPar[2],optPar[3],optPar[4],max(CPU)))

#Generate the figure of the fitting result for Firefox dataset 3.0
png(
  filename = "firefox1N.png",
  type = "cairo", 
  res = 300, # 300dpi 
  width = 2400, height = 2400,  #size
  bg = "transparent" # background color
)

par(cex=1.5)
plot(CPU,detected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),pch=17,xlab=" ",ylab=" ",main="Firefox 3.0")
curve(M(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="black",lwd=2,add=T)
par(new=T)
par(cex=1.5)
plot(CPU,corrected,xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",pch=16,xlab="cumulative testing time",ylab="cumulative faults")
curve(R(optPar[1],optPar[2],optPar[3],optPar[4],x),xlim=c(0,max(CPU)),ylim=c(0,max(detected)),col="gray",lwd=2,add=T)
par(cex=1.2)
legend(0.56*max(CPU),0.3*max(detected),c("Actual Cumu. Detec.","Actual Cumu. Correc.","Fit. Detec. Func.","Fit. Correc. Func."),pch=c(17,16,NaN,NaN),lwd=c(NaN,NaN,2,2),col=c("1","gray","1","gray"))

dev.off()

