set.seed(2232)
# Consider the toy dataset of 20 observations

Y=c(-0.39,0.12,0.94,1.67,1.76,2.44,3.72,4.28,4.92,5.53,0.06,0.48,1.01,1.68,1.80,3.25,4.12,4.60,5.28,6.22)


# Make Histrogram
hist(Y,breaks=c(-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5),col='lightblue')

N=length(Y)
y=as.matrix(Y,nrow=N)

# Get Basic summary statistics
# Initalize first guess
# Mixing parameter pi

ybar=mean(Y)

pi=0.5

mu1=sample(Y,1)
mu2=sample(Y,1)

var1=sum(((Y-ybar)^2)/N)
var2=sum(((Y-ybar)^2)/N)

init_guess=c(pi,mu1,mu2,var1,var2)

# Expectation Maximization

makeEM=function(Y,parameters,iter){
  
  #Responsiblity 
  
  responsibility=matrix(0,nrow=N)
  for(j in 1:iter){
  
  for(i in 1:N){
  
  model2=parameters[1]* dnorm(Y[i],mean=parameters[3],sd=sqrt(parameters[5]))
  
  model1=(1-parameters[1])*dnorm(Y[i],mean=parameters[2],sd=sqrt(parameters[4]))
  
  responsibility[i]=model2/(model1+model2)

  }
    
  R=responsibility
  
  mu1=sum((1-R)*Y)/sum(1-R)
  
  mu2=sum(R*Y)/sum(R)
  
  var1=sum((1-R)*(Y-mu1)^2)/sum(1-R)
  
  var2=sum(R*(Y-mu2)^2)/sum(R)
  
  pi=sum(R)/N

  parameters=c(pi,mu1,mu2,var1,var2)

  }
  return(parameters)

}

# Run on 30 iterations

runEM=makeEM(Y,init_guess,iter=30)

print(paste('The final parameters are: pi=',round(runEM[1],3),'mu1=',round(runEM[2],3),'mu2=',round(runEM[3],3),'var1=',round(runEM[4],3),'var2=',round(runEM[5],3)))

pi=runEM[1]
mu1=runEM[2]
mu2=runEM[3]
var1=runEM[4]
var2=runEM[5]

library(distr)

myMix <- UnivarMixingDistribution(Norm(mean=mu1, sd=sqrt(var1)), Norm(mean=mu2, sd=sqrt(var2)), mixCoeff=c(1-pi, pi))
                                
rmyMix <- r(myMix)

# Sample a million random variates
x <- rmyMix(1e6)

# Histogram of mixture
hist(x[x>-0.5 & x<6.5], breaks=50, col="red", main="",xlab='y',freq=FALSE,ylab='Density')



















