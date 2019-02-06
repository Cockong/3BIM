library(phaseR)
library(deSolve)


rm(list=ls()) #reinitialiser variables

par(mfrow=c(2,1)) # 2 fig


Eq1 = function(t,x,parameters)
{
  with(as.list(c(parameters,x)),{
    dS = - beta*S*I  
    dI =   beta*S*I - gamma*I - mu*I 
    dR =   mu*I 
    dD =   gamma*I
    list(c(dS,dI,dR,dD))
  })
}

beta=0.00225 #infectiosité
gamma=0.01 #mortalité
mu=0.04 #recouvrement
N=764 #pop initiale
I=1 #infecté initiale
parameters=c(beta=beta,gamma=gamma,N=N,mu=mu)
ini=c(S=N-I,I=I,D=0,R=0)
time=seq(0,40,0.1)


sol=lsoda(func=Eq1,y=ini,times=time,parms = parameters) #resoud

plot(x=time,y=sol[,2],xlab="x",ylab="Population",type="l",ylim=c(0,700),col="green")
points(x = time,sol[,3],col="red")
points(x = time,sol[,4],col="blue",type="l")
points(x = time,sol[,5],col="black",type="l")




#rajout des accroissements naturels

Eq2 = function(t,x,parameters) #modèle complexe
{
  with(as.list(c(parameters,x)),{
    dS = - beta*S*I  + alpha*(S+I+R) -lambda*S
    dI =   beta*S*I - gamma*I - mu*I -lambda*I
    dR =   mu*I - lambda*R
    dD =   gamma*I + lambda*(R+I+S)
    list(c(dS,dI,dR,dD))
  })
}

beta=0.00225 #infectiosité
gamma=0.01 #mortalité
mu=0.04 #recouvrement
N=764 #pop initiale
I=1 #infecté initiale
lambda=0.01
alpha=0.02
parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha)
ini=c(S=N-I,I=I,D=0,R=0)
time=seq(0,40,0.1)

sol=lsoda(func=Eq2,y=ini,times=time,parms = parameters) #resoud

plot(x=time,y=sol[,2],xlab="x",ylab="Population",type="l",ylim=c(0,700),col="green")
points(x = time,sol[,3],col="red")
points(x = time,sol[,4],col="blue",type="l")
points(x = time,sol[,5],col="black",type="l")
