library(phaseR)
library(deSolve)


rm(list=ls()) #reinitialiser variables


Eq3 = function(t,x,parameters) #modèle complexe
{
  with(as.list(c(parameters,x)),{
    dS = - beta*S*I  + alpha*(S+I+R) -lambda*S - theta*(S+I+R)*S/k
    dI =   beta*S*I - gamma*I - mu*I -lambda*I - theta*(S+I+R)*I/k
    dR =   mu*I - lambda*R - theta*(S+I+R)*R/k
    dD =   gamma*I + lambda*(R+I+S) + theta*(S+I+R)^2/k
    list(c(dS,dI,dR,dD))
  })
}

beta=0.00225 #infectiosité
gamma=0.25 #mortalité
mu=0.04 #recouvrement
N=764 #pop initiale
I=1 #infecté initiale
lambda=0.001 #morts naturelles
alpha=0.002 #naissances naturelles
theta=alpha-lambda #pression démographique
k=800 #capacité limite
parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha,theta=theta,k=k)
ini=c(S=N-I,I=I,R=0,D=0)
inisansmaladie=c(S=N,I=0,R=0,D=0)
maxtime=5000
time=seq(0,maxtime,1)


solm=lsoda(func=Eq3,y=inisansmaladie,times=time,parms = parameters) #resoud
morts_sans_maladie=solm[length(solm[,5]),5]


sol=lsoda(func=Eq3,y=ini,times=time,parms = parameters) #resoud

epidemie="endémique"
for (j in 1:maxtime)
{
  if (sol[j,3]<1) #moins d1 infecté
  {
    morts_pendant_maladie=sol[j,5]
    fin_epidemie=j
    epidemie="eteinte"
    break
  }
}

S=as.numeric(sol[fin_epidemie,2])
I=0
R=as.numeric(sol[fin_epidemie,4])
D=as.numeric(sol[fin_epidemie,5])

continuer=c(S=S,I=I,R=R,D=D)
time2=seq(0,maxtime-time[fin_epidemie]-1,1) # ne pas prendre 0.1 prendre une longueur max

sol2=lsoda(func=Eq3,y=continuer,times=time2,parms = parameters) #resoud

morts_avec_maladie=sol2[length(sol2[,5]),5]


smaladie=c(sol[1:j,2],sol2[,2])
imaladie=c(sol[1:j,3],sol2[,3])
rmaladie=c(sol[1:j,4],sol2[,4])
dmaladie=c(sol[1:j,5],sol2[,5])



par(mfrow=c(2,1)) # 2 fig

plot(x=time,y=smaladie,xlab="x",ylab="Population",type="l",ylim=c(0,900),col="green",main = "Graphique avec maladie")
points(x = time,imaladie,col="red",type="l")
points(x = time,rmaladie,col="blue",type="l")
points(x = time,dmaladie,col="black",type="l")


plot(x=time,y=solm[,2],xlab="x",ylab="Population",type="l",ylim=c(0,900),col="green",main = "Graphique sans maladie")
points(x = time,solm[,3],col="red",type="l")
points(x = time,solm[,4],col="blue",type="l")
points(x = time,solm[,5],col="black",type="l")
