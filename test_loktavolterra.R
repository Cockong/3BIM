library(phaseR)
library(deSolve)


rm(list=ls()) #reinitialiser variables

Eq3 = function(t,x,parameters) #modèle complexe
{
  with(as.list(c(parameters,x)),{
    dS = - beta*S*I  + alpha*(S+I+R) -lambda*S - theta*(S+I+R)*S/k + proies_necessaires*probabilite_rencontre*S*P
    dI =   beta*S*I - gamma*I - mu*I -lambda*I - theta*(S+I+R)*I/k + proies_necessaires*probabilite_rencontre*I*P
    dR =   mu*I - lambda*R - theta*(S+I+R)*R/k + proies_necessaires*probabilite_rencontre*R*P
    dD =   gamma*I + lambda*(R+I+S) + theta*(S+I+R)^2/k
    dP =   sigma*P*(1-P/kp )  - probabilite_rencontre*(S+I+R)*P #on estime quune rencontre = 1 mort
    list(c(dS,dI,dR,dD,dP))
  })
}





simulation_logistique = function(ini,inisansmaladie,parameters,maxtime,pasdetemps,
                                 ploting = TRUE, ylim=10*c(0,max(parameters['N'],parameters['k'])))
{
  #simulation sans infectés initiaux
  
  time=seq(0,maxtime,pasdetemps)
  solm=lsoda(func=Eq3,y=inisansmaladie,times=time,parms = parameters) #resoud
  morts_sans_maladie=solm[length(solm[,"D"]),"D"] #5eme colonne est celle des morts
  
  
  
  #simulation avec des infectés initiaux
  
  sol=lsoda(func=Eq3,y=ini,times=time,parms = parameters) #resoud
  
  
  #On cherche a regardersi l'épidémie s'éteint : le nombre d'infecté <1
  epidemie="endémique"
  fin_epidemie=NA
  for (j in 1:maxtime)
  {
    if (sol[j,"I"]<1) #moins d1 infecté => maladie éteinte
    {
      morts_pendant_maladie=sol[j,"D"]
      fin_epidemie=j #indice du temps à laquelle l'épidémie s'éteint
      epidemie="eteinte"
      break
    }
  }
  
  
  
  if (epidemie !="endémique")
  {
    
    
    #On relance une simulation pour continuer l'évolution après l'épidémie
    S=as.numeric(sol[fin_epidemie,"S"])
    I=0                               #round(as.numeric(sol[fin_epidemie,"S"]))
    R=as.numeric(sol[fin_epidemie,"R"])
    D=as.numeric(sol[fin_epidemie,"D"])
    P=as.numeric(sol[fin_epidemie,"P"])
    
    continuer=c(S=S,I=I,R=R,D=D,P=P)
    time2=seq(0,maxtime-time[fin_epidemie],pasdetemps) # temps moins important
    
    sol2=lsoda(func=Eq3,y=continuer,times=time2,parms = parameters) #resoud
    
    morts_avec_maladie=sol2[length(sol2[,"D"]),"D"]
    
    #On combine les résultats pendant et après l'épidémie
    smaladie=c(sol[1:fin_epidemie-1,"S"],sol2[,"S"]) #on tronque pour ne pas compter 2 fois le jour final de l'épidémie
    imaladie=c(sol[1:fin_epidemie-1,"I"],sol2[,"I"])
    rmaladie=c(sol[1:fin_epidemie-1,"R"],sol2[,"R"])
    dmaladie=c(sol[1:fin_epidemie-1,"D"],sol2[,"D"])
    pmaladie=c(sol[1:fin_epidemie-1,"P"],sol2[,"P"])
    
  }
  else
  {
    #On récupère les données épidémiologiques à la fin de l'épidémie
    smaladie=sol[,"S"]
    imaladie=sol[,"I"]
    rmaladie=sol[,"R"]
    dmaladie=sol[,"D"]
    pmaladie=sol[,"P"]
    morts_avec_maladie=sol[length(sol[,"D"]),"D"]
    morts_pendant_maladie=morts_avec_maladie
  }
  
  #affichage des courbes
  #par(mfrow=c(2,1)) # 2 fig
  if (ploting)
  {
    plot(x=time,y=smaladie,xlab="Temps",ylab="Population",type="l",ylim=ylim,col="green",main = "Graphique avec maladie")
    points(x = time,imaladie,col="red",type="l")
    points(x = time,rmaladie,col="blue",type="l")
    points(x = time,dmaladie,col="black",type="l")
    points(x = time,smaladie+imaladie+rmaladie,col="purple",type="l")
    legend("topright", legend=c("S","I","R","D","k","N"), col=c("green","red","blue","black","brown","purple"),lty=1)
    abline(h=parameters['k'],col="brown")
    
    plot(x=time,y=solm[,"S"],xlab="Temps",ylab="Population",type="l",ylim=ylim,col="green",main = "Graphique sans maladie")
    points(x = time,solm[,"I"],col="red",type="l")
    points(x = time,solm[,"R"],col="blue",type="l")
    points(x = time,solm[,"D"],col="black",type="l")
    points(x = time,smaladie+imaladie+rmaladie,col="purple",type="l")
    legend("topright", legend=c("S","I","R","D","k","Nmaladie"), col=c("green","red","blue","black","brown","purple"),lty=1)
    abline(h=parameters['k'],col="brown")

    
        
    plot(x=time,y=pmaladie,xlab="Temps",ylab="Population",type="l",ylim=c(0,3000),col="green",main = "Comparaison Population de proies")
    points(x = time,solm[,"P"],col="red",type="l")
    legend("topright", legend=c("Pmaladie","Psansmaladie"), col=c("green","red"),lty=1)
    abline(h=parameters['kp'],col="brown")
  }
  
  
  # résultat intéresant : le nombre d'individus ayant vécus
  nbr_sans_maladie = sum(solm[length(sol[,"D"]),2:5])  
  nbr_avec_maladie = sum(smaladie[length(smaladie)],imaladie[length(imaladie)],rmaladie[length(rmaladie)],dmaladie[length(dmaladie)])
  
  
  # autre résultat intéresant : la somme des "années" de vies de tous les individus
  années_de_vie_sans_maladie=( sum(solm[,"S"])+sum(solm[,"I"])+sum(solm[,"R"]) ) 
  années_de_vie_avec_maladie=( sum(smaladie)+sum(imaladie)+sum(rmaladie) ) 
  
  #simulation de l'espérance de vie : temps de vie total / nombre dindividus ayant vécu
  esperance_sans_maladie = années_de_vie_sans_maladie / nbr_sans_maladie
  esperance_avec_maladie = années_de_vie_avec_maladie / nbr_avec_maladie
  
  # renvoi des résultats
  list(c(epidemie = epidemie, 
         fin_epidemie = fin_epidemie, 
         morts_epargnees = morts_sans_maladie - morts_avec_maladie, 
         temps_de_vie_gagnees = années_de_vie_avec_maladie - années_de_vie_sans_maladie,
         individus_gagnes = nbr_avec_maladie - nbr_sans_maladie, 
         esperance_gagnée = esperance_avec_maladie - esperance_sans_maladie,
         esperance_avec_maladie = esperance_avec_maladie,
         esperance_sans_maladie = esperance_sans_maladie,
         nbr_avec_maladie = nbr_avec_maladie,
         nbr_sans_maladie = nbr_sans_maladie))
}



#Données espérance positive
beta=0.00225 #infectiosité
gamma=0.04 #mortalité
mu=0.095 #recouvrement
N=300 #pop initiale
I=1 #infecté initiale
lambda=0.001 #morts naturelles
alpha=0.002 #naissances naturelles
theta=alpha-lambda #pression démographique
k=800 #capacité limite

#paramètre des proies
proies_necessaires =4
probabilite_rencontre=10^-6
M=N*4 #population initiale proie
kp=N*12 #capacité limite proie
sigma= 0.002 #taux d'accroissement naturelle des proies


parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha,theta=theta,k=k,
             proies_necessaires =proies_necessaires,probabilite_rencontre=probabilite_rencontre,kp=kp,sigma=sigma)
ini=c(S=N-I,I=I,R=0,D=0,P=M)
inisansmaladie=c(S=N,I=0,R=0,D=0,P=M)
maxtime=2000
pasdetemps=0.1




#1 On montre que le modèle est correcte : 

#en initialisant sansinfecté dans les deux cas on obtient bien les mème résultats
#De plus on observe bien els résultats attends pour une population suivant le modèle logistique
par(mfrow=c(2,1))

testblanc=simulation_logistique(ini=inisansmaladie ,inisansmaladie=inisansmaladie,
                                parameters=parameters,maxtime=maxtime,pasdetemps=pasdetemps)
print(testblanc)



#2) On test pour différentes valeurs de K :

#D'après nos résultats on a systématiquement un nombre de morts 
#plus faible dans le cas d'une pandémie que quand il n'y en a pas

# Il faut noter que si l'épidémie est endémique, 
# la capacité limite effective de la population est réduite

#Comme attendu la somme des années de vie est en faveur de la situation sans maladie
#Dans ce cas la on a en effet une population qui est plus longtemps importante

#Pour les mêmes raisons Nous avons beaucoup d'individus qui n'ont jamais vu le jour dans le cas d'une maladie

#voici les résultats sur l'espérance de vie telle que nous l'avons simuler (on suppose que chaque individu vit aussi longtemps) :
# generalement espérance de vie avec maladie < espérance de vie sans maladie 
#V ATTENTION CE RESULTAT VARIE EN FONCTION DES TESTS notemment quand mortalité maladie<< taux de recouvrement et que maladie non endémique
#D'après résultat taux de natalité trop élevé 
par(mfrow=c(2,2))

parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha,theta=theta,k=N*2) #capacité limite > Pop départ
testksup=simulation_logistique(ini=ini,inisansmaladie=inisansmaladie,parameters=parameters,
                               maxtime=maxtime,pasdetemps=pasdetemps)
print(testksup)


parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha,theta=theta,k=N) #capacité limite = Pop départ
testkegal=simulation_logistique(ini=ini,inisansmaladie=inisansmaladie,parameters=parameters,
                                maxtime=maxtime,pasdetemps=pasdetemps)
print(testkegal)


parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda,alpha=alpha,theta=theta,k=N/2) #capacité limite < que  Pop départ
testkinf=simulation_logistique(ini=ini,inisansmaladie=inisansmaladie,parameters=parameters,
                               maxtime=maxtime,pasdetemps=pasdetemps)
print(testkinf)







