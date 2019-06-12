library(phaseR)
library(deSolve)


rm(list=ls()) #reinitialiser variables

Eq3 = function(t,x,parameters) #modèle complexe
{
  with(as.list(c(parameters,x)),{
    dS = - beta*S*I -lambda*S  + proies_necessaires*probabilite_rencontre*(S+I+R)*P
    dI =   beta*S*I - gamma*I - mu*I -lambda*I  
    dR =   mu*I - lambda*R  
    dD =   gamma*I + lambda*(R+I+S) 
    dP =   sigma*P*(1-P/kp )  - probabilite_rencontre*(S+I+R)*P #on estime quune rencontre = 1 mort
    list(c(dS,dI,dR,dD,dP))
  })
}





simulation_logistique = function(ini,inisansmaladie,parameters,maxtime,pasdetemps,
                                 ploting = TRUE)
{
  #simulation sans infectés initiaux
  
  time=seq(0,maxtime,pasdetemps)
  solm=lsoda(func=Eq3,y=inisansmaladie,times=time,parms = parameters) #resoud
  populationm="en vie"
  extinctionm=NA
  # On vérifie que la population de proie ne s'éteint pas
  for (j in 1:maxtime/pasdetemps)
  {
    if (solm[j,"P"]<1) #moins d1 proie => population éteinte
    {
      extinctionm=j #indice du temps à laquelle l'espèce s'éteind
      populationm = "eteinte"
      break
    }

  }

  if (populationm =="eteinte")
  {
    
    #On relance une simulation pour continuer l'évolution après l'extinction
    S=as.numeric(solm[extinctionm,"S"])
    I=as.numeric(solm[extinctionm,"I"]) 
    R=as.numeric(solm[extinctionm,"R"])
    D=as.numeric(solm[extinctionm,"D"])
    P=0
    
    continuer=c(S=S,I=I,R=R,D=D,P=P)
    time2=seq(0,maxtime-time[extinctionm],pasdetemps) # temps moins important
    
    solm2=lsoda(func=Eq3,y=continuer,times=time2,parms = parameters) #resoud
    
    morts_sans_maladie=solm2[length(solm2[,"D"]),"D"]
    
    #On combine les résultats pendant et après l'épidémie
    s=c(solm[1:extinctionm-1,"S"],solm2[,"S"]) #on tronque pour ne pas compter 2 fois le jour final de l'épidémie
    i=c(solm[1:extinctionm-1,"I"],solm2[,"I"])
    r=c(solm[1:extinctionm-1,"R"],solm2[,"R"])
    d=c(solm[1:extinctionm-1,"D"],solm2[,"D"])
    p=c(solm[1:extinctionm-1,"P"],solm2[,"P"])
    
  }
  else
  {
    #On récupère les données épidémiologiques à la fin de l'épidémie
    s=solm[,"S"]
    i=solm[,"I"]
    r=solm[,"R"]
    d=solm[,"D"]
    p=solm[,"P"]
    morts_sans_maladie=solm[length(solm[,"D"]),"D"] #5eme colonne est celle des morts
  }
  
  #simulation avec des infectés initiaux
  
  sol=lsoda(func=Eq3,y=ini,times=time,parms = parameters) #resoud
  
  
  #On cherche a regardersi l'épidémie s'éteint : le nombre d'infecté <1
  epidemie="endémique"
  fin_epidemie=NA
  for (j in 1:maxtime/pasdetemps)
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
  
  # On vérifie que la population de prédateurs ne s'éteint pas
  population="en vie"
  extinction=NA
  
  for (j in 1:maxtime/pasdetemps)
  {
    if (pmaladie[j]<1) #moins d1 infecté => maladie éteinte
    {
      extinction=j #indice du temps à laquelle l'espèce s'éteind
      population = "eteinte"
      break
    }
  }
  if (population =="eteinte")
  {
    
    #On relance une simulation pour continuer l'évolution après l'épidémie
    S=as.numeric(smaladie[extinction])
    I=as.numeric(imaladie[extinction]) #round(as.numeric(sol[fin_epidemie,"S"]))
    R=as.numeric(rmaladie[extinction])
    D=as.numeric(dmaladie[extinction])
    P=0
    
    continuer=c(S=S,I=I,R=R,D=D,P=P)
    time2=seq(0,maxtime-time[extinction],pasdetemps) # temps moins important
    
    sol3=lsoda(func=Eq3,y=continuer,times=time2,parms = parameters) #resoud
    
    morts_avec_maladie=sol3[length(sol3[,"D"]),"D"]
    
    #On combine les résultats pendant et après l'épidémie
    smaladie=c(smaladie[1:extinction-1],sol3[,"S"]) #on tronque pour ne pas compter 2 fois le jour final de l'épidémie
    imaladie=c(imaladie[1:extinction-1],sol3[,"I"])
    rmaladie=c(rmaladie[1:extinction-1],sol3[,"R"])
    dmaladie=c(dmaladie[1:extinction-1],sol3[,"D"])
    pmaladie=c(pmaladie[1:extinction-1],sol3[,"P"])
    
  }
  else
  {
    extinction =0
    population ="en vie"
  }
  
  
  
  
  #affichage des courbes
  if (ploting)
  {
    par(mfrow=c(3,1)) # 3 fig
    ylim=c(0,max(smaladie,imaladie,rmaladie))
    plot(x=time,y=smaladie,xlab="Temps",ylab="Population",type="l",ylim=ylim,col="green",main = "Graphique avec maladie")
    points(x = time,imaladie,col="red",type="l")
    points(x = time,rmaladie,col="blue",type="l")
    points(x = time,dmaladie,col="black",type="l")
    points(x = time,smaladie+imaladie+rmaladie,col="purple",type="l")
    legend("topright", legend=c("S","I","R","D","k","N"), col=c("green","red","blue","black","brown","purple"),lty=1)

    plot(x=time,y=s,xlab="Temps",ylab="Population",type="l",ylim=ylim,col="green",main = "Graphique sans maladie")
    points(x = time,i,col="red",type="l")
    points(x = time,r,col="blue",type="l")
    points(x = time,d,col="black",type="l")
    points(x = time,smaladie+imaladie+rmaladie,col="purple",type="l")
    legend("topright", legend=c("S","D","Nmaladie"), col=c("green","black","purple"),lty=1)

    
    plot(x=time,y=pmaladie,xlab="Temps",ylab="Population",type="l",ylim=c(0,max(pmaladie,p)),col="green",main = "Comparaison Population de proies")
    points(x = time,p,col="red",type="l")
    legend("topright", legend=c("Pmaladie","Psansmaladie"), col=c("green","red"),lty=1)
    abline(h=parameters['kp'],col="brown")
  }
  
  
  # résultat intéresant : le nombre d'individus ayant vécus
  nbr_sans_maladie = sum(s[length(smaladie)],i[length(smaladie)],r[length(smaladie)],d[length(smaladie)])  
  nbr_avec_maladie = sum(smaladie[length(smaladie)],imaladie[length(imaladie)],rmaladie[length(rmaladie)],dmaladie[length(dmaladie)])
  
  
  # autre résultat intéresant : la somme des "années" de vies de tous les individus
  années_de_vie_sans_maladie=( sum(s)+sum(i)+sum(r) ) 
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
         nbr_sans_maladie = nbr_sans_maladie,
         population_proie_sans_maladie = populationm,
         population_proie_avec_maladie = population,
         extinction_des_proies_sans_maladie = extinctionm,
         extinction_des_proies_avec_maladie = extinction))
}



#Données espérance positive
beta=0.00225 #infectiosité
gamma=0.04 #mortalité
mu=0.095 #recouvrement
N=300 #pop initiale
I=1 #infecté initiale
lambda=0.001 #morts naturelles

#paramètre des proies
proies_necessaires =4
probabilite_rencontre=1.75*10^-6 #10^-6
M=N*4 #population initiale proie
kp=M*2 #capacité limite proie
sigma= 0.002 #taux d'accroissement naturelle des proies 0.002


parameters=c(beta=beta,gamma=gamma,N=N,mu=mu,lambda=lambda, proies_necessaires =proies_necessaires,probabilite_rencontre=probabilite_rencontre,kp=kp,sigma=sigma)
ini=c(S=N-I,I=I,R=0,D=0,P=M)
inisansmaladie=c(S=N,I=0,R=0,D=0,P=M)
maxtime=4000
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

testksup=simulation_logistique(ini=ini,inisansmaladie=inisansmaladie,parameters=parameters,
                               maxtime=maxtime*2,pasdetemps=pasdetemps)
print(testksup)









