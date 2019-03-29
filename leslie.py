#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:55:33 2019

@author: jlouison
"""
import classe as cl
import numpy as np
import matplotlib.pyplot as plt


#Population initiale pour les 4 classes d'individus
N1=50
N2=70
N3=120
N4=56


# initialisation de l'épidémie 
enfant=cl.enfants(N1,N1,0,0)
ado=cl.ados(N2,N2,0,0)
adulte=cl.adultes(N3,N2,0,0)
vieu=cl.vieux(N4,N2,0,0)

# Fonction 
#Prend en entrée les population de chaque classe d'age et retourne la population restante après un pas de temps de l'épidémie
def epidemie(N1,N2,N3,N4):
    enfant.disease()
    ado.disease()
    adulte.disease()
    vieu.disease()
    N1,N2,N3,N4=enfant.N1,ado.N2,adulte.N3,vieu.N4       # Si arrondissement, utiliser round pour N1,N2,N3 et N4
    return(N1,N2,N3,N4)

#facteurs de croissance des populations
f1,f2,f3=0.7,0.2,0.5
d1,d2,d3=0.1,0.3,0.2


leslie=np.array([[0,f1,f2,f3],[d1,0,0,0],[0,d2,0,0],[0,0,0,d3]])
pop=np.array(list(epidemie(N1,N2,N3,N4)))


#initialisation des listes
POP=[list(pop)]
POP_tot=[sum(list(pop))]


n=100   #Nombre de pas

# Utilisation de la matrice de Leslie, croissance de population
for i in range(n):
    pop=pop+leslie.dot(pop)
    pop=list(epidemie(pop[0],pop[1],pop[2],pop[3]))
    POP.append(pop)
    POP_tot.append(sum(pop))


#graphique

time=[i for i in range(n+1)]

plt.xlabel('temps')
plt.ylabel('population')

plt.plot(time,list(np.array(POP)[:,0]))
plt.plot(time,list(np.array(POP)[:,1]))
plt.plot(time,list(np.array(POP)[:,2]))
plt.plot(time,list(np.array(POP)[:,3]))

plt.plot(time,POP_tot,'r--')



