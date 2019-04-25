# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 15:59:04 2019

@author: tijod
"""


import numpy as np
import matplotlib.pyplot as plt


#Population initiale pour les 4 classes d'individus
N1=50
N2=70
N3=120
N4=56


#Initialisation de la matrice de population
pop=np.array([[N1,N2,N3,N4],[0,0,0,0],[0,0,0,0]])


# Paramètres 
beta=0.00225 #infectiosité
gamma=0.3 #mortalité
mu=0.095 #recouvrement
lambda0=0.02 #morts naturelles
alpha=0.1 #naissances naturelles
theta=alpha-lambda0 #pression démographique
k=800 #capacité limite
c1,c2,c3=0.2,0.09,0.1
maxtime=100   
pasdetemps=0.1

# Fonctions
#Prend en entrée les population de chaque classe d'age et retourne la population restante après un pas de temps de l'épidémie
def epidemie(matrix):
    S_=sum(matrix[0,:])
    I_=sum(matrix[1,:])
    R_=sum(matrix[2,:])
    leslie=np.array([[-beta*I_+alpha-lambda0-theta*S_/k,alpha-theta*S_/k,alpha-theta],[beta*I_-theta*I_/k,-gamma-mu-lambda0-theta*I_/k,-theta*I_/k],[-theta*R_/k,mu-theta*R_/k,lambda0-theta*R_/k]]) 
    matrix=matrix+np.dot(leslie,matrix)
    return(matrix)

# Cette fonction prend en compte la le passage des individus d'une classe d'âge à une autre.
def demographie(matrix):
    leslie2=np.array([[-c1,0,0,0],[c1,-c2,0,0],[0,c2,-c3,0],[0,0,c3,0]])
    matrix=np.transpose(matrix)+np.dot(leslie2,np.transpose(matrix))
    return(np.transpose(matrix))


#initialisation des listes
POP=[[sum(pop[:,0]),sum(pop[:,1]),sum(pop[:,2]),sum(pop[:,3])]]
POP_tot=[sum(POP[0])]


# Utilisation de la matrice de Leslie, croissance de population
for i in range(maxtime):
    pop=epidemie(pop)
    pop=demographie(pop)
    POP.append([sum(pop[:,0]),sum(pop[:,1]),sum(pop[:,2]),sum(pop[:,3])])
    POP_tot.append(sum(POP[i+1]))


#graphique

time=[i for i in range(maxtime+1)]

plt.xlabel('temps')
plt.ylabel('population')

plt.plot(time,np.array(POP)[:,0],label="enfant")
plt.plot(time,np.array(POP)[:,1],label="ados")
plt.plot(time,np.array(POP)[:,2],label="adultes")
plt.plot(time,np.array(POP)[:,3],label="vieux")

plt.plot(time,POP_tot,'r--',label="population totale")
plt.legend()