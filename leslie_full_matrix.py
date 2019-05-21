# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:20:21 2019

@author: tijod
"""

import numpy as np
import matplotlib.pyplot as plt


#Population initiale pour les 4 classes d'individus
N1=500    # Nombre d'enfants initial
N2=540   # Nombre d'ados initial
N3=458   # Nombre d'adultes initial
N4=354    # Nombre de vieux initial


#Initialisation de la matrice de population

  #colonnes : Enfants | Ados| Adultes| Vieux
pop=np.array([[N1,N2,N3,N4],   # Sains
              [0,0,0,0],       # Infectés
              [0,0,0,0],       # Immunisés
              [0,0,0,0]])      # Morts
 


# Paramètres 
beta=0.00225   # Taux d'infectiosité
gamma=0.3   # Taux de mortalité
mu=0.095 # Taux de recouvrement
lambda0=0.02  # Taux de mortalité naturelles
alpha=0.1   #taux de naissances naturelles
theta=alpha-lambda0   #pression démographique
k=800   #capacité limite
croissance1,croissance2,croissance3=0.2,0.09,0.1    # taux de croissance des enfants, ados et adultes
maxtime=100    #Durée en jours, mois ou années?
pasdetemps=0.1   

# Matrice de Leslie statiques
leslie=np.array([[-croissance1,croissance1,0,0],            # Les coefficients c1,c2,c3 de la matrice représentent les taux de croissance pour les enfants, ados et adultes
                 [0,-croissance2,croissance2,0],
                 [0,0,-croissance3,croissance3],
                 [0,0,0,0]])

Mnatalité=np.array([[0,0,0,0],
                   [0,0,0,0],
                   [alpha,0,0,0],
                   [alpha,0,0,0]])

Mbébésain=np.array([[1,1,1,0],
                   [0,0,0,0],
                   [0,0,0,0],
                   [0,0,0,0]])

sup_death=np.array([[1,0,0,0],
                    [0,1,0,0],
                    [0,0,1,0],
                    [0,0,0,0]])

# Fonctions
#Prend en entrée les population de chaque classe d'age et retourne la population restante après un pas de temps de l'épidémie



def demographie(matrix,POP_tot):
    N=POP_tot[-1]
    S_=sum(matrix[0,:])     # nombre total de personnes saines dans la population
    Mépidémie=np.array([[0,-beta*S_,0,0],
                        [0,beta*S_-gamma-mu,0,0],
                        [0,mu,0,0],
                        [0,gamma,0,0]]) 
    
    Mmortalité=np.array([[-lambda0-theta*N/k,0,0,0],
                       [0,-lambda0-theta*N/k,0,0],
                       [0,0,-lambda0-theta*N/k,0],
                       [lambda0+theta*N/k,lambda0+theta*N/k,lambda0+theta*N/k,0]])
    matrix=matrix+np.dot(Mépidémie,matrix) +np.dot(np.dot(sup_death,matrix),leslie)+np.dot(np.dot(Mbébésain,matrix),Mnatalité)+np.dot(Mmortalité,matrix)       
    return(matrix)


def simulation(pop):
    #initialisation des listes
    evolution_population_avec_maladie=[pop]
    POP=[[sum(pop[:3,0]),sum(pop[:3,1]),sum(pop[:3,2]),sum(pop[:3,3])]]
    POP_tot=[sum(POP[0])]
    Death=[sum(pop[3,:])]
    #Death=[sum()]
    
    
    # Utilisation de la matrice de Leslie, croissance de population
    for i in range(maxtime):
        pop=demographie(pop,POP_tot)
        POP.append([sum(pop[:3,0]),sum(pop[:3,1]),sum(pop[:3,2]),sum(pop[:3,3])])
        POP_tot.append(sum(POP[i+1]))
        Death.append(sum(pop[3,:]))
        evolution_population_avec_maladie.append(pop)
    
    
    #graphique classe d'âge
    
    time=[i for i in range(maxtime+1)]
    
    plt.xlabel('temps')
    plt.ylabel('population')
    
    plt.plot(time,np.array(POP)[:,0],label="enfant")
    plt.plot(time,np.array(POP)[:,1],label="ados")
    plt.plot(time,np.array(POP)[:,2],label="adultes")
    plt.plot(time,np.array(POP)[:,3],label="vieux")
    plt.plot(time,np.array(Death),'b--',label="morts")
    plt.axhline(y=k, color='k' )
    
    plt.plot(time,POP_tot,'r--',label="population totale")
    plt.legend()
    plt.close()
    # graphique SIR
    
    plt.xlabel('temps')
    plt.ylabel('population')
    
    print(np.shape(evolution_population_avec_maladie[:][0][1]))
    print(np.shape(time))
    plt.plot(time,evolution_population_avec_maladie[:][0][1],label="enfant infecte")
    plt.plot(time,np.array(POP)[:,1],label="ados")
    plt.plot(time,np.array(POP)[:,2],label="adultes")
    plt.plot(time,np.array(POP)[:,3],label="vieux")
    plt.plot(time,np.array(Death),'b--',label="morts")
    plt.axhline(y=k, color='k' )
    
    plt.plot(time,POP_tot,'r--',label="population totale")
    plt.legend()
    plt.close()
    
simulation(pop)