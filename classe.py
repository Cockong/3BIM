#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:01:21 2019

@author: jlouison
"""


# Facteurs épidémiologiques
inf=0.5
recov=0.2
death=0.3

class enfants:
    def __init__(self,N1,S,I,R):
        self.N1=N1
        self.S=S
        self.I=I
        self.R=R
        
    def disease(self):
        S_=self.S
        self.S=self.S-inf*self.S
        I_=self.I
        self.I=self.I+inf*S_-recov*self.I-death*self.I
        self.R=self.R+recov*I_
        self.N1=self.N1-death*I_  
        

class ados:
    def __init__(self,N2,S,I,R):
        self.N2=N2
        self.S=S
        self.I=I
        self.R=R
        
    def disease(self):
        S_=self.S
        self.S=self.S-inf*self.S
        I_=self.I
        self.I=self.I+inf*S_-recov*self.I-death*self.I
        self.R=self.R+recov*I_
        self.N2=self.N2-death*I_  
        
        
        
class adultes:
    def __init__(self,N3,S,I,R):
        self.N3=N3
        self.S=S
        self.I=I
        self.R=R
        
    def disease(self):
        S_=self.S
        self.S=self.S-inf*self.S
        I_=self.I
        self.I=self.I+inf*S_-recov*self.I-death*self.I
        self.R=self.R+recov*I_
        self.N3=self.N3-death*I_  
        
        
        
class vieux:
    def __init__(self,N4,S,I,R):
        self.N4=N4
        self.S=S
        self.I=I
        self.R=R
        
    def disease(self):
        S_=self.S
        self.S=self.S-inf*self.S
        I_=self.I
        self.I=self.I+inf*S_-recov*self.I-death*self.I
        self.R=self.R+recov*I_
        self.N4=self.N4-death*I_  
        
        
        
    

    
    
        
        
        