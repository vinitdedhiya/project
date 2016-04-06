# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 09:30:03 2016

@author: Vinit Dedhiya
"""

import scipy
import scipy.optimize 
from scipy.optimize import fsolve
import matplotlib
import pylab as py

'''User Defined Values:'''
T=input("Enter the Temperature of the System(K)=") #K 
P=input("Enter the Pressure of the System(bar)=") #bar
R=0.08206

'''Defining the Arrays:'''
Pc=[None]*2 #Critical Pressure
Tc=[None]*2 #Critical Temperature
w=[None]*2 #Pitzer Acentric Factor
Tr=[None]*2 #Reduced Temperature
k1=[None]*2
k2=[None]*2
alpha=[None]*2
a=[None]*2 #Attraction Parameter
b=[None]*2 #Co-Volume Prarmeter
c=[None]*2 #Semi Co-Volume Parameter
Vsol=[None]*2 #Molar Volume
Z=[None]*2 #Compressibilty Factor
lnphi=[None]*2 #Fugacity Coefficient

#DATA:
'''WATER + CARBON TETRACHLORIDE'''
Pc[0]=221.2; Pc[1]=45.6 #bar
Tc[0]=647.3; Tc[1]=556.4 #K
w[0]=0.343; w[1]=0.191


#HKM2 PARAMETERS:
m=-0.5
n=-0.5


for i in range(0,2):

    Tr[i]=T/Tc[i] #Reduced Temp. 

    #DEFINING THE PARAMETERS:
    k1[i]=(0.0821 + 0.3042*w[i] - 0.0730*w[i]**2)
    k2[i]=(3.058 + 1.5479*Tr[i])*(1-Tr[i]*scipy.exp(k1[i]))
    alpha[i]=scipy.exp(k2[i])

    #Consider:
    h1=.5; h2=0.25; h3=0.125

    a[i]=h1*R**2*Tc[i]**2/Pc[i]*alpha[i]
    b[i]=h2*R*Tc[i]/Pc[i]
    c[i]=h3*R*Tc[i]/Pc[i]

def vol(V):
    return P-(R*T/(V-b[i]))+(a[i]/(V+n*b[i])/(V+m*c[i]))

for i in range(0,2):    
    Vguess=R*T/P #Guess Value = Ideal Gas Volume
    Vsol[i]=fsolve(vol,Vguess)

print "Molar Volume of Water in litres=",Vsol[0] #litres
print "Molar Volume of CCl4 in litres=",Vsol[1] #litres

Z[0]=P*Vsol[0]/(R*T)
Z[1]=P*Vsol[1]/(R*T)
print "Compressibility of Water =",Z[0] 
print "Compressibility of CCl4 =",Z[1] 

for i in range(0,2):
    lnphi[i]=scipy.log(Z[i]-P*b[i]/R/T)*-1+Z[i]-1+(a[i]/R/T/(m*c[i]-n*b[i]))*scipy.log((Vsol[i]+n*b[i])/(Vsol[i]+m*c[i]))

print "Fugacity Coefficient of Water =",scipy.exp(lnphi[0]) 
print "Fugacity Coefficient of CCl4 =",scipy.exp(lnphi[1])