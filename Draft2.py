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
x1=input("Enter the Liquid Phase Molar Fraction of Component 1=")
R=0.08206

'''Defining the Arrays:'''
Pc=[None]*2 #Critical Pressure
Tc=[None]*2 #Critical Temperature
w=[None]*2 #Pitzer Acentric Factor
Tr=[None]*2 #Reduced Temperature
k1=[None]*2
k2=[None]*2
h1=[None]*2
h2=[None]*2
h3=[None]*2
alpha=[None]*2
a=[None]*2 #Attraction Parameter
b=[None]*2 #Co-Volume Prarmeter
c=[None]*2 #Semi Co-Volume Parameter
Vsol=[None]*2 #Molar Volume
Z=[None]*2 #Compressibilty Factor
lnphi=[None]*2 #Fugacity Coefficient
Zc=[None]*2

#DATA:
'''WATER + CARBON TETRACHLORIDE'''
Pc[0]=221.2; Pc[1]=45.6 #bar
Tc[0]=647.3; Tc[1]=556.4 #K
w[0]=0.343; w[1]=0.191
Zc[0]=0.3175-0.0364*w[0]-0.0245*w[0]**2
Zc[1]=0.3175-0.0364*w[1]-0.0245*w[1]**2

#HKM2 PARAMETERS:
m=-0.5
n=-0.5

def parameters(h2):
    return h2**3+(2-n**2+n-3*Zc[i])*h2**2+(3*Zc[i]+(1+n)*(1-3*Zc[i]))*h2-Zc[i]**3
    
for i in range(0,2):
    h2[i]=fsolve(parameters,0) #phib
    h3[i]=(1-3*Zc[i]+(1-n)*h2[i])/m #phic
    h1[i]=3*Zc[i]**2-n*m*h2[i]*h3[i]+(n*h2[i]+m*h3[i])*(h2[i]+1) #hia
   
    Tr[i]=T/Tc[i] #Reduced Temp. 
#DEFINING THE PARAMETERS:
    k1[i]=(0.0821 + 0.3042*w[i]-0.0730*w[i]**2)
    k2[i]=(3.058 + 1.5479*Tr[i])*(1-Tr[i]*scipy.exp(k1[i]))
    alpha[i]=scipy.exp(k2[i])

    a[i]=h1[i]*R**2*Tc[i]**2/Pc[i]*alpha[i]
    b[i]=h2[i]*R*Tc[i]/Pc[i]
    c[i]=h3[i]*R*Tc[i]/Pc[i]

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