# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 12:50:58 2016

@author: Vinit
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
x2=1-x1
R=0.08206

'''Defining the Arrays:'''
Pc=[None]*2 #Critical Pressure
Tc=[None]*2 #Critical Temperature
w=[None]*2 #Pitzer Acentric Factor
Tr=[None]*2 #Reduced Temperature
k1=[None]*2
k2=[None]*2
h1=[None]*2 #omega a
h2=[None]*2 #omega b
h3=[None]*2 #omega c
alpha=[None]*2
a=[None]*2 #Attraction Parameter
b=[None]*2 #Co-Volume Prarmeter
c=[None]*2 #Semi Co-Volume Parameter
Vgas=[None]*2 #Gas Phase Molar Volume
Vliq=[None]*2 #Liquid Phase Molar Volume
Zg=[None]*2 #Compressibilty Factor for Gas
Zl=[None]*2 #Compressibilty Factor for Liquid
lnphiG=[None]*2 #Fugacity Coefficient for Gas
lnphiL=[None]*2 #Fugacity Coefficient for Liquid
Zc=[None]*2 #Critical Compressibilty Factor

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
    return h2**3+(2-n**2+n-3*Zc[i])*h2**2+(3*Zc[i]**2+(1+n)*(1-3*Zc[i]))*h2-Zc[i]**3
    
for i in range(0,2):
    h2[i]=fsolve(parameters,0) #omegab
    h3[i]=(1-3*Zc[i]+(1-n)*h2[i])/m #omegac
    h1[i]=3*Zc[i]**2-n*m*h2[i]*h3[i]+(n*h2[i]+m*h3[i])*(h2[i]+1) #omegaa
   
    Tr[i]=T/Tc[i] #Reduced Temp. 
#DEFINING THE PARAMETERS:
    k1[i]=(0.0821 + 0.3042*w[i]-0.0730*w[i]**2)
    k2[i]=(3.058 + 1.5479*Tr[i])*(1-pow(Tr[i],k1[i]))
    alpha[i]=scipy.exp(k2[i])

    a[i]=h1[i]*R**2*Tc[i]**2*alpha[i]/Pc[i]
    b[i]=h2[i]*R*Tc[i]/Pc[i]
    c[i]=h3[i]*R*Tc[i]/Pc[i]

#Explicit Cubic Equation for Pressure as a Function of Volume:
def vol(V):
    return P-(R*T/(V-b[i]))+(a[i]/((V+n*b[i])*(V+m*c[i])))

'''GAS PHASE'''
for i in range(0,2):    
    Vguess=R*T/P #Guess Value = Ideal Gas Volume
    Vgas[i]=fsolve(vol,Vguess)

#Gas Phase Molar Volumes
print ""
print "GAS PHASE Molar Volume of Water in litres=",Vgas[0] #litres
print "GAS PHASE Molar Volume of CCl4 in litres=",Vgas[1] #litres

#Compressibilities
Zg[0]=P*Vgas[0]/(R*T)
Zg[1]=P*Vgas[1]/(R*T)
print ""
print "GAS PHASE Compressibility of Water =",Zg[0] 
print "GAS PHASE Compressibility of CCl4 =",Zg[1] 

#Fugacity Coefficients
for i in range(0,2):
    lnphiG[i]=-(scipy.log(Zg[i]-P*b[i]/R/T))+Zg[i]-1+(a[i]/(R*T*(m*c[i]-n*b[i])))*scipy.log((Vgas[i]+n*b[i])/(Vgas[i]+m*c[i]))

print ""
print "GAS PHASE Fugacity Coefficient of Water =",scipy.exp(lnphiG[0]) 
print "GAS PHASE Fugacity Coefficient of CCl4 =",scipy.exp(lnphiG[1])

'''LIQUID PHASE'''
for i in range(0,2):    
    Vliq[i]=fsolve(vol,b[i])

#Liquid Phase Molar Volumes
print ""
print "LIQUID PHASE Molar Volume of Water in litres=",Vliq[0] #litres
print "LIQUID PHASE Molar Volume of CCl4 in litres=",Vliq[1] #litres

#Compressabilities
Zl[0]=P*Vliq[0]/(R*T)
Zl[1]=P*Vliq[1]/(R*T)
print ""
print "LIQUID PHASE Compressibility of Water =",Zl[0] 
print "LIQUID PHASE Compressibility of CCl4 =",Zl[1] 

#Fugacity Coefficients
for i in range(0,2):
    lnphiL[i]=-(scipy.log(Zl[i]-P*b[i]/R/T))+Zl[i]-1+(a[i]/(R*T*(m*c[i]-n*b[i])))*scipy.log((Vliq[i]+n*b[i])/(Vliq[i]+m*c[i]))
'''
print ""
print "LIQUID PHASE Fugacity Coefficient of Water =",scipy.exp(lnphiL[0]) 
print "LIQUID PHASE Fugacity Coefficient of CCl4 =",scipy.exp(lnphiL[1])
'''

'''COMPONENT FUGACITIES IN RESPECTIVE FLUID PHASE'''

'''Gas Phase'''
U00=((n*b[0]-m*c[0])*Vgas[0]+n*m*(c[0]*b[0]-b[0]*c[0]))/((Vgas[0]+n*b[0])*(Vgas[0]+m*c[0]))
lnphiG1=-(scipy.log((Vgas[0]-b[0])/Vgas[0]))+(b[0]/(Vgas[0]-b[0]))+((2*(x1*a[0]+x2*a[0]))-a*(m*c[0]-n*b[0]))*(scipy.log((Vgas[0]+n*b[0])/(Vgas[0]+m*c[0])))/(R*T*(m*c[0]-n*b[0])**2)+((a[0]/(R*T*(m*c[0]-n*b[0])))*U00)-scipy.log(Zg[0])
print""
print "phiG1=",scipy.exp(lnphiG1[0])
print "phiG2=",scipy.exp(lnphiG1[1])

'''Liquid Phase'''
V00=((n*b[0]-m*c[0])*Vliq[0]+n*m*(c[0]*b[0]-b[0]*c[0]))/((Vliq[0]+n*b[0])*(Vliq[0]+m*c[0]))
lnphiL1=-(scipy.log((Vliq[0]-b[0])/Vliq[0]))+(b[0]/(Vliq[0]-b[0]))+((2*(x1*a[0]+x2*a[0]))-a*(m*c[0]-n*b[0]))*(scipy.log((Vliq[0]+n*b[0])/(Vliq[0]+m*c[0])))/(R*T*(m*c[0]-n*b[0])**2)+((a[0]/(R*T*(m*c[0]-n*b[0])))*V00)-scipy.log(Zl[0])
'''
print""
print "phiL1=",scipy.exp(lnphiL1[0])
print "phiL2=",scipy.exp(lnphiL1[1])
'''

'''STABILITY BASED ON THE PHI-PHI MODEL'''
'''
y1=x1*scipy.exp(lnphiG1[0])/scipy.exp(lnphiL1[0]) #Composition of component in Vapour phase
print ""
print "Vapour Phase Composition of Component 1=",y1
'''

"""References:
    1. A new three-parameter cubic equation of state for calculation of 
        physical properties and vapour-liquid equilibria.
        
        - A. Haghtalab, M. J. Kamali, S. H. Mazloumi, P. Mahmoodi
        
    2. A simple and uniﬁed algorithm to solve ﬂuid phase
       equilibria using either the gamma-phi or the phi-phi
       approach for binary and ternary mixtures

        - Romain PRIVAT∗,a, Jean-Noël JAUBERTa and Yannick PRIVATb
