# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:58:06 2016

@author: Vinit
"""

import scipy
import scipy.optimize 
from scipy.optimize import fsolve

'''User Defined Values:'''
T=input("Enter the Temperature of the System(K)=") #K 
P=input("Enter the Pressure of the System(bar)=") #bar
#Initial Guess Values
x1=0.01 #Liquid Phase
y1=0.99 #Vapour Phase

R=0.08206 #Universal Gas Constant (lit.atm/K.mol)

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

'''Solving for Omega b'''
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

    a[i]=h1[i]*R**2*Tc[i]**2*alpha[i]/Pc[i] #Attraction Parameter
    b[i]=h2[i]*R*Tc[i]/Pc[i] #Co-Volume Parameter
    c[i]=h3[i]*R*Tc[i]/Pc[i] #Semi Co-Volume Parameter
    
#Explicit Cubic Equation for Pressure as a Function of T,V and Parameters:
def vol(V):
    return P-(R*T/(V-b[i]))+(a[i]/((V+n*b[i])*(V+m*c[i])))

'''VAPOUR PHASE'''
for i in range(0,2):    
    Vguess=R*T/P #Guess Value = Ideal Gas Volume
    Vgas[i]=fsolve(vol,Vguess)

#Vapour Phase Molar Volumes
print ""
print "VAPOUR PHASE Molar Volume of Water in litres at given P and T=",Vgas[0] #litres
print "VAPOUR PHASE Molar Volume of CCl4 in litres at given P and T=",Vgas[1] #litres

'''LIQUID PHASE'''
for i in range(0,2):    
    Vliq[i]=fsolve(vol,0.1) #Guess Values is so choosen to get Lowest Postive Root as Vliq
    
#Liquid Phase Molar Volumes
print ""
print "LIQUID PHASE Molar Volume of Water in litres at given P and T=",Vliq[0] #litres
print "LIQUID PHASE Molar Volume of CCl4 in litres at given P and T=",Vliq[1] #litres

'''COMPRESSIBILITIES'''
#Vapour Phase Compressibilities
Zg[0]=P*Vgas[0]/(R*T)
Zg[1]=P*Vgas[1]/(R*T)

print ""
print "VAPOUR PHASE Compressibility of Pure Water =",Zg[0] 
print "VAPOUR PHASE Compressibility of Pure CCl4 =",Zg[1] 

#Liquid Phase Compressabilities
Zl[0]=P*Vliq[0]/(R*T)
Zl[1]=P*Vliq[1]/(R*T)

print ""
print "LIQUID PHASE Compressibility of Pure Water =",Zl[0] 
print "LIQUID PHASE Compressibility of Pure CCl4 =",Zl[1] 

'''PURE COMPONENT FUGACITIES'''
#Fugacity Coefficients for Pure Fluid Components
for i in range(0,2):
    lnphiG[i]=-(scipy.log(Zg[i]-P*b[i]/R/T))+Zg[i]-1+(a[i]/(R*T*(m*c[i]-n*b[i])))*scipy.log((Vgas[i]+n*b[i])/(Vgas[i]+m*c[i]))
    lnphiL[i]=-(scipy.log(Zl[i]-P*b[i]/R/T))+Zl[i]-1+(a[i]/(R*T*(m*c[i]-n*b[i])))*scipy.log((Vliq[i]+n*b[i])/(Vliq[i]+m*c[i]))


'''ALGORITHM FOR ARRIVING AT 'X1' AND 'Y1' AT TH GIVEN VALUES OF PRESS. AND TEMP.'''

"""Using the phi-phi Model for Stability Analysis"""

counter=0
error=1

'''For Mixture Properties:'''
a1=a[0]; a2=a[1]; b1=b[0]; b2=b[1]; c1=c[0]; c2=c[1];

Vg1=Vgas[0]; Vg2=Vgas[1]; Vl1=Vliq[0]; Vl2=Vliq[1];

Zg1=Zg[0]; Zg2=Zg[1]; Zl1=Zl[0]; Zl2=Zl[1]

#Liquid Phase Mixture Properties
bml=x1*b1+(1-x1)*b2
cml=x1*c1+(1-x1)*c2
a11=a1; a22=a2
a12=((a1*a2)**0.5)*0.5; a21=a12 #Considering k12=0.5
aml=(x1**2*a11+(1-x1)**2*a2+2*x1*(1-x1)*a12)

#Vapour Phase Mixture Properties
bmg=y1*b1+(1-y1)*b2
cmg=y1*c1+(1-y1)*c2
amg=(y1**2*a11+(1-y1)**2*a2+2*y1*(1-y1)*a12)


while error>10**-10:
    
    '''COMPONENT FUGACITIES IN RESPECTIVE FLUID PHASE'''

    '''COMPONENT 1 - WATER'''    
    '''Vapour Phase'''
    U00=((n*b1-m*c1)*Vg1+n*m*(cmg*b1-bmg*c1))/((Vg1+n*bmg)*(Vg1+m*cmg))
    lnphiG1=-(scipy.log((Vg1-bmg)/Vg1))+(b1/(Vg1-bmg))+((2*y1*a1)-amg*(m*c1-n*b1))*(scipy.log((Vg1+n*bmg)/(Vg1+m*cmg)))/(R*T*(m*cmg-n*bmg)**2)+((amg/(R*T*(m*cmg-n*bmg)))*U00)-scipy.log(Zg1)

    '''Liquid Phase'''
    V00=((n*b1-m*c1)*Vl1+n*m*(cml*b1-bml*c1))/((Vl1+n*bml)*(Vl1+m*cml))
    lnphiL1=-(scipy.log((Vl1-bml)/Vl1))+(b1/(Vl1-bml))+((2*x1*a1)-aml*(m*c1-n*b1))*(scipy.log((Vl1+n*bml)/(Vl1+m*cml)))/(R*T*(m*cml-n*bml)**2)+((aml/(R*T*(m*cml-n*bml)))*V00)-scipy.log(Zl1)

    phiG1=scipy.exp(lnphiG1);    phiL1=scipy.exp(lnphiL1)
       
    '''COMPONENT 2 - CCl4'''
    '''Vapour Phase'''
    U10=((n*b2-m*c2)*Vg2+n*m*(cmg*b2-bmg*c2))/((Vg2+n*bmg)*(Vg2+m*cmg))
    lnphiG2=-(scipy.log((Vg2-bmg)/Vg2))+(b2/(Vg2-bmg))+((2*(1-y1)*a2)-amg*(m*c2-n*b2))*(scipy.log((Vg2+n*bmg)/(Vg2+m*cmg)))/(R*T*(m*cmg-n*bmg)**2)+((amg/(R*T*(m*cmg-n*bmg)))*U10)-scipy.log(Zg2)
    
    '''Liquid Phase'''
    V10=((n*b2-m*c2)*Vl2+n*m*(cml*b2-bml*c2))/((Vl2+n*bml)*(Vl2+m*cml))
    lnphiL2=-(scipy.log((Vl2-bml)/Vl2))+(b2/(Vl2-bml))+((2*(1-x1)*a2)-aml*(m*c2-n*b2))*(scipy.log((Vl2+n*bml)/(Vl2+m*cml)))/(R*T*(m*cml-n*bml)**2)+((aml/(R*T*(m*cml-n*bml)))*V10)-scipy.log(Zl2)
    
    phiG2=scipy.exp(lnphiG2);   phiL2=scipy.exp(lnphiL2)
    
    '''PHI-PHI MODEL FOR STABILITY ANALYSIS'''    
    error=abs((y1*phiG1-x1*phiL1)+((1-y1)*phiG2-(1-x1)*phiL2))  #Error Function to be Minimized to Required Tolerance Level
    
    '''Modification to the Initial Guess Value'''
    x1=(phiG1*(phiG2-phiL2))/(phiL1*phiG2-phiG1*phiL2)
    y1=x1*phiL1/phiG1
    
    counter=counter+1
    
'''
print "Vapour Phase Fugacity of Component 1 in the Mixture = phiV1 =",phiG1
print "Vapour Phase Fugacity of Component 2 in the Mixture = phiV2 =",phiG2
print "Liquid Phase Fugacity of Component 1 in the Mixture = phiL1 =",phiL1
print "Liquid Phase Fugacity of Component 2 in the Mixture = phiL2 =",phiL2
''' 

print ""    
print "Liquid Phase Mole Fraction of Component 1 (Water) = x1 =",x1
print "Vapour Phase Mole Fraction of Component 1 (Water) = y1 =",y1
print "Final Error =",error

#Remaining Parameters:
x2=1-x1
y2=1-y1

"""
References:

    1. A new three-parameter cubic equation of state for calculation of 
        physical properties and vapour-liquid equilibria.
        
        - A. Haghtalab, M. J. Kamali, S. H. Mazloumi, P. Mahmoodi
        
    2. A simple and uniﬁed algorithm to solve ﬂuid phase
       equilibria using either the gamma-phi or the phi-phi
       approach for binary and ternary mixtures

        - Romain Privat, Jean-Noël Jaubert and Yannick Privat.
"""