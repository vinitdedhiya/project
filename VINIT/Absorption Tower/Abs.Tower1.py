# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:55:23 2016

@author: Vinit Dedhiya
"""

'''DATA'''

Temp=298.0 #Kelvin
Lin=100.0 #Unit : 'mol/s' - Input Flow Rate on Liquid Side
Gin=100.0 #Unit : 'mol/s' - Input Flow Rate on Gas Side

""" 
1=METHANE ; 2=CARBON DIOXIDE
"""
kla1=0.001 # Overall Mass Transfer Coefficient for Methane
kla2=0.002 # Overall Mass Transfer Coefficient for Carbon-dioxide
kga_w=0.0005 # Overall Gas side Mass Transfer Coefficient for Water

n=100 # Grid Size
n_list=range(n-1)


y1=[None]*n # Gas Phase Molar Fraction Array for METHANE
y2=[None]*n# Gas Phase Molar Fraction Array for CARBON DIOXIDE
yw=[None]*n # Gas Phase Water(VAPOUR) Composition

"""Fixed Values"""

y1[0]=0.5  # Input Gas Stream Composition of Methane
y2[0]=0.5  # Input Gas Stream Composition of Carbon Dioxide
yw[0]=0.0  # Input Gas Stream Composition of Water - WATER VAPOUR

"""Data"""
P=1.0 # Operating Pressure in 'atm'
Psat=3170.0/101325.0 # psat of water at 298K 'atm'
z=10.0 # Height of Column in 'metres'
dz=z/(n-1) # Step Size
A=10 # Area of Cross Section for Mass Transfer  
H1=714.28# Henry's Constant for Methane in 'lit*atm/mol'  
H2=29.4 # Henry's Constant for Carbon Dioxide in 'lit*atm/mol'
"""References are mentioned at the end of the program """

Cl_1=[None]*n
Cl_2=[None]*n

Pw=[None]*n
P1=[None]*n
P2=[None]*n

#flow Rates
L1=[None]*n
L2=[None]*n
G1=[None]*n
G2=[None]*n
Lw=[None]*n
Gw=[None]*n

L1final=0.0
L2final=0.0
Lwfinal=100.0

k=0 # No of Steps
error=1

"""Guess Values"""


L1guess=0.5
L2guess=0.5
Lwguess=99.0

G1[0]=50 
G2[0]=50
Gw[0]=0

""" Using Euler's Explicit Method to solve the Differential Equations"""

while error>0.000001:
    L1[0]=L1guess
    L2[0]=L2guess
    Lw[0]=Lwguess
    
   
    for i in n_list: 
        
               
             
        Cl_1[0]=(L1[0]/Lw[0])*1000/18
        Cl_2[0]=(L2[0]/Lw[0])*1000/18
        
        P1[i]=y1[i]*P
        P2[i]=y2[i]*P
        Pw[i]=yw[i]*P
        
        L1[i+1] = L1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A*dz
        G1[i+1] = G1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A*dz
        
        L2[i+1] = L2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A*dz
        G2[i+1] = G2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A*dz
        
        Lw[i+1] = Lw[i] - kga_w*(Pw[i] - Psat)*A*dz
        Gw[i+1] = Gw[i] - kga_w*(Pw[i] - Psat)*A*dz
        
        y1[i+1]=G1[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        y2[i+1]=G2[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        yw[i+1]=Gw[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        
        Cl_1[i+1]=(L1[i+1]/Lw[i+1])*1000/18
        Cl_2[i+1]=(L2[i+1]/Lw[i+1])*1000/18
        
        
        
        
    e1 = L1[n-1]-L1final # Errors
    e2 = L2[n-1]-L2final
    e3 = Lw[n-1]-Lwfinal
        
    L1guess=L1guess-e1/10 #correction to the guess Values
    L2guess=L2guess-e2/10
    Lwguess=Lwguess-e3/10
        
       
        
    error=abs(e1)+abs(e2)+abs(e3); # Minimizing the Sum of Errors
                
    k=k+1; 
    
z=[None]*n    
for i in range(n):
    z[i]=i*dz
    
import pylab as py

x=z
Y1=Cl_1
Y2=Cl_2

py.plot(x, Y1, '-r', label='Concentration of CH4')
py.plot(x, Y2, '-b', label='Concentration of CO2')
py.xlabel('Height of Absorber (m)')
py.ylabel('Concentration (M)')
py.title('Concentration (M) vs Height of Absorber (m)')
py.legend()
py.show()

    
    
"""
print G1
print G2
print Gw
print L1
print L2
print Lw
print y1
print y2
print yw
print Cl_1
print Cl_2     

    References:
    1. Perrys handbook 1999
    page no 2113,2114,2115
    2. For Henrys constants
        Henry-3.0 pdf,R.Sander: Henrys law constant
"""
        
        
    
    
