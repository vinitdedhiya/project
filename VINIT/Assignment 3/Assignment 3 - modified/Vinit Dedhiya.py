# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:55:23 2016

@author: Vinit Dedhiya
"""
import scipy
Temp=298.0 #Kelvin
Lin=100.0 #litres per hour
Gin=100.0 #litres per hour

""" 
1=METHANE ; 2=CARBON DIOXIDE
"""

kla1=0.00003 # Overall Mass Transfer Coefficient for Methane
kla2=0.0007 # Overall Mass Transfer Coefficient for Carbon-dioxide
kga_w=0.00003 # Overall Gas side Mass Transfer Coefficient for Water

n=10  # Grid Size
n_list=range(n)

x1=[None]*n # Liquid Phase Mole Fraction Array for METHANE
x2=[None]*n # Liquid Phase Mole Fraction Array for CARBON DIOXIDE
y1=[None]*n # Gas Phase Mole Fraction Array for METHANE
y2=[None]*n # Gas Phase Mole Fraction Array for CARBON DIOXIDE
xw=[None]*n # Liquid Phase Water Composition
yw=[None]*n # Gas Phase Water(VAPOUR) Composition


x1[n-1]=0.0  # Input Liquid Stream Composition of Methane
x2[n-1]=0.0 # Input Liquid Stream Composition of Carbon Dioxide
y1[0]=0.5  # Input Gas Stream Composition of Methane
y2[0]=0.5  # Input Gas Stream Composition of Carbon Dioxide
xw[n-1]=1.0  # Input Liquid Stream Composition of Water - PURE WATER
yw[0]=0.0  # Input Gas Stream Composition of Water - WATER VAPOUR

"""
Data
"""
P=101325.0 # Operating Pressure in 'Pa'
Psat=3170.0 # psat of water at 298K Pa
z=10.0 # Height of Column in 'metres'
dz=z/(n-1) # Step Size
A=5 # Area of Cross Section for Mass Transfer  
H1=39200.0 # Henry's Constant for Methane    
H2=1635.0 # Henry's Constant for Carbon Dioxide
"""References are mentioned at the end of the program """

Cl_1=[None]*n
Cl_2=[None]*n
Cl_1guess=0.05 # Output Conc. of Methane in Liquid Stream
Cl_2guess=0.05 # Output Conc. of Carbon Dioxide in Liquid Stream
Pw=[None]*n
P1=[None]*n
P2=[None]*n

Cl_1[0]=Cl_1guess
Cl_2[0]=Cl_2guess


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

"""
Guess Values
"""
L1guess=5.0
L2guess=5.0
Lwguess=90.0
G1[0]=y1[0]*Gin; 
G2[0]=y2[0]*Gin;
Gw[0]=0
""" Using Euler's Explicit Method to solve the Differential Equations"""

while error>0.01:
    L1[0]=L1guess
    L2[0]=L2guess
    Lw[0]=Lwguess
    
    for i in n_list:  
        P1[i]=y1[i]*P
        P2[i]=y2[i]*P
        Pw[i]=yw[i]*P
        
        L1[i+1] = L1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A
        G1[i+1] = G1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A
        
        L2[i+1] = L2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A
        G2[i+1] = G2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A
        
        Lw[i+1] = Lw[i] - kga_w*(Pw[i] - Psat)*A
        Gw[i+1] = Gw[i] - kga_w*(Pw[i] - Psat)*A
        
        x1[i]=L1[i]/(L1[i]+L2[i]+Lw[i])
        x2[i]=L2[i]/(L1[i]+L2[i]+Lw[i])
        xw[i]=Lw[i]/(L1[i]+L2[i]+Lw[i])
        
        
    e1 = L1[n-1]-L1final # Errors
    e2 = L2[n-1]-L2final
    e3 = Lw[n-1]-Lwfinal
        
    L1guess=L1guess-e1/5 #correction to the guess Values
    L2guess=L2guess-e2/5
    Lwguess=Lwguess-e3/5
        
       
        
    error=e1**2+e2**2+e3**2; # Minimizing the Sum of Errors
                
    k=k+1; 
    
print G1[n-1]
print G2[n-1]
print Gw[n-1]        
"""
    References:
    1. Perrys handbook 1999
    page no 2113,2114,2115
    2. For Henrys constants
        Henry-3.0 pdf,R.Sander: Henrys law constant
"""
        
        
    
    
