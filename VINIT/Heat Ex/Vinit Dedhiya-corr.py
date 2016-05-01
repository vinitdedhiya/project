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
m1=0.016 #Molar mass of Methane in 'kg/mol'
m2=0.044 #Molar mass of Carbon dioxide in 'kg/mol'
mw=0.018 #Molar mass of Water in 'kg/mol'
r1=0.000656 #Density of Methane in 'kg/lit'
r2=0.00198 #Density of Carbon dioxide in 'kg/lit'
rw=1 #Density of Water in 'kg/lit'

kla1=0.002 # Overall Mass Transfer Coefficient for Methane
kla2=0.001 # Overall Mass Transfer Coefficient for Carbon-dioxide
kga_w=0.0005 # Overall Gas side Mass Transfer Coefficient for Water

n=10  # Grid Size
n_list=range(n-1)

x1=[None]*n # Liquid Phase Volume Fraction Array for METHANE
x2=[None]*n # Liquid Phase Volume Fraction Array for CARBON DIOXIDE
y1=[None]*n # Gas Phase Volume Fraction Array for METHANE
y2=[None]*n# Gas Phase Volume Fraction Array for CARBON DIOXIDE
xw=[None]*n # Liquid Phase Water Composition
yw=[None]*n # Gas Phase Water(VAPOUR) Composition

"""Fixed Values"""
x1[n-1]=0.0  # Input Liquid Stream Composition of Methane
x2[n-1]=0.0 # Input Liquid Stream Composition of Carbon Dioxide
y1[0]=0.5  # Input Gas Stream Composition of Methane
y2[0]=0.5  # Input Gas Stream Composition of Carbon Dioxide
xw[n-1]=1.0  # Input Liquid Stream Composition of Water - PURE WATER
yw[0]=0.0  # Input Gas Stream Composition of Water - WATER VAPOUR

"""Data"""
P=1.0 # Operating Pressure in 'atm'
Psat=3170.0/101325.0 # psat of water at 298K Pa
z=10.0 # Height of Column in 'metres'
dz=z/(n-1) # Step Size
A=10 # Area of Cross Section for Mass Transfer  
H1=714.28# Henry's Constant for Methane in 'lit*atm/mol'   
H2=29.4 # Henry's Constant for Carbon Dioxide in 'lit*atm/mol'
"""References are mentioned at the end of the program """

Cl_1=[None]*n
Cl_2=[None]*n
#Cl_1[n-1]=0.0
#Cl_2[n-1]=0.0

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
x1[0]=0.05 # Output Volume fraction of Methane in Liquid Stream
x2[0]=0.05 # Output Volume fraction of Carbon Dioxide in Liquid Stream
xw[0]=1-x1[0]-x2[0] #Output Volume fraction of Water in Liquid Stream

L1guess=x1[0]*Lin
L2guess=x2[0]*Lin
Lwguess=xw[0]*Lin

G1[0]=y1[0]*Gin 
G2[0]=y2[0]*Gin
Gw[0]=yw[0]*Gin

""" Using Euler's Explicit Method to solve the Differential Equations"""

while error>0.001:
    L1[0]=L1guess
    L2[0]=L2guess
    Lw[0]=Lwguess
    
   
    for i in n_list: 
        
        x1[i]=L1[i]/(L1[i]+L2[i]+Lw[i])
        x2[i]=L2[i]/(L1[i]+L2[i]+Lw[i])
        xw[i]=Lw[i]/(L1[i]+L2[i]+Lw[i])
        
        y1[i]=G1[i]/(G1[i]+G2[i]+Gw[i])
        y2[i]=G2[i]/(G1[i]+G2[i]+Gw[i])
        yw[i]=Gw[i]/(G1[i]+G2[i]+Gw[i])
        
        Cl_1[i]=(L1[i]/Lw[i])*18/1000
        Cl_2[i]=(L2[i]/Lw[i])*18/1000
        
        P1[i]=y1[i]*P
        P2[i]=y2[i]*P
        Pw[i]=yw[i]*P
        
        L1[i+1] = L1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A
        G1[i+1] = G1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A
        
        L2[i+1] = L2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A
        G2[i+1] = G2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A
        
        Lw[i+1] = Lw[i] - kga_w*(Pw[i] - Psat)*A
        Gw[i+1] = Gw[i] - kga_w*(Pw[i] - Psat)*A
        
        
        
        
    e1 = L1[n-1]-L1final # Errors
    e2 = L2[n-1]-L2final
    e3 = Lw[n-1]-Lwfinal
        
    L1guess=L1guess-e1/10 #correction to the guess Values
    L2guess=L2guess-e2/10
    Lwguess=Lwguess-e3/10
        
       
        
    error=e1**2+e2**2+e3**2; # Minimizing the Sum of Errors
                
    k=k+1; 
    
"""
print G1
print G2
print Gw
print L1
print L2
print Lw
print x1  
print x2
print xw
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
        
        
    
    
