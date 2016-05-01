# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 21:55:23 2016

@author: Vinit Dedhiya
"""
import scipy
Temp=298 #Kelvin
Lin=100 #litres per hour
Gin=100 #litres per hour

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


x1[n-1]=0  # Input Liquid Stream Composition of Methane
x2[n-1]=0  # Input Liquid Stream Composition of Carbon Dioxide
y1[0]=0.5  # Input Gas Stream Composition of Methane
y2[0]=0.5  # Input Gas Stream Composition of Carbon Dioxide
xw[n-1]=1.0  # Input Liquid Stream Composition of Water - PURE WATER
yw[0]=0  # Input Gas Stream Composition of Water - WATER VAPOUR

"""
Data
"""
P=101325.0 # Operating Pressure in 'Pa'
Psat=3170.0 # psat of water at 298K Pa
z=10 # Height of Column in 'metres'
dz=z/(n-1) # Step Size
A=5 # Area of Cross Section for Mass Transfer  
H1=39200 # Henry's Constant for Methane
H2=1635 # Henry's Constant for Carbon Dioxide

Cg_1=[None]*n
Cg_2=[None]*n
Cl_1=[None]*n
Cl_2=[None]*n
Cl_1guess=0.05 # Output Conc. of Methane in Liquid Stream
Cl_2guess=0.05 # Output Conc. of Carbon Dioxide in Liquid Stream
Pw=[None]*n

L1=[None]*n
L2=[None]*n
G1=[None]*n
G2=[None]*n
Lw=[None]*n
Gw=[None]*n

k=0 # No of Steps
error=1

"""
Guess Values
"""
L1guess=5
L2guess=5
Lwguess=90

while error>0.01:
    
    for i in n_list:
        L1[0]=L1guess
        L2[0]=L2guess
        Lw[0]=Lwguess
        
        Ta[i+1] = Ta[i] - (P*U*dx*(Ta[i]-Tbresult[i]))/(m_a*Cpa);
        Tbresult[i+1] = Tbresult[i] - (P*U*dx*(Ta[i]-Tbresult[i]))/(m_b*Cpb);
        e = Tbresult[n-1] - Tb[n-1];
        Tbguess = Tbguess - (e/10);
        
        k=k+1; 
    
        Tb = Tbresult;
        
"""
    References:
    1. Perrys handbook 1999
    page no 2113,2114,2115
    2. For Henrys constants
        Henry-3.0 pdf,R.Sander: Henrys law constant
"""
        
        
    
    
