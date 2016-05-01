# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:42:25 2016

@author: Sparsh Ganju
"""

import scipy
Temp=298.0 #Kelvin
Lin=100.0 #moles per s
Gin=100.0 #moles per s

""" 
1=METHANE ; 2=CARBON DIOXIDE
"""

kla1=0.002 # Overall Mass Transfer Coefficient for Methane
kla2=0.001 # Overall Mass Transfer Coefficient for Carbon-dioxide
kga_w=0.0005 # Overall Gas side Mass Transfer Coefficient for Water

n=10  # Grid Size
n_list=range(n-1)


"""Data"""
P=1.0 # Operating Pressure in 'atm'
Psat=3170.0/101325.0 # psat of water at 298K Pa
z=10.0 # Height of Column in 'metres'
dz=z/(n-1) # Step Size
A=10 # Area of Cross Section for Mass Transfer  
H1=714.28# Henry's Constant for Methane in 'lit*atm/mol'   
H2=29.4 # Henry's Constant for Carbon Dioxide in 'lit*atm/mol'
"""References are mentioned at the end of the program """


#flow Rates
L1=[None]*n
L2=[None]*n
G1=[None]*n
G2=[None]*n
Lw=[None]*n
Gw=[None]*n


G1[0]=50.0
G2[0]=50.0
Gw[0]=0.0

L1final=0.0
L2final=0.0
Lwfinal=100.0

k=0 # No of Steps
error=1

"""Guess Values"""

L1guess=0.5
L2guess=0.5
Lwguess=99.0


""" Using Euler's Explicit Method to solve the Differential Equations"""

while abs(error)>0.000001:
    L1[0]=L1guess
    L2[0]=L2guess
    Lw[0]=Lwguess
    
   
    for i in n_list: 
        
        
        
             
        L1[i+1] = L1[i] - (kla1*((P*G1[i]/(G1[i]+G2[i]+Gw[i]))/H1)-((L1[i]/Lw[i])*1000/18))*A*dz
        G1[i+1] = G1[i] - (kla1*((P*G1[i]/(G1[i]+G2[i]+Gw[i]))/H1)-((L1[i]/Lw[i])*1000/18))*A*dz
        
        L2[i+1] = L2[i] - (kla2*((P*G2[i]/(G2[i]+G2[i]+Gw[i]))/H2)-((L2[i]/Lw[i])*1000/18))*A*dz
        G2[i+1] = G2[i] - (kla2*((P*G2[i]/(G2[i]+G2[i]+Gw[i]))/H2)-((L2[i]/Lw[i])*1000/18))*A*dz
        
        Lw[i+1] = Lw[i] - (kga_w*((P*Gw[i]/(G1[i]+G2[i]+Gw[i])) - Psat))*A*dz
        Gw[i+1] = Gw[i] - (kga_w*((P*Gw[i]/(G1[i]+G2[i]+Gw[i])) - Psat))*A*dz
        
        
        
        
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
        
        
    
    
