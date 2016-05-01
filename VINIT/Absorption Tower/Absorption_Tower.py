# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 22:55:35 2016

@author: Vinit Dedhiya
"""

"""Counter-Current Gas Absorption Tower"""

""" Objective: To obtain the Concentration Profile for the Absorption of 2
               Gases - CARBON DIOXIDE and METHANE along Height of the 
               Absorption Tower
"""
                
#DATA:
'''KEY:   1=METHANE ; 2=CARBON DIOXIDE ; w=WATER ; [0]=BOTTOM ; [n-1]=TOP '''
                
Temp=298.0 #Temperature in 'Kelvin'
P=1.0 #Operating Pressure in 'atm'

Lin=100.0 #Unit : 'mol/s' - Input Molar Flow Rate of Liquid from the TOP
Gin=100.0 #Unit : 'mol/s' - Input Molar Flow Rate of Gas from the BOTTOM

n=10 #Grid size
n_list=range(n-1)
Z=10.0 #Height of Tower in 'metres(m)'
A=10.0 #Area of Cross-Section for Maass Transfer in 'm^2'
dz=Z/(n-1) #Step Size

'''Mass Transfer Coefficients'''
kla1=0.002 # Overall Mass Transfer Coefficient for Methane on Liquid side
kla2=0.004 # Overall Mass Transfer Coefficient for Carbon-dioxide on Liquid side
kga_w=0.0005 # Overall Gas side Mass Transfer Coefficient for Water on Gas side
"""Reference:
             Perry's Chemical Engineers' Handbook (8th Edition)
             Table 19.10-19.11, Page 19.42 (pg.2118) """

'''Vapour Pressure of Water'''                
Psaturation=3170.0 #Saturation Pressure of Water at 298K in 'Pa'
#Conversion into the units used in the program:
Psat=3170.0/101325.0 #Psat of Water in 'atm'
"""Reference:
            The Antoine Co-efficients for Water have been taken from:
            Perry's Chemical Engineers' Handbook (8th Edition)
            Table 12.3, Page 12.5 (pg.1325) """

'''Henry's Costants'''
H1=714.28# Henry's Constant for Methane in 'lit*atm/mol'  
H2=29.41 # Henry's Constant for Carbon Dioxide in 'lit*atm/mol'
"""Reference:
            Compilation of Henry’s Law Constants for Inorganic and
            Organic Species of Potential Importance in Environmental Chemistry
            -Rolf Sander - (Version 3)"""
            
#VARIABLES:
            
'''Molar Flow Rates'''
#LIQUID:
L1=[None]*n #liquid Phase Molar Flow Rate of METHANE in 'mol/s'
L2=[None]*n #liquid Phase Molar Flow Rate of CARBON DIOXIDE in 'mol/s'
Lw=[None]*n #Molar Flow Rate of Liquid WATER in 'mol/s'

#GAS:
G1=[None]*n #Gas Phase Molar Flow Rate of METHANE in 'mol/s'
G2=[None]*n #Gas Phase Molar Flow Rate of CARBON DIOXIDE in 'mol/s'
Gw=[None]*n #Molar Flow Rate of WATER VAPOUR(gas) in 'mol/s'

'''Mole Fractions in the Gas Steam'''
y1=[None]*n #Gas Phase Molar Fraction Array for METHANE
y2=[None]*n #Gas Phase Molar Fraction Array for CARBON DIOXIDE
yw=[None]*n #Gas Phase WATER VAPOUR Composition

'''Concentration of Gases in the Liquid Stream'''
Cl_1=[None]*n #Concentration of dissolved METHANE in liquid stream in 'mol/lit'
Cl_2=[None]*n #Concentration of dissolved CARON DIOXIDE in liquid stream in 'mol/lit'

'''Partial Pressures in the Gas Stream'''
P1=[None]*n #Partial Pressure of METHANE in 'atm'
P2=[None]*n #Partial Pressure of CARBON DIOXIDE in 'atm'
Pw=[None]*n #Partial Pressure of WATER VAPOUR in 'atm'

#KNOWN VALUES:

y1[0]=0.5  # Input Gas Stream Composition of Methane
y2[0]=0.5  # Input Gas Stream Composition of Carbon Dioxide
yw[0]=0.0  # Input Gas Stream Composition of Water - WATER VAPOUR

'''In terms of Gas Flow Rates:'''
G1[0]=50 # 'mol/s'
G2[0]=50 # 'mol/s'
Gw[0]=00 # 'mol/s'

'''Liquid stream - VALUES AT THE TOP:'''
L1final=0.0 # 'mol/s'
L2final=0.0 # 'mol/s'
Lwfinal=100 # 'mol/s'

"""BRIEF DESCRIPTION OF THE METHOD USED:"""
'''-Euler's Explicit method is used to numerically solve the given system of
    Differential Equations 
    
   -We start solving by using a set of arbitrary Guess Values and calculate 
    the values at the other end
    
   -However, these values are already known(fixed) as this is a BVP
    So, we can calculate the Error by taking the difference between
    the known values and the calculated values
    
   -This Error value is used to change the Initial Guess Value in order to
    get a Final Value within the Range of Accuracy required
    
   -The Process is repeated until the final solution converges by changing
    the Guess Value after every Iteration [SHOOTING METHOD] '''
    
'''Guess Values:'''
L1guess=0.5 # 'mol/s' - Liquid Phase Flow Rate of METHANE at the BOTTOM in 'mol/s'
L2guess=0.5 # 'mol/s' - Liquid Phase Flow Rate of CARBON DIOXIDE at the BOTTOM in 'mol/s'
Lwguess=99 # 'mol/s' - WATER Flow Rate at the BOTTOM in 'mol/s'

"""EULER'S EXPLICIT METHOD to solve the System of Differential Eqations"""
k=0 #Number of Iterations
error=1 #Total Error Value


while error>0.0000001:
    L1[0]=L1guess
    L2[0]=L2guess
    Lw[0]=Lwguess
    
   
    for i in n_list: 
        
               
        '''Conversion Factor: 1000/18 = (55.55 M)'''     
        Cl_1[0]=(L1[0]/Lw[0])*1000/18 #Initial Concentration of METHANE in 'mol/lit'
        Cl_2[0]=(L2[0]/Lw[0])*1000/18 #Initial Concentration of CARBON DIOXIDE in 'mol/lit'
        
        P1[0]=y1[0]*P #Initial Partial Pressure of METHANE in 'atm'
        P2[0]=y2[0]*P #Initial Partial Pressure of CARBON DIOXIDE in 'atm'
        Pw[0]=yw[0]*P #Initial Partial Pressure of WATER VAPOUR in 'atm'
        
        '''Solving DIfferential Equations for each of the 3 Components'''
        #METHANE:
        L1[i+1] = L1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A*dz
        G1[i+1] = G1[i] - kla1*((P1[i]/H1)-(Cl_1[i]))*A*dz
        
        #CARBON DIOXIDE:
        L2[i+1] = L2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A*dz
        G2[i+1] = G2[i] - kla2*(P2[i]/H2-(Cl_2[i]))*A*dz
        
        #WATER:
        Lw[i+1] = Lw[i] - kga_w*(Pw[i] - Psat)*A*dz
        Gw[i+1] = Gw[i] - kga_w*(Pw[i] - Psat)*A*dz
        
        #Gas Phase Mole Fractions:
        y1[i+1]=G1[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        y2[i+1]=G2[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        yw[i+1]=Gw[i+1]/(G1[i+1]+G2[i+1]+Gw[i+1])
        
        Cl_1[i+1]=(L1[i+1]/Lw[i+1])*1000/18 #Concentration of METHANE in 'mol/lit'
        Cl_2[i+1]=(L2[i+1]/Lw[i+1])*1000/18 #Concentration of CARBON DIOXIDE in 'mol/lit'
        
        P1[i+1]=y1[i+1]*P #Partial Pressure of METHANE in 'atm'
        P2[i+1]=y2[i+1]*P #Partial Pressure of CARBON DIOXIDE in 'atm'
        Pw[i+1]=yw[i+1]*P #Partial Pressure of WATER VAPOUR in 'atm'
        
    '''ERRORS:'''       
    e1 = L1[n-1]-L1final 
    e2 = L2[n-1]-L2final
    e3 = Lw[n-1]-Lwfinal
        
    '''Correction to the Guess Values:'''
    L1guess=L1guess-e1/10 
    L2guess=L2guess-e2/10
    Lwguess=Lwguess-e3/10
        
            
    error=abs(e1)+abs(e2)+abs(e3); # Minimizing the Total Error
                
    k=k+1; #Counter
    
"""PLOTTING OF GRAPH - Concentration Profile"""

z=[None]*n  #Step-wise Height of Tower   
for i in range(n):
    z[i]=i*dz
    
'''Using PYLAB to plot the Concentration Profiles'''    
import pylab as py 

x=z
Y1=Cl_1 #Profile of METHANE
Y2=Cl_2 #Profile of CARBON DIOXIDE

py.plot(x, Y1, '-r', label='Concentration of CH4')
py.plot(x, Y2, '-b', label='Concentration of CO2')
py.xlabel('Height of Tower (m)')
py.ylabel('Concentration (mol/lit)')
py.title('Concentration (mol/lit) vs Height of Tower (m)')
py.legend()
py.show()

print 'The Final Error is ',error;
print 'Number of Iterations =',k;

"""
import numpy as np
result= np.column_stack((Cl_1,Cl_2));
print 'The values of Concentrations from Bottom to Top are:',result
"""