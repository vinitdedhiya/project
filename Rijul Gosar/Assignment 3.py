#Mass transfer of CO2 and methane into water

import scipy
import numpy
import matplotlib.pyplot as plt

Temp=298.0 #Kelvin
L_rate= 100.0 #moles per hour
G_in=50.0 #moles per hour

kla_met= 0.00003 #m^2/s overall mass transfer coefficient for methane from Perry's Handbook pg 2113
kla_co2= 0.0007 #m^2/s overall mass transfer coeeficient for CO2 from Perry's Handbook pg 2113
kga=0.00003 #m^2/s overall gas side mass transfer coefficient for water from Perry's Handbook pg 2113

n=10 #grid points

xmetg=numpy.empty([n])
xco2g=numpy.empty([n])
xco2l=numpy.empty([n])
xmetl=numpy.empty([n])
xwl=numpy.empty([n])
xwg=numpy.empty([n])

xmetg[0]=0.5 #input of methane  
xco2g[0]=0.5 #input of CO2


xmetfinal=0.0 #input of methane in liquid

xco2final=0.0 #input of co2 in liquid


xwg[0]=0.0 #water vapour input in gas
xwfinal=1.0 #water in from top of column

psat=3290.0 #pascal
#http://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
H_met=7142.86 #Pa*m^3/mol; Henry's Constant for methane from Perry's Handbook pg 2113
H_co2=29.4118 #Pa*m^3/mol; Henry's Constant for carbon Dioxide from Perry's Handbook pg 2113
A=5.0 #Area of cross section
z=5.0 #height of column
p=101325.0 #operation pressure
dz=z/(n-1)

itr=0
err=1.0
err2=1.0
err3=1.0
err1=1.0

Cl_met=numpy.empty([n]) #Concentration of methane in liquid phase
Cl_co2=numpy.empty([n]) #Concentration of Carbon dioxide in liquid phase

#liquid phase concentration is the same amount transferred from the gas phase


G_met=numpy.empty([n]) #Gas phase flow rate of methane
G_co2=numpy.empty([n]) #Gas phase flow rate of Carbon Dioxide
L_met=numpy.empty([n]) #Liquid phase flow rate of Methane
L_co2=numpy.empty([n]) #Liquid phase flow rate of Carbon Dioxide
Gw=numpy.empty([n]) #Water vapour flow rate
Lw=numpy.empty([n]) #Water flow rate in liquid phase

xmetl[0]=0.1 #output of methane in liquid-assumption
xco2l[0]=0.2 #output of co2 in liquid-assumption
xwl[0]=0.7 #water output in liquid
         

z=numpy.empty([n])

for i in range(n):
    z[i]=i*dz
    
x=z

while abs(err)>0.001:
    
    G_met[0]=xmetg[0]*G_in #Methane gas flow rate at entry
    G_co2[0]=xco2g[0]*G_in #Carbon dioxide gas flow rate at entry
    L_met[0]=xmetl[0]*L_rate #Methane in liquid at bottom of tower
    L_co2[0]=xco2l[0]*L_rate #Carbon dioxide in liquid at bottom of tower
    Gw[0]=0 #Water vapour in gaseous phase at entry
    Lw[0]=xwl[0]*L_rate #Water in liquid at bottom  

    
    Cl_met[0]=xmetl[0]/xwl[0]*10**6/18 #Concentration of methane at bottom of tower
    Cl_co2[0]=xco2l[0]/xwl[0]*10**6/18 #Concentration of co2 at bottom of tower

#Iteration to estimate final mole fractions
    
    for i in range(n-1):
        
        G_met[i+1]=G_met[i]+kla_met*A*dz*(Cl_met[i]-xmetg[i]*p/H_met)
        G_co2[i+1]=G_co2[i]+kla_co2*A*dz*(Cl_co2[i]-xco2g[i]*p/H_co2)
        
        L_met[i+1]=L_met[i]+kla_met*A*dz*(Cl_met[i]-xmetl[i]*p/H_met)
        L_co2[i+1]=L_co2[i]+kla_co2*A*dz*(Cl_co2[i]-xco2l[i]*p/H_co2)
        
        Gw[i+1]=Gw[i]-kga*A*dz*(-psat+xwg[i]*p)
        Lw[i+1]=Lw[i]-kga*A*dz*(-psat+xwl[i]*p)
        
        xmetl[i+1]=L_met[i+1]/(L_met[i+1]+L_co2[i+1]+Lw[i+1])
        xco2l[i+1]=L_co2[i+1]/(L_met[i+1]+L_co2[i+1]+Lw[i+1])
        xwl[i+1]=Lw[i]/(L_met[i+1]+L_co2[i+1]+Lw[i+1])
        
        xmetg[i+1]=G_met[i+1]/(G_met[i+1]+G_co2[i+1]+Gw[i+1])
        xco2g[i+1]=G_co2[i+1]/(G_met[i+1]+G_co2[i+1]+Gw[i+1])
        xwg[i+1]=Gw[i+1]/(G_met[i+1]+G_co2[i+1]+Gw[i+1])
        
        Cl_met[i+1]=xmetl[i+1]/(xwl[i+1]*10**6/18)
        Cl_co2[i+1]=xco2l[i+1]/(xwl[i+1]*10**6/18)
        
    err1=xmetl[n-1]-xmetfinal
    err2=xco2l[n-1]-xco2final
    err3=xwl[n-1]-xwfinal
       
    err= err1**2+err2**2+err3**2
    
    
    xmetl[0]=xmetl[0]-err1/20.0 #correction for guess value
    xco2l[0]=xco2l[0]-err2/20.0 #correction to guess value
    xwl[0]=xwl[0]-err3/20.0 #correction to guess value

    #print itr
#print err1, err2, err3
#print G_met[n-1], L_met[n-1], G_co2[n-1], L_co2[n-1], Gw[n-1], Lw[n-1]

print 'xmetg='
print xmetg, "\n"
print 'xco2g='
print xco2g, "\n"
print 'xwg='
print xwg


plt.subplot(311)
plt.title('Mole fraction (in gas phase) vs Height')
#plt.xlabel('Height')
plt.ylabel('CO2')
plt.plot(x,xco2g,'r')
plt.show()

plt.subplot(312)
#plt.title('Mole fraction vs Height')
#plt.xlabel('Height')
plt.ylabel('Water')
plt.plot(x,xwg,'b')
plt.show()

plt.subplot(313)
#plt.title('Mole fraction vs Height')
plt.xlabel('Height')
plt.ylabel('Methane')
plt.plot(x,xmetg,'g')
plt.show()
