# -*- coding: utf-8 -*-
"""
Created on Tue Jan 05 01:51:19 2016

@author: Sagar
"""

import scipy as sc
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pylab
Mh=1
Mc=2
Thin=450 #K
Tcin=350#K
U=300
A=100
n=10
Thguess=n*[Thin]
Tcguess=n*[Tcin]
Tguess=sc.array(Thguess+Tcguess)
def Cph(T):
    Cp=4000+10*T+10**(-2)*T**(2)#4184+10**(-4)*T+10**(-6)*T**2+10**(-9)*T**3
    return Cp
def Cpc(T):
    Cp=3000+5*T+2*10**(-2)*T**(2)#4184+10**(-4)*T+10**(-6)*T**2+10**(-9)*T**3
    return Cp    
#the array whose elements you want to modify is the first argument to the function
def residuals(T,U,A,Thin,Tcin,Mh,Mc):
    n=len(T)
    Th=T[:n/2]
    Tc=T[n/2:]
    dA=A/((n/2)-1)
    residualHL=U*(Thin-Tc[0])/(Mh*Cph(Thin))+(Th[1]-Thin)/dA
    residualCL=U*(Thin-Tc[0])/(Mc*Cpc(Tc[0]))+(Tc[1]-Tc[0])/dA
    residualHR=U*(Th[-1]-Tcin)/(Mh*Cph(Th[-1]))+(Th[-1]-Th[-2])/dA
    residualCR=U*(Tc[-1]-Tcin)/(Mh*Cpc(Tcin))+(Tcin-Tc[-2])/dA
    residualH=sc.zeros(n/2)
    residualC=sc.zeros(n/2)
    residualH[0]=residualHL;residualH[-1]=residualHR
    residualC[0]=residualCL;residualC[-1]=residualCR
    residualH[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mh*Cph(Th[1:-1]))+(Th[2:]-Th[1:-1])/dA
    #residualH[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mh*Cph(Th[1:-1]))+(Th[2:]-Th[0:-2])/dA for central difference
    residualC[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mc*Cpc(Tc[1:-1]))+(Tc[2:]-Tc[1:-1])/dA
    return sc.concatenate((residualH,residualC))

soln=opt.leastsq(residuals,Tguess,args=(U,A,Thin,Tcin,Mh,Mc))
Tsoln=soln[0]
Thsoln=Tsoln[:20/2];Thsoln[0]=Thin
Tcsoln=Tsoln[20/2:];Tcsoln[-1]=Tcin
print "Thsoln=",Thsoln 
print "Tcsoln=",Tcsoln
pylab.plot(range(0,100,10),Thsoln,label='Thsoln')
pylab.plot(range(0,100,10),Tcsoln,label='Tcsoln')
plt.legend(loc='lower left')
pylab.title('Temperature vs no. of grids')
pylab.show()    
    
    

