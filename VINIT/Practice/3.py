# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 02:01:44 2016

@author: Vinit Dedhiya
"""
import pylab
#Data
ma=1.0 #Mass Flow Rate of A
mb=2.0 #Mass Flow Rate of B
Tain=400.0 #Inlet flid Temp. of A
Tbin=300.0 #Inlet fluid Temp. of B
U=100.0 #Overall Heat Transfer Coefficient
P=1.0 #Perimeter of Rod
L=10.0 #Length of Rod
N=101#Grid Size List
n=range(N)
dx=L/(N-1)
k=0
error=10
while abs(error)>0.00001:
    Tax=Tain
    Tbx=Tbin+0.1*k
    for i in n:
        Ta=Tax;
        Tb=Tbx;
        def Ca(t):
            Ca= 4000.0 + 0.1*t+ 0.01*t*t
            return Ca
        def Cb(t):
            Cb=3000.0+0.2*t+0.05*t*t
            return Cb
        
        cpa=Ca(Ta)
        cpb=Cb(Tb)
        
        f=-((P*U)/(ma*cpa)*(Ta-Tb))
        g=-((P*U)/(mb*cpb)*(Ta-Tb))
        Tax=f*dx+Ta
        Tbx=g*dx+Tb
    error=Tb-Tbin
    k=k+1;
    print "Tb_out=",Tb;
    print "Ta_out=",Ta;
pylab.plot(range(0,100,10),Tb,label='Tb')
pylab.title('Temperature vs no. of grids')
pylab.show()    
