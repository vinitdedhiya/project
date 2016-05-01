# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 00:05:58 2016

@author: Vinit Dedhiya
"""
class ThinRod:
#Data
    ma=1 #Mass Flow Rate of A
    mb=2 #Mass Flow Rate of B
    Tain=400 #Inlet flid Temp. of A
    Tbin=300 #Inlet fluid Temp. of B
    U=100 #Overall Heat Transfer Coefficient
    P=1 #Perimeter of Rod
    L=10 #Length of Rod
    N=5 #Grid Size
    
    def Hx(self):
        ma=self.ma
        mb=self.mb
        Tain=self.Tain
        Tbin=self.Tbin
        U=self.U
        P=self.P
        L=self.L
        N=self.N
        dx=L/(N-1)
        n=range(N)
        k=0
        error=1
        while abs(error)>0.0001:
            Tax=Tain
            Tbx=Tbin+0.0001*k
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
        
                Tax=-((P*U)/(ma*cpa)*(Ta-Tb))*dx+Ta
                Tbx=-((P*U)/(mb*cpb)*(Ta-Tb))*dx+Tb
                
            error=Tb-Tbin
            k=k+1;
        print "Tb_out=",Tb;
        print "Ta_out=",Ta;
        
    
    