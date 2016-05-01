# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import scipy

class HeatExchanger:
      ma=1.0 #mass flow rate of A
      mb=2.0 #mass flow rate of B
      U=100.0 #overall heat transfer coefficient
      P=1.0 #Perimeter of heat exchanger
      L=10.0 #length 
      Tain=400 #inlet temperature of A
      Tbin=300 
      n=10 #grid points
          
      def Heat(self):
          ma=self.ma
          mb=self.mb
          U=self.U
          P=self.P
          L=self.L
          Tain=self.Tain
          #Tbin=self.Tbin
          n=self.n
          deltal=L/(n-1)
          n_list=range(n)
          itr = 0;
          Tguess = 340;
          err=1
          while err>0.0001:
              Ta1=Tain
              print itr
              Tb=Tguess
              for i in n_list:
                  Ta=Ta1
                  
                  Cpa= 4000.0 + 0.1*Ta1+ 0.01*Ta1*Ta1
                  Cpb= 3000.0 +0.2*Tb + 0.05*Tb*Tb
                  dta=-((U*P*L)*(Ta-Tb)/(ma*Cpa))
                  dtb=((U*P*L)*(Ta-Tb)/(mb*Cpb))
                  Ta1=(dta*deltal)+Ta
                  Tb=-(dtb*deltal)+Tb
              err=scipy.absolute(Tb-300)   
              itr = itr + 1;
              Tguess = Tguess - (err/10.0);
              #print Tguess,Tb
          print "outlet temperature of B", Tguess,dta,dtb,deltal,;
          print "outlet temperature of A", Ta;
          
          def deltaHa(T):
              deltaHa=(4000*(400-T)+0.05*(400*400-T*T)+(0.01*(400**3-T**3)/3))*ma
              return deltaHa
          def deltaHb(T):
              deltaHb=mb*(3000*(T-300)+0.1*(T*T-300*300)+(0.05*(T*T*T-300*300*300)/3))
              return deltaHb
          Ha=deltaHa(Ta)
          Hb=deltaHb(Tguess)
          print "Enthalpy of A",Ha;
          print "Enthalpy of B",Hb;
          return Tb