# -*- coding: utf-8 -*-
"""
Created on Sun Jan 03 23:23:27 2016

@author: Sparsh Ganju
"""
class HX:
    n = 10; #No. of nodes including ends

    def counter_HX(self):    
        #import scipy
        U = 100.0; #Overall htc
        P = 10.0; #Perimeter
        L = 10.0; #Length
        n = self.n; #No. of nodes including ends
        dx = L/(n-1);
        m_a = 1; #Mass Flow rate of A
        m_b = 2; #Mass flow rate of B
    
        Ta = [None]*n;
        Ta[0] = 400.0; #Inlet temp. of fluid A
    
        Tb = [None]*n;
        Tb[n-1] = 300.0;
        
        Tbresult = [None]*n;
        #Tbguess = [None]*1000;
        Tbguess = 350;
        
        e = 1;
        itr = 0; #iteration counter
        while abs(e)>0.0001:
            for i in range(n-1):
                Tbresult[0] = Tbguess;
                Cpa = 4000+0.1*Ta[i]+0.01*(Ta[i]**2);
                Cpb = 3000+0.2*Tbresult[i]+0.05*(Tbresult[i]**2);
                Ta[i+1] = Ta[i] - (P*U*dx*(Ta[i]-Tbresult[i]))/(m_a*Cpa);
                Tbresult[i+1] = Tbresult[i] - (P*U*dx*(Ta[i]-Tbresult[i]))/(m_b*Cpb);
            e = Tbresult[n-1] - Tb[n-1];
            Tbguess = Tbguess - (e/10);
            itr = itr + 1; #print itr;
    
        Tb = Tbresult;
        #print "The temperature of fluid A along the length of the HX is: ",Ta;
        #print "The temperature of fluid B along the length of the HX is: ",Tb;
    
        dh_A = m_a*(4000*(Ta[0]-Ta[n-1])+0.05*(Ta[0]**2-Ta[n-1]**2)+(0.01/3)*(Ta[0]**3-Ta[n-1]**3));
        dh_B = m_b*(3000*(Tb[0]-Tb[n-1])+0.1*(Tb[0]**2-Tb[n-1]**2)+(0.05/3)*(Tb[0]**3-Tb[n-1]**3));
        error_dh = dh_A-dh_B;
        error = abs(error_dh*100/dh_A);
        return error;
        
"""
For executing class on console and plotting a graph of error vs grid size the
following commands should be typed into the console:
"""
if __name__=="__main__":
    list_n = [2,5,10,20,30,40,50,100,250,500,750,1000]
    list_error = []
    inst1 = HX()
    for n in list_n:
        inst1.n = n
        error = inst1.counter_HX()
        list_error += [error]
        print n, error
        
    import matplotlib.pyplot as plt
    plt.plot(list_n,list_error,'b')
    plt.axis([-100,1000,-0.1,0.6])
    plt.ylabel('Percent Error')
    plt.xlabel('Grid Size')
    plt.show()
