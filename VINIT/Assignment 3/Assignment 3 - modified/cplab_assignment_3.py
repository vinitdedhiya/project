class asbsorption:
    Kgc=0.1 #mass transfer coeff of CO2 kmol/hr.m3.kPa
    Kgm=0.1 #mass transfer coeff of CO2 kmol/hr.m3.kPa
    Kgw=0.1 #mass transfer coeff of water kmol/hr.m3.kPa
    H1=29.4 #listed henry's const for co2 1/atm
    H2=714.0  # listed henry'sconst for methane 1/atm
    Gtotal=100.0#basis is 100moles/hr
    Ltotal=100.0# basis in 100moles/hr
    Gm=50.0 # given flow rate of methane in moles/hr at the bottom
    Gc=50.0 # given flow rate of CO2 in moles/hr at the bottom
    Gw=0.0 # given flow rate of water in moles/hr at the bottom
    Z=10.0 # height of column in metres
    A=1 # area of cross section in m2
    Lw0=96.0 # assumed flow rate of water in moles/hr at bottom
    Lc0=2.0 # assumed flow rate of CO2 in moles/hr at bottom
    Lm0=2.0 # assumed flow rate of methane in moles/hr at bottom
    Lmn=0.0 # flow rate of methane at top
    Lcn=0.0 # flow rate of CO2 at the top
    Lwn=100.0 # flow rate of water at top
    P=101325.0#Pa
    Psat=3170.0# psat of water at 298K Pa
     
     
    def Hc(x):
         Hc=(x/101325)
         return Hc #convert to si units
      
    '''def Hm(y):
         Hm=(y/101325)
         return Hm #convert to si units
     print H2'''
    
    n=10
    dz=Z/(n-1)
    n_list=range(n)
    e=1
     
    while e>0.001:
         Lw=Lw0
         Lc=Lc0
         Lm=Lm0
         
         for i in n_list:
             xc=Lc/(Lc+Lw+Lm)
             yc=Gc/(Gc+Gm+Gw)
             fc=Kgc*A*(xc-(P*yc)/(Hc(H1)))
             Lc=Lc-fc*dz
             Gc=Gc-fc*dz
             
             xm=Lm/(Lc+Lw+Lm)
             ym=Gm/(Gc+Gm+Gw)
             fm=Kgc*A*(xm-(P*ym)/(Hc(H2)))
             Lm=Lm-fm*dz
             Gm=Gm-fm*dz
             
             xw=Lw/(Lc+Lw+Lm)
             yw=Gw/(Gc+Gm+Gw)
             fw=Kgw*A*(P*yw-Psat)
             Lw=Lw+fw*dz
             Gw=Gw+fw*dz
          
         e=abs(Lw-Lwn)+abs(Lc-Lcn)+abs(Lm-Lmn) 
         print e
         Lm0=Lm-0.1*(Lm-Lmn)
         Lc0=Lc-0.1*(Lc-Lmn)
         Lw0=Lw-0.1*(Lw-Lwn)
        
    print(Lw0)# water flow rate at top
    print(Lc0)# CO2 flow rate at top
    print(Lm0) #methane flow rate at top
    print(Gm) # methane flow rate(gas) at top
    print(Gc) # CO2 flow rate(gas) at top
    print(Gw) # water flow rate(vapor) at top
"""
    References:
    1. Perrys handbook 1999
    page no 2114,2115
    2. For Henrys constants
        Henry-3.0 pdf,R.Sander: Henrys law constant
"""
          
          
          
             
             
             
             
             
         
         
     
     
        
         
     
    
     
         
         
         
         