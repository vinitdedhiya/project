''' The following program calculates the various stream compositions of a 
counter current gas absorption tower operating under steady state.'''

'''As a general rule, SI units have been used for all calculations.'''

'''The various input parameters are listed below with their respective units in 
comments.'''
class absorber():
    L = 10.0; #Length of absorber (m)
    A = 1.0; #Cross- sectional area of absorber (m2)
    T = 298.15; #Temperature (K)
    P = 101325.0; #Pressure (Pa)
    n = 100; #Number of elements in pipe
    
      
    def counter_abs(self):    
        L = self.L; #Length of absorber (m)
        A = self.A; #Cross- sectional area of absorber (m2)
        T = self.T; #Temperature (K)
        P = self.P; #Pressure (Pa)
        n = self.n; #Number of elements in pipe
        
        dz = L/(n-1); #Length of each element (m)
        V_water = 18*1e-3; #Volume of 1 kmol liquid water
        
        '''OVERALL MASS TRANSFER COEFFICIENTS (1/s) * INTERFACIAL AREA (m2/m3)
        The kla values have been obtained from pg. 19-42, Table 19-10, Perry's Chemical 
        Engineers' Handbook (8th Edition).
        The values given below are for representative purposes only. Actual values may 
        be calculated from the various correlations available.'''
        kla_co2 = 0.002; #for CO2 based on liquid side
        kla_meth = 0.001; #for Methane based on liquid side
        kga = 0.0005; #for water vapour based on gas side
        
        
        
        '''HENRY'S CONSTANTS (atm)
        The Henry's constants have been obtained from pg. 2-130, Table 2-123 (for CO2) 
        and Table 2-124 (for CH4), Perry's Chemical Engineers' Handbook (Edition 8).'''
        from H2SI import H2SI
        H_meth = 36600; #atm
        H_co2 = 1635 #atm
        #conversion to SI units
        H_meth = H2SI(H_meth); 
        H_co2 = H2SI(H_co2); 
        
        
        
        '''The antoine coefficients for water have been taken from pg. 682, Appendix B, 
        Table B2, Chemical Engineering Thermodynamics, Smmith Van Ness (Edition 7). 
        The saturation vapour pressure returned is in Pascal.'''
        from P_sat_water import P_sat_water
        P_sat_water = P_sat_water(T); #Pascal
        
        
        
        '''OVERALL MOLAR FLOW RATE OF GAS 
        The units of flow rate are taken to be kmol/s. The gas phase contains CO2, CH4 and water vapour.
        The mole fractions of the various components of the gas phase have been subsequently mentioned.'''
        G_gas = [None]*n; #Initialising an array of molar flow rates
        G_gas[0] = 100; #Gas Flow Rate at bottom of absorber (kmol/s);
        
        
        
        '''The various parameters (element-wise) for the various cpmponents have been 
        defined in the following lines of code.
        These parameters include: 
        1) Gas molar Flow Rate (kmol/s)
        2) Mole Fraction in gas phase
        3) Partial Pressures in gas phase (Pa)
        3) Concentration in water (kmol/m3)
        4) Molar flow rate of dissolved component'''
        
        '''METHANE'''
        G_meth = [None]*n; #Initialising an array of molar flow rates of CH4 (kmol/s)
        X_meth_gas = [None]*n; #Mole fraction of CH4 in gas phase
        X_meth_gas[0] = 0.5; #Inlet mole fraction of CH4 in gas phase
        P_meth = [None]*n; #Partial Pressure of CH4 (Pa)
        C_meth = [None]*n; #Concentration of CH4 in liquid phase (kmol/m3)
        C_meth[n-1] = 0; #Concentration of CH4 in liquid phase at entry (kmol/m3)
        L_meth = [None]*n; #Molar flow rate of dissolved CH4 in water (kmol/s)
        
        '''CARBON DIOXIDE'''
        G_co2 = [None]*n; #Initialising an array of molar flow rates of CH4 (kmol/s)
        X_co2_gas = [None]*n; #Mole fraction of CO2 in gas phase
        X_co2_gas[0] = 0.5; #Inlet mole fraction of CO2 in gas phase
        P_co2 = [None]*n; #Partial Pressure of CO2 (Pa)
        C_co2 = [None]*n; #Concentration of CO2 in liquid phase (kmol/m3)
        C_co2[n-1] = 0; #Concentration of C02 in liquid phase at entry (kmol/m3)
        L_co2 = [None]*n; #Molar flow rate of dissolved C02 in water (kmol/s)
        
        '''WATER AND WATER VAPOUR'''
        '''There is no such thing as a concentration driving force for water.'''
        G_water = [None]*n; #Initialising an array of molar flow rates of water vapour (kmol/s)
        X_water = [None]*n; #Mole fraction of water vapour in gas phase
        X_water[0] = 0; #Inlet mole fraction of water vapour in gas phase
        P_water = [None]*n; #Initialising an array of partial pressures of water vapour (Pa)
        L_water = [None]*n; #Initialising an array of molar flow rates of liquid water (kmol/s)
        L_water[n-1] = 100; #Input molar flow rate of liquid water (kmol/s)
        
        
        
        '''ALGORITHM FOR SOLVING THE DIFFERENTIAL EQUATIONS OF RATE OF MASS TRANSFER
        1) Euler's explicit method is used to discretize the equations. This is a BVP.
        2) Guess values of exit liquid flow rates for all the components are taken.
        3) The input values are back calculated and the error is obtained.
        4) This obtained is used to obtain a new value of guess value.
        5) The above procedure is repeated until the solution has converged (SHOOTING METHOD).'''
        
        
        '''GUESS VALUES FOR EXIT LIQUID FLOW RATES'''
        L_meth_guess = 0.05; #Guess value for dissolved CH4 (kmol/s)
        L_co2_guess = 0.01; #Guess value for dissolved CO2 (kmol/s)
        L_water_guess = 95; #Guess value for liquid water (kmol/s)
        
        '''ARRAYS FOR STORING BACK CALCULATED VALUES OF LIQUID MOLAR FLOW RATES'''
        L_meth_result = [None]*n; #For dissolved CH4 (kmol/s)
        L_co2_result = [None]*n; #For dissolved CO2 (kmol/s)
        L_water_result = [None]*n; #For liquid water (kmol/s)
        
        '''ARRAYS FOR STORING THE CALCULATED VALUES OF DISSOLVED CONCENTRATIONS'''
        C_co2_result = [None]*n; #For dissolved CO2 (kmol/m3)
        C_meth_result = [None]*n; #For dissolved CH4 (kmol/m3)
        
        '''DEFINING GUESS ERROR LIMITS BASED ON TOP LIQUID MOLAR FLOW RATES'''
        e1 = 1; #For dissolved CO2
        e2 = 1; #For dissolved CH4
        e3 = 1; #For liquid water
        err = abs(e1) + abs(e2) + abs(e3); #Total error
        
        itr = 0; #iteration counter
        while err>0.0000001:
            for i in range(n-1):
                L_water_result[0] = L_water_guess;
                L_meth_result[0] = L_meth_guess;
                L_co2_result[0] = L_co2_guess;
                
                '''Calculation of exit concentration of dissolved CO2'''   
                C_co2_result[0] = L_co2_result[0]/(L_water_result[0]*V_water);
                '''Calculation of exit concentration of dissolved CH4'''
                C_meth_result[0] = L_meth_result[0]/(L_water_result[0]*V_water);
                
                '''Molar gas flow rates (kmol/s)'''        
                G_co2[i] = G_gas[i]*X_co2_gas[i]; #CO2
                G_meth[i] = G_gas[i]*X_meth_gas[i]; #CH4
                G_water[i] = G_gas[i]*X_water[i]; #Water vapour
                
                '''Partial Pressures (Pa)'''
                P_meth[i] = P*X_meth_gas[i]; #CH4
                P_co2[i] = P*X_co2_gas[i]; #CO2        
                P_water[i] = P*X_water[i]; #Water vapour
                
                '''Solving mass transfer equation for CO2  (both liquid as well as gas)'''      
                G_co2[i+1] = G_co2[i] + kla_co2*A*dz*(C_co2_result[i]- (P_co2[i]/H_co2)); #Gas
                L_co2_result[i+1] = L_co2_result[i] + kla_co2*A*dz*(C_co2_result[i]- (P_co2[i]/H_co2)); #Liquid
                
                '''Solving mass transfer equation for CH4 (both liquid as well as gas)'''
                G_meth[i+1] = G_meth[i] + kla_meth*A*dz*(C_meth_result[i]- (P_meth[i]/H_meth)); #Gas
                L_meth_result[i+1] = L_meth_result[i] + kla_meth*A*dz*(C_meth_result[i]- (P_meth[i]/H_meth)); #Liquid
                
                '''Solving mass transfer equation for Water and Water Vapour'''
                G_water[i+1] = G_water[i] + kga*A*dz*(P_sat_water-P_water[i]); #Water
                L_water_result[i+1] = L_water_result[i] + kga*A*dz*(P_sat_water-P_water[i]); #Water Vapour
                
                '''Calculating Values of dissolved concentrations of CO2 and CH4'''
                C_co2_result[i+1] = L_co2_result[i+1]/(L_water_result[i+1]*18*1e-3); #CO2
                C_meth_result[i+1] = L_meth_result[i+1]/(L_water_result[i+1]*18*1e-3); #CH4
                
                '''Total molar gas flow rates at i+1 element'''
                G_gas[i+1] = G_co2[i+1] + G_meth[i+1] + G_water[i+1];
                
                '''Calculating gas phase mole fractions for i+1 element'''
                X_co2_gas[i+1] = G_co2[i+1]/G_gas[i+1]; #CO2
                X_meth_gas[i+1] = G_meth[i+1]/G_gas[i+1]; #CH4
                X_water[i+1] = G_water[i+1]/G_gas[i+1]; #Water vapour
               
            '''Calculating errors in outlet liquid flow rates'''
            e1 = C_co2_result[n-1]-C_co2[n-1]; #CO2
            e2 = C_meth_result[n-1]-C_meth[n-1]; #CH4
            e3 = L_water_result[n-1]-L_water[n-1]; #Liquid Water
            err = abs(e1) + abs(e2) + abs(e3); #Total error (always positive)
            
            '''Error Generating/ Changing Functions'''
            L_co2_guess = L_co2_guess - (e1/10); #CO2
            L_meth_guess = L_meth_guess - (e2/10); #CH4
            L_water_guess = L_water_guess - (e3/10); #Liquid Water 
            
            itr = itr + 1;
        
        '''FINAL LIQUID MOLAR FLOW RATES'''
        from transpose import transpose
        L_co2 = transpose(L_co2_result);
        L_meth = transpose(L_meth_result);
        L_water = transpose(L_water_result);
        
        '''FINAL DISSOLVED CONCENTRATIONS'''
        C_co2 = transpose(C_co2_result);
        C_meth = transpose(C_meth_result);
        
        '''EXIT CONCENTRATIONS'''
        import numpy as np
        exit = np.column_stack((C_co2_result,C_meth));
        #The height of the column is assumed from the ground level, i.e., element 1 is closest to ground.
        print 'The concentrations of CO2 and CH4 dissolved in water from bottom to top respectively are (in kmol/m3): \n'
                
        '''PLOTTING VARIOUS GRAPHS'''
        z = [None]*n;
        for i in range(n):
            z[i] = i*dz;
        
        '''PLOTTING DATA'''
        import matplotlib.pyplot as plt
        plt.figure('Concentration Profiles')
        plt.subplots_adjust(left = 0.2,right = 0.95,wspace = 1, hspace = 0.5)
        
        plt.subplot(211)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.title('Concentration of CO2 vs. Height');
        plt.xlabel('Height (m)');
        plt.ylabel('Concentration of CO2 (kmol/m3)');
        plt.plot(z,C_co2_result,'b')
        plt.show()
        
        plt.subplot(212)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.title('Concentration of CH4 vs. Height');
        plt.xlabel('Height (m)');
        plt.ylabel('Concentration of CH4 (kmol/m3)');
        plt.plot(z,C_meth_result,'r')
        plt.show()
        
        return exit, G_water