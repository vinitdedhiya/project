# -*- coding: utf-8 -*-
"""
Created on Tue Feb 02 21:01:47 2016

@author: Archit Datar
"""
import numpy, scipy
import matplotlib.pyplot as plt
import pandas as pd
import Error_estimation as est

def f (Xdata,pguess):
    Y= pguess[0]*Xdata**pguess[1]
    return Y

location= 'C:\Users\Omkar Mehta\Desktop\Error Estimation\SynthData.xlsx'
df= pd.read_excel(location,0)

Xdata= df['Re']
Xdata=numpy.array(Xdata)
Ydata= df['f']
Ydata= numpy.array(Ydata)
Errdata= df['Error']
Errdata= numpy.array(Errdata)
pguess= [0.0,0.0]
pguess= numpy.array(pguess)

fig=plt.figure()
ax=fig.add_subplot(111)
fig.show()

fig2= plt.figure()
ax2=fig2.add_subplot(111)
fig2.show()


popt, pcov,perr,p95,p_p,chisquare,resanal= est.fitdata(f,Xdata,Ydata,Errdata,pguess,ax=ax,ax2=ax2)

values= numpy.array([popt, pcov,perr,p95,p_p,chisquare,resanal])
names= numpy.array(['popt', 'pcov','perr','p95','p_p','chisquare','resanal'])

df1= pd.DataFrame(values,names)

#writer= pd.ExcelWriter('writing.xlsx')
#df1.to_excel(writer)
print df1