# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 09:41:51 2016

@author: rijulgosar
"""

import numpy as np
import matplotlib.pyplot as plt

Z= np.random.randint(0,3,(6,6))
print Z

#itr=0
bare=0
grass=0
prey=0
#pred=3

fig=plt.figure()
ax=fig.add_subplot(111)

for i in range (1,5):
    for j in range(i,5):
   
#for bare earth       
      if Z[i-1,j-1]==0:
          bare=bare+1
      if Z[i-1,j]==0:
          bare=bare+1
      if Z[i-1,j+1]==0:
          bare=bare+1
      if Z[i,j-1]==0:
          bare=bare+1
      if Z[i,j+1]==0:
          bare=bare+1
      if Z[i+1,j-1]==0:
          bare=bare+1
      if Z[i+1,j]==0:
          bare=bare+1
      if Z[i+1,j+1]==0:
          bare=bare+1
       
      #grass
      if Z[i-1,j-1]==1:
          grass=grass+1
      if Z[i-1,j]==1:
          grass=grass+1
      if Z[i-1,j+1]==1:
          grass=grass+1
      if Z[i,j-1]==1:
          grass=grass+1
      if Z[i,j+1]==1:
          grass=grass+1
      if Z[i+1,j-1]==1:
          grass=grass+1
      if Z[i+1,j]==1:
          grass=grass+1
      if Z[i+1,j+1]==1:
          grass=grass+1
       
      #prey 
      if Z[i-1,j-1]==2:
          prey=prey+1
      if Z[i-1,j]==2:
          prey=prey+1
      if Z[i-1,j+1]==2:
          prey=prey+1
      if Z[i,j-1]==2:
          prey=prey+1
      if Z[i,j+1]==2:
          prey=prey+1
      if Z[i+1,j-1]==2:
          prey=prey+1
      if Z[i+1,j]==2:
          prey=prey+1
      if Z[i+1,j+1]==2:
          prey=prey+1
      
      
      if Z[i,j]==0 and grass>0:
          Z[i,j]=1
    
      if Z[i,j]==0 and prey>=1 and grass>=2:
          Z[i,j]=2
          
      if Z[i,j]==2 and grass==1:
          Z[i,j]=0
    
      if Z[i,j]==1 and prey>1:
          Z[i,j]==0
          
    #print Z
    
    ax.imshow(Z,interpolation='nearest',cmap=plt.cm.RdBu_r)
    plt.pause(0.2)
    plt.show()
    
    