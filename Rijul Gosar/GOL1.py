# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 08:59:37 2016

@author: Eeshani Godbole
"""
import matplotlib.pyplot as plt
import numpy
z=numpy.random.randint(2,size=(10,10))
"""z=[[0,0,0,0,0,0],
   [0,0,0,1,1,0],
   [0,0,1,1,0,0],
   [0,0,0,0,1,0],
   [0,1,0,1,0,0],
   [0,0,0,0,0,0]]"""
   
def fun(z):
    size=len(z),len(z[0])
    for x in range(1,size[0]-1):
        for y in range(1,size[1]-1):
            n=[z[x-1][y-1],z[x][y-1],z[x+1][y-1],z[x-1][y],z[x+1][y],z[x-1][y+1],z[x][y+1],z[x+1][y+1]]
            
            if z[x][y]==0 and n.count(0)==3:
                z[x][y]=1
            if z[x][y]==1 and n.count(1)>3:
                z[x][y]=0
            if z[x][y]==1 and n.count(1)<2:
                z[x][y]=0
            if z[x][y]==1 and n.count(1)==2:
                z[x][y]=1
            if z[x][y]==1 and n.count(1)==3:
                z[x][y]=1
        plt.imshow(z,interpolation="none")
        plt.show()
        plt.pause(0.00000000001)
        print(z)

    return z    
fun(z)