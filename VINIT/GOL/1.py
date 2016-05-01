# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 16:38:51 2016

@author: Vinit Dedhiya
"""
import matplotlib.pyplot as plt
import numpy

def iterate(Z):
    # find number of neighbours that each square has
    N = numpy.zeros(Z.shape)
    N[1:, 1:] += Z[:-1, :-1]
    N[1:, :-1] += Z[:-1, 1:]
    N[:-1, 1:] += Z[1:, :-1]
    N[:-1, :-1] += Z[1:, 1:]
    N[:-1, :] += Z[1:, :]
    N[1:, :] += Z[:-1, :]
    N[:, :-1] += Z[:, 1:]
    N[:, 1:] += Z[:, :-1]
    # a live cell is killed if it has fewer 
    # than 2 or more than 3 neighbours.
    part1 = ((Z == 1) & (N < 4) & (N > 1)) 
    # a new cell forms if a square has exactly three members
    part2 = ((Z == 0) & (N == 3))
    return (part1 | part2).astype(int)

Z = numpy.array([[0,0,0,0,0,0],
                 [0,0,0,1,0,0],
                 [0,1,0,1,0,0],
                 [0,0,1,1,0,0],
                 [0,0,0,0,0,0],
                 [0,0,0,0,0,0]])
print 'Initial state:'
print Z[1:-1,1:-1]
for i in range(4):
    Z = iterate(Z)
print 'Final state:'
print Z[1:-1,1:-1]
imgplot=plt.imshow(Z)
plt.pause(0.1)

'''
def iterate(Z):
    shape = len(Z), len(Z[0])
    N  = [[0,]*(shape[0]+2)  for i in range(shape[1]+2)]
    # Compute number of neighbours for each cell
    for x in range(1,shape[0]-1):
        for y in range(1,shape[1]-1):
            N[x][y] = Z[x-1][y-1]+Z[x][y-1]+Z[x+1][y-1] \
                    + Z[x-1][y]            +Z[x+1][y]   \
                    + Z[x-1][y+1]+Z[x][y+1]+Z[x+1][y+1]
    # Update cells
    for x in range(1,shape[0]-1):
        for y in range(1,shape[1]-1):
            if Z[x][y] == 0 and N[x][y] == 3:
                Z[x][y] = 1
            elif Z[x][y] == 1 and not N[x][y] in [2,3]:
                Z[x][y] = 0
    return Z

def display(Z):
    shape = len(Z), len(Z[0])
    for x in range(1,shape[0]-1):
        for y in range(1,shape[1]-1):
            print Z[x][y],
        print
    print

Z = [[0,0,0,0,0,0],
     [0,0,0,1,0,0],
     [0,1,0,1,0,0],
     [0,0,1,1,0,0],
     [0,0,0,0,0,0],
     [0,0,0,0,0,0]]
print 'Initial state:'
display(Z)
for i in range(4):
    iterate(Z)
print 'Final state:'
display(Z)
'''
