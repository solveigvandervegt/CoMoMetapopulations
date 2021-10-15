#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 14:19:06 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Gillespie_SIR_function import single_model_run

# model parameters
# beta: probability of infection
b = 0.4
# mu: probability of recovery
m = 1/2.5
# lambda: probability of losing immunity
l = 0.005
# p: probability of moving
p = 0.01


# algorithm parameters
n = 10 # number of nodes
c = 3
final_timepoint = 3*365 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes

t,x,em = single_model_run(final_timepoint,total_pop,n,np.array([b,m,l,p]))

# plot results overall
store_total = np.zeros([c,len(t)])
for i in [0,1,2]:
    temp = x[:,i:-1:3]
    store_total[i] = np.sum(temp,axis = 1)

plt.plot(t,store_total[0],color='darkred',lw=2)
plt.plot(t,store_total[1],color='darkgreen',lw=2)
plt.plot(t,store_total[2],color='darkblue',lw=2)
plt.legend(["Total S","Total I","Total R"],fontsize=12)
plt.tick_params(direction='in',size=6) 
plt.title("Compartments summed over all patches")

# plot results per node
plt.subplots(figsize = (14,15))
for i in range(1,n+1):
    plt.subplot(5,2,i)
    plt.plot(t,[item[(i-1)*c] for item in x],color='darkred',lw=2)
    plt.plot(t,[item[(i-1)*c+1] for item in x],color='darkgreen',lw=2)
    plt.plot(t,[item[(i-1)*c+2] for item in x],color='darkblue',lw=2)
    strS = "Patch %d: S" % i
    strI = "Patch %d: I" % i
    strR = "Patch %d: R" % i
    plt.legend([strS,strI,strR],fontsize=12)
    plt.tick_params(direction='in',size=6) 

