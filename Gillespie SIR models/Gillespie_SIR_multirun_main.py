#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 15:59:54 2021

@author: vandervegt
"""

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
b = 0.3
# mu: probability of recovery
m = 1/10
# lambda: probability of losing immunity
l = 2/365
# p: probability of moving
p = 0.01


# algorithm parameters
n = 10 # number of nodes
c = 3
final_timepoint = 10*365 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes
number_of_runs = 20 # number of times we want to run the stochastic model
p_edge = 0.75

store_totaltimetraces = []
store_fractionInodes = []

for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
    t,x,em = single_model_run(final_timepoint,total_pop,n,p_edge,np.array([b,m,l,p]))
    print("Finished simulation, processing output...")
    # store necessary information from this run
    # total S, I and R over time
    store_total = np.zeros([c,len(t)])
    for i in [0,1,2]:
        temp = x[:,i:-1:3]
        store_total[i] = np.sum(temp,axis = 1)
    store_totaltimetraces.append([t,store_total])
    # number of infected notes over time
    ttemp_I = x[:,1:-1:3] # select all I columns of the output matrix from the simulation
    store_count = np.count_nonzero(ttemp_I,axis = 1)
    store_fractionInodes.append([t,np.divide(store_count,n)])
    print("Finished processing output")
    # other factors for future work: no. of outbreaks in a patch, persistence,etc.

#%%

plt.subplots(figsize = (14,6))
# plot all time traces
plt.subplot(1,2,1)
for i in range(number_of_runs):
    ttemp = store_totaltimetraces[i]
    plt.plot(ttemp[0],ttemp[1][0],color='darkred',lw=1)
    plt.plot(ttemp[0],ttemp[1][1],color='darkgreen',lw=1)
    plt.plot(ttemp[0],ttemp[1][2],color='darkblue',lw=1)
plt.legend(["Total S","Total I","Total R"],fontsize=12)
plt.tick_params(direction='in',size=6) 
plt.title("Compartments summed over all patches")

# plot number of infected nodes over time
plt.subplot(1,2,2)
for i in range(number_of_runs):
    ttemp = store_fractionInodes[i]
    plt.plot(ttemp[0],ttemp[1],color='grey',lw=1)
plt.tick_params(direction='in',size=6) 
plt.title("Number of infected compartments over time")

plotname = "SIRmodel_"+str(n)+"nodes_"+str(total_pop)+"agents.eps"
plt.savefig(plotname)