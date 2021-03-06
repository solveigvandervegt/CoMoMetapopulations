#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Gillespie_SEIR_function_effupdate import single_model_run_SEIR_eff
from network_function import build_powerlaw_network
from reactionvector_function import build_reaction_vector_list
#import random
#import networkx as nx
import datetime

begin = datetime.datetime.now()

# model parameters
# beta: probability of infection
b = 0.8
# gamma: probability of becoming infectious
g = 1
# mu: probability of recovery
m = 1/5
# lambda: probability of losing immunity
l = 2/365
# p: probability of moving
p = 0.01


# algorithm parameters
n = 20 # number of nodes
c = 4 # number of compartments
final_timepoint = 1*365 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes
number_of_runs = 20 # number of times we want to run the stochastic model

store_totaltimetraces = []
store_fractionInodes = []        
#%%
for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
    # build network
    print('Building network...')
    adjacency_matrix = build_powerlaw_network(n)
    # list of reactions
    print('Building list of reactions...')
    rxn = build_reaction_vector_list(n,c,adjacency_matrix)
    print('Starting simulation...')
    #t,x,em = single_model_run_SEIR(final_timepoint,total_pop,n,np.array([b,g,m,l,p]),adjacency_matrix)
    t,x = single_model_run_SEIR_eff(final_timepoint,total_pop,n,np.array([b,g,m,l,p]),adjacency_matrix,rxn)
    print("Finished simulation, processing output...")
    # store necessary information from this run
    # total S, E, I and R over time
    store_total = np.zeros([c,len(t)])
    for i in [0,1,2,3]:
        temp = x[:,i:-1:c]
        store_total[i] = np.sum(temp,axis = 1)
    store_totaltimetraces.append([t,store_total])
    # number of infected notes over time
    ttemp_I = x[:,2:-1:c] # select all I columns of the output matrix from the simulation
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
    plt.plot(ttemp[0],ttemp[1][1],color='gold',lw=1)
    plt.plot(ttemp[0],ttemp[1][2],color='darkgreen',lw=1)
    plt.plot(ttemp[0],ttemp[1][3],color='darkblue',lw=1)
plt.legend(["Total S","Total E","Total I","Total R"],fontsize=12)
plt.tick_params(direction='in',size=6) 
plt.title("Compartments summed over all patches")

# plot number of infected nodes over time
plt.subplot(1,2,2)
for i in range(number_of_runs):
    ttemp = store_fractionInodes[i]
    plt.plot(ttemp[0],ttemp[1],color='grey',lw=1)
plt.tick_params(direction='in',size=6) 
plt.title("Number of infected compartments over time")

plotname = "SEIRmodel_"+str(n)+"nodes_"+str(total_pop)+"agents.eps"
plt.savefig(plotname)

print(datetime.datetime.now()-begin)