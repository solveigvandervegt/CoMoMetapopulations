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
from Rt_function import Rt
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
number_of_runs = 5 # number of times we want to run the stochastic model


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
    #if i == 0:
    #    store_all = x
    #else:
    #    np.append(store_all,[x],axis=0)

    # process outputs of simulation to compute Rt values per simulation
    # local
    store_local_Rt = np.zeros((n,len(t))) 
    for ii in range(n):
        for jj in range(len(t)):
            store_local_Rt[ii][jj] = Rt(x[jj][ii*4],b,m)
    
    #regional: at most one edge from center node
    store_regional_Rt = np.zeros((n,len(t)))
    for ii in range(n):
        for jj in range(len(t)):
            St = x[jj][ii*4]
            for kk in adjacency_matrix[ii]:
                if kk == 1:
                    St = St + x[jj][int(kk)*4]
            store_regional_Rt[ii][jj] = Rt(St,b,m)
    
    #global
    store_global_Rt = np.zeros(len(t)) 
    for jj in range(len(t)):
        temp = x[jj,0:-1:c]
        St = np.sum(temp)
        store_global_Rt[jj] = Rt(St,b,m)
    
    #store all Rt values appropriately
    if i ==0:
        local_Rts = [store_local_Rt]
        regional_Rts = [store_regional_Rt]
        global_Rts = [store_global_Rt]
        time_series = [t]
    else:
        np.append(local_Rts,[store_local_Rt],axis=0)
        np.append(regional_Rts,[store_regional_Rt],axis=0)
        np.append(global_Rts,[store_global_Rt],axis=0)
        np.append(time_series,[t],axis=0)

#%% 
# plot Rt outputs over time
color = plt.cm.rainbow(np.linspace(0,1,number_of_runs))
plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    plt.plot(time_series[ii],global_Rts[ii],color = color[ii])
plt.title('Global Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    for jj in range(n):
        plt.plot(time_series[ii],local_Rts[ii][jj],color = color[ii])
plt.title('Local Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    for jj in range(n):
        plt.plot(time_series[ii],regional_Rts[ii][jj],color = color[ii])
plt.title('Regional Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show
