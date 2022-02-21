#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Gillespie_SEIR_function_diffnetwork import single_model_run_SEIR
from Gillespie_SEIR_function_effupdate import single_model_run_SEIR_eff
import random
import networkx as nx
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

store_totaltimetraces = []
store_fractionInodes = []

# build network
print('Building network...')
# new degree generation (for old, see singlerun)
list_of_degrees = nx.generators.random_graphs.random_powerlaw_tree_sequence(n,gamma = 3,tries=5000)
list_of_degrees = [x+1 for x in list_of_degrees] # minimum degree 2 of each node
if np.sum(list_of_degrees)%2 != 0:
    r_int = random.randint(0,n-1)
    list_of_degrees[r_int] += 1

list_of_stubs =[]
for i in range(n):
    list_of_stubs.extend(list_of_degrees[i]*[i])


adjacency_matrix = np.zeros((n,n)) # initiate edges matrix
list_of_edges = [] # initiate list of all edges in network
while len(list_of_stubs)>0:
    if (len(list_of_stubs) == 2) & (list_of_stubs[0] == list_of_stubs[1]): # cannot connect to own node
        # break up previously made edge, and make two new ones
        adjacency_matrix[last_edge[0]][last_edge[1]] = 0
        adjacency_matrix[last_edge[1]][last_edge[0]] = 0
        adjacency_matrix[last_edge[0]][list_of_stubs[0]] = 1
        adjacency_matrix[last_edge[1]][list_of_stubs[1]] = 1
        adjacency_matrix[list_of_stubs[0]][last_edge[0]]= 1
        adjacency_matrix[list_of_stubs[1]][last_edge[1]] = 1
        break
    edge = random.sample(list_of_stubs,2) # create edge from two stubs
    if (edge[0] != edge[1]) & ~(edge in list_of_edges) & ~([edge[1], edge[0]] in list_of_edges): # check if not connecting to self and edge doesn't already exist
        adjacency_matrix[edge[0]][edge[1]] = 1 # connect nodes in edges matrix
        adjacency_matrix[edge[1]][edge[0]] = 1
        list_of_stubs.remove(edge[0]) # remove stubs from list
        list_of_stubs.remove(edge[1])
        last_edge = edge
        list_of_edges.append(edge) 
        

# list of reactions
print('Building list of reactions...')
n_rxn = 4 # number of reaction within compartment
rxn = np.zeros((n*(n_rxn+n*c), n*c))
for i in range(n):
    # compartment reactions
    StoE = np.repeat(0,n*c)
    StoE[i*c] = -1
    StoE[i*c+1] = 1
    rxn[i*(n_rxn+(n*c))] = StoE
    EtoI = np.repeat(0,n*c)
    EtoI[i*c+1] = -1
    EtoI[i*c+2] = 1
    rxn[i*(n_rxn+(n*c))+1] = EtoI
    ItoR = np.repeat(0,n*c)
    ItoR[i*c+2] = -1
    ItoR[i*c+3] = 1
    rxn[i*(n_rxn+(n*c))+2] = ItoR
    RtoS = np.repeat(0,n*c)
    RtoS[i*c+3] = -1
    RtoS[i*c] = 1
    rxn[i*(n_rxn+(n*c))+3] = RtoS
    # movement reactions
    #count = 0
    for j in range(n):
        if adjacency_matrix[i][j] == 1:
            Sitoj = np.repeat(0,n*c)
            Sitoj[i*c] = -1
            Sitoj[j*c] = 1
            rxn[i*(n_rxn+(n*c)) + c + j*c] = Sitoj
            Eitoj = np.repeat(0,n*c)
            Eitoj[i*c + 1] = -1
            Eitoj[j*c + 1] = 1
            rxn[i*(n_rxn+(n*c)) + c + j*c + 1] = Eitoj
            Iitoj = np.repeat(0,n*c)
            Iitoj[i*c + 2] = -1
            Iitoj[j*c + 2] = 1
            rxn[i*(n_rxn+(n*c)) + c + j*c + 2] = Iitoj
            Ritoj = np.repeat(0,n*c)
            Ritoj[i*c + 3] = -1
            Ritoj[j*c + 3] = 1
            rxn[i*(n_rxn+(n*c)) + c + j*c + 3] = Ritoj

#%%
for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
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