#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Gillespie_SEIR_function_diffnetwork import single_model_run_SEIR
import random
import networkx as nx

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
total_pop = 10**5 # approximate total population over all nodes
number_of_runs = 20 # number of times we want to run the stochastic model

store_totaltimetraces = []
store_fractionInodes = []

# build network
# matrix of egdes between nodes
min_edges = 2 # minimum number of edges for a node
max_edges = np.floor(np.sqrt(n)) # maximum number of edges for a node
max_edges = max_edges.astype(int)
# get poisson distribution of degrees
pdf = [x**(-3) for x in range(min_edges,max_edges+1)]
cdf = np.cumsum(pdf)
n_edges_per_node = np.random.uniform(0,cdf[-1],n)
list_of_stubs = []
# get number of edges per node, save and save in list of stubs
for i in range(n):
    c_edges = 0
    # find minimum number of degrees for which U(0,1) is more than poisson prob.
    while n_edges_per_node[i] > cdf[c_edges]:
        c_edges += 1
    n_edges_per_node[i] = c_edges + min_edges # save number of edges
    list_of_stubs.extend([i]*n_edges_per_node[i].astype(int)) # save number of stubs for network building
if len(list_of_stubs)%2 != 0: #if the number of edges is not even, we need to fix that
    r_int = random.randint(0,n-1)
    while list_of_stubs.count(list_of_stubs[r_int]) <= min_edges: # check that we don't decrease degree belwo minimum
        r_int = random.randint(0,n-1)
    list_of_stubs.remove(r_int)
edges_mat = np.zeros((n,n)) # initiate edges matrix
list_of_edges = [] # initiate list of all edges in network
print('Building network...')
while len(list_of_stubs)>0:
    if (len(list_of_stubs) == 2) & (list_of_stubs[0] == list_of_stubs[1]): # cannot connect to own node
        # break up previously made edge, and make two new ones
        edges_mat[last_edge[0]][last_edge[1]] = 0
        edges_mat[last_edge[1]][last_edge[0]] = 0
        edges_mat[last_edge[0]][list_of_stubs[0]] = 1
        edges_mat[last_edge[1]][list_of_stubs[1]] = 1
        edges_mat[list_of_stubs[0]][last_edge[0]]= 1
        edges_mat[list_of_stubs[1]][last_edge[1]] = 1
        break
    edge = random.sample(list_of_stubs,2) # create edge from two stubs
    if (edge[0] != edge[1]) & ~(edge in list_of_edges) & ~([edge[1], edge[0]] in list_of_edges): # check if not connecting to self and edge doesn't already exist
        edges_mat[edge[0]][edge[1]] = 1 # connect nodes in edges matrix
        edges_mat[edge[1]][edge[0]] = 1
        list_of_stubs.remove(edge[0]) # remove stubs from list
        list_of_stubs.remove(edge[1])
        last_edge = edge
        list_of_edges.append(edge) 

#%%
for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
    t,x,em = single_model_run_SEIR(final_timepoint,total_pop,n,np.array([b,g,m,l,p]),edges_mat)
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