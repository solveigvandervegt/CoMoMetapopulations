#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 14:19:06 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Tau_SEIR_function_diffnetwork import single_model_run_SEIR
import networkx as nx

# model parameters
# beta: probability of infection
b = 0.15
# gamma: probability of becoming infectious
g = 1/3
# mu: probability of recovery
m = 1/10
# lambda: probability of losing immunity
l = 2/365
# p: probability of moving
p = 0.1


# algorithm parameters
n = 5 # number of nodes
c = 4
final_timepoint = 5 # final time point for simulations
total_pop = 100 # approximate total population over all nodes
p_edge = 0.75
tau = 0.1

t,x,em = single_model_run_SEIR(final_timepoint,total_pop,n,p_edge,np.array([b,g,m,l,p]),tau)

# plot results overall
store_total = np.zeros([c,len(t)])
for i in [0,1,2,3]:
    temp = x[:,i:-1:c]
    store_total[i] = np.sum(temp,axis = 1)

#%%
plt.plot(t,store_total[0],color='darkred',lw=1)
plt.plot(t,store_total[1],color='gold',lw=1)
plt.plot(t,store_total[2],color='darkgreen',lw=1)
plt.plot(t,store_total[3],color='darkblue',lw=1)
plt.legend(["Total S","Total E","Total I","Total R"],fontsize=12)
plt.tick_params(direction='in',size=6) 
plt.title("Compartments summed over all patches")
#%%
# plot results per node
plt.subplots(figsize = (14,15))
for i in range(1,10+1):
    plt.subplot(5,2,i)
    i = i
    plt.plot(t,[item[(i-1)*c] for item in x],color='darkred',lw=2)
    plt.plot(t,[item[(i-1)*c+1] for item in x],color='gold',lw=2)
    plt.plot(t,[item[(i-1)*c+2] for item in x],color='darkgreen',lw=2)
    plt.plot(t,[item[(i-1)*c+3] for item in x],color='darkblue',lw=2)
    strS = "Patch %d: S" % i
    strE = "Patch %d: E" % i
    strI = "Patch %d: I" % i
    strR = "Patch %d: R" % i
    plt.legend([strS,strE,strI,strR],fontsize=12)
    plt.tick_params(direction='in',size=6) 
#%%
# plot network graph
G = nx.from_numpy_array(em, create_using=nx.DiGraph)
labels = {}
for i in range(n):
    labels[i] = i+1
# G.degree 
plt.figure(figsize=(9, 6))
n_color = np.asarray([G.degree[n] for n in G.nodes])
nx.draw(G, labels = labels, node_color = n_color)