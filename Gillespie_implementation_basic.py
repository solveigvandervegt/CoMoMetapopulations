#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 14:19:06 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
import networkx as nx

# basic parameter
# beta: probability of infection
b = 0.15
# mu: probability of recovery
m = 0.07
# lambda: probability of losing immunity
l = 0.01
# p: probability of moving
p = 0.05


# try for n nodes, fully connected
# initial conditions of the form S1, I1, R1, S2, I2, R2, etc.
n = 10 # number of nodes
c = 3 # number of compartments
approx_total_number_individuals = 5000 
avg_pop_size = approx_total_number_individuals/n
x = np.repeat(0,n*c) # initialized compartments for all nodes
x[0:-1:3] = np.rint(np.random.normal(1,0.5,10)*avg_pop_size) # initial all susceptible populations
x[0:3] = [x[0]*0.95, x[0]*0.05, 0] # initial infections in first node
x[x < 0] = 0
x_store = [x]
t_store = [0]
t = 0
t_final = 365*3

# draw random numbers now
rand_to_generate = 10**6
random_numbers = np.random.uniform(0,1,rand_to_generate)
random_number_count = 0

# matrix of egdes between nodes
random_mat = np.random.random((n,n))
edges_mat = np.zeros((n,n))
for i in range(n):
    for j in range(i-1):
        if (random_mat[i][j] <= 0.4 and i != j):
            edges_mat[i][j] = 1
            edges_mat[j][i] = 1 # assuming the connections go both directions


# list of reactions
rxn = np.zeros((n*(3+(n-1)*c), n*c))
for i in range(n):
    # compartment reactions
    StoI = np.repeat(0,n*c)
    StoI[i*c] = -1
    StoI[i*c+1] = 1
    rxn[i*(3+((n-1)*c))] = StoI
    ItoR = np.repeat(0,n*c)
    ItoR[i*c+1] = -1
    ItoR[i*c+2] = 1
    rxn[i*(3+((n-1)*c))+1] = ItoR
    RtoS = np.repeat(0,n*c)
    RtoS[i*c+2] = -1
    RtoS[i*c] = 1
    rxn[i*(3+((n-1)*c))+2] = RtoS
    # movement reactions
    count = 0
    for j in range(n):
        if edges_mat[i][j] == 1:
            Sitoj = np.repeat(0,n*c)
            Sitoj[i*c] = -1
            Sitoj[j*c] = 1
            rxn[i*(3+((n-1)*c)) + c + count*c] = Sitoj
            Iitoj = np.repeat(0,n*c)
            Iitoj[i*c + 1] = -1
            Iitoj[j*c + 1] = 1
            rxn[i*(3+((n-1)*c)) + c + count*c + 1] = Iitoj
            Ritoj = np.repeat(0,n*c)
            Ritoj[i*c+2] = -1
            Ritoj[j*c+2] = 1
            rxn[i*(3+((n-1)*c)) + c + count*c + 2] = Ritoj
            count += 1
    


# start loop over time
while t <= t_final:
    # compute propensities
    prop = np.array([])
    for i in range(n):
        # patch 1: S to I, I to R, R to S, move 1 to 2 (for each comp), move 1 to 3 (for each comp)
        N = np.sum(x[i*c:(i+1)*c])
        if N == 0:
            pinfect = np.zeros(c)
            pmove = np.zeros((n-1)*c)
        else:
            pinfect = np.array([(1-(1-b/N)**x[i*c+1])*x[i*c], m*x[i*c+1], l*x[i*c+2]])
            pmove = (n-1)*[(p/2)*(x[i*c]/N), (p/2)*(x[i*c+1]/N), (p/2)*(x[i*c+2]/N)]
            pmove = np.array(pmove)
        p_temp = np.concatenate((pinfect, pmove))
        # total propensities: concatenate lists
        prop = np.concatenate((prop,p_temp))
        
    prop_cumsummed = np.cumsum(prop)
    
    if random_number_count + 1 >= rand_to_generate:
        random_numbers = np.random.uniform(0,1,rand_to_generate)
        random_number_count = 0
    
    # compute time step
    dt = - np.log(random_numbers[random_number_count])*(1/prop_cumsummed[-1])
    random_number_count += 1
    
    t = t + dt
    t_store.append(t)
    
    # choose reaction
    rn = random_numbers[random_number_count]
    random_number_count += 1
    
    prop_cumsummed = np.divide(prop_cumsummed,prop_cumsummed[-1])
    
    for count in range(len(prop)):
        if rn <= prop_cumsummed[count]:
            break
        
    # find reaction in list from above
    current_rxn = rxn[count]
    # add reaction to population vector
    x = [sum(temp) for temp in zip(x, current_rxn)]
    x_store.append(x)
    
# plot results
x_store = np.array(x_store)
# plot results overall
store_total = np.zeros([c,len(t_store)])
for i in [0,1,2]:
    temp = x_store[:,i:-1:3]
    store_total[i] = np.sum(temp,axis = 1)

plt.plot(t_store,store_total[0],color='darkred',lw=2)
plt.plot(t_store,store_total[1],color='darkgreen',lw=2)
plt.plot(t_store,store_total[2],color='darkblue',lw=2)
plt.legend(["Total S","Total I","Total R"],fontsize=12)
plt.tick_params(direction='in',size=6) 
plt.title("Compartments summed over all patches")

# plot results per node
plt.subplots(figsize = (14,15))
for i in range(1,n+1):
    plt.subplot(5,2,i)
    plt.plot(t_store,[item[(i-1)*c] for item in x_store],color='darkred',lw=2)
    plt.plot(t_store,[item[(i-1)*c+1] for item in x_store],color='darkgreen',lw=2)
    plt.plot(t_store,[item[(i-1)*c+2] for item in x_store],color='darkblue',lw=2)
    strS = "Patch %d: S" % i
    strI = "Patch %d: I" % i
    strR = "Patch %d: R" % i
    plt.legend([strS,strI,strR],fontsize=12)
    plt.tick_params(direction='in',size=6) 


# plot network graph
G = nx.from_numpy_array(edges_mat, create_using=nx.DiGraph)
labels = {}
for i in range(n):
    labels[i] = i+1
# G.degree 
plt.figure(figsize=(9, 6))
nx.draw(G, labels = labels)