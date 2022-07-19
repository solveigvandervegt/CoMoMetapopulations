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
import math

#begin = datetime.datetime.now()

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
n = 100 # number of nodes
c = 4 # number of compartments
final_timepoint = 2*365 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes
number_of_runs = 1 # number of times we want to run the stochastic model


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
    # select earliest timepoint in each day
    t_day = [t[0]]
    x_day = [x[0]]
    for ii in range(1,len(t)):
        if math.floor(t[ii]) != math.floor(t[ii-1]):
            t_day.append(t[ii])
            x_day.append(x[ii])

    # local
    store_local_Rt = [[0]*len(t_day) for _ in range(n)]#np.zeros((n,len(t_day)))
    for ii in range(n):
        for jj in range(len(t_day)):
            store_local_Rt[ii][jj] = Rt(x_day[jj][ii*4]/sum(x_day[jj][ii*4:ii*4+4]),b,m)

    #regional: at most one edge from center node
    store_regional_Rt = [[0]*len(t_day) for _ in range(n)]#np.zeros((n,len(t_day)))
    for ii in range(n):
        for jj in range(len(t_day)):
            St = x_day[jj][ii*4]
            regional_total_pop = sum(x_day[jj][ii*4:ii*4+4])
            for kk in range(len(adjacency_matrix[ii])):
                if adjacency_matrix[ii][kk] == 1:
                    St = St + x_day[jj][int(kk)*4]
                    regional_total_pop = regional_total_pop + sum(x_day[jj][int(kk)*4:int(kk)*4+4])
            store_regional_Rt[ii][jj] = Rt(St/regional_total_pop,b,m)

    #global
    global_total_pop = sum(x[0])
    store_global_Rt = [0]*len(t_day) #np.zeros(len(t_day))
    for jj in range(len(t_day)):
        temp = x_day[jj][0:-1:c]
        St = np.sum(temp)
        store_global_Rt[jj] = Rt(St/global_total_pop,b,m)

    #store all Rt values appropriately
    if i ==0:
        local_Rts = [store_local_Rt]
        regional_Rts = [store_regional_Rt]
        global_Rts = [store_global_Rt]
        time_series = [t_day]
    else:
        #np.append(local_Rts,[store_local_Rt],axis=0)
        #np.append(regional_Rts,[store_regional_Rt],axis=0)
        #np.append(global_Rts,[store_global_Rt],axis=0)
        local_Rts.append(store_local_Rt)
        regional_Rts.append(store_regional_Rt)
        global_Rts.append(store_global_Rt)
        time_series.append(t_day)

#%%
# plot Rt outputs over time
color = plt.cm.rainbow(np.linspace(0,1,number_of_runs))
plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    plt.plot(time_series[ii],global_Rts[ii],color = color[ii],linewidth=0.5)
plt.title('Global Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_test_global.eps"
plt.savefig(plotname)

plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    for jj in range(n):
        plt.plot(time_series[ii],local_Rts[ii][jj],color = color[ii],linewidth=0.5)
plt.title('Local Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_test_local.eps"
plt.savefig(plotname)

plt.figure(figsize = (9,3))
for ii in range(number_of_runs):
    for jj in range(n):
        plt.plot(time_series[ii],regional_Rts[ii][jj],color = color[ii],linewidth=0.5)
plt.title('Regional Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_test_regional.eps"
plt.savefig(plotname)

#%%
plt.figure(figsize = (9,9))
for ii in range(number_of_runs):
    for jj in range(n):
        plt.plot(global_Rts[ii],local_Rts[ii][jj],color = color[ii],linewidth=0.5)
plt.plot([x for x in range(5)],[x for x in range(5)],'k--',linewidth=0.5)
plt.title('Global vs Local Rt')
plt.xlabel('Global Rt')
plt.ylabel('Local Rt')
plt.xlim([0,4])
plt.ylim([0,4])
plt.show

plotname = "SEIRmodel_Rt_test.eps"
plt.savefig(plotname)
