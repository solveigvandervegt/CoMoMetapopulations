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
import pickle

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
n = 1000 # number of nodes
c = 4 # number of compartments
final_timepoint = 5 * 365 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes
number_of_runs = 1 # number of times we want to run the stochastic model


def dump_file(obj, obj_name, basefile):
    with open(basefile + obj_name + ".bin", "wb") as f:
        pickle.dump(obj, f)

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
        local_Rts.append(store_local_Rt)
        regional_Rts.append(store_regional_Rt)
        global_Rts.append(store_global_Rt)
        time_series.append(t_day)

    basefile = "simulations/network_sim_" + str(n) + "_" + str(final_timepoint) + "_" + str(total_pop) + "_"
    dump_file(local_Rts, "local_Rts", basefile)
    dump_file(regional_Rts, "regional_Rts", basefile)
    dump_file(global_Rts, "global_Rts", basefile)
    dump_file(time_series, "time_series", basefile)



pickle_in = open(basefile + "global_Rts.bin","rb")
test = pickle.load(pickle_in)
print(test)
