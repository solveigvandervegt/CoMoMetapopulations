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
#from Rt_function import Rt
from local_Rt_function import compute_local_Rt
from regional_Rt_function import compute_regional_Rt
from global_Rt_function import compute_global_Rt
#import random
#import networkx as nx
#import datetime
import math
import pickle

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
final_timepoint = 2 * 365 # final time point for simulations
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

    # process outputs of simulation to compute Rt values per simulation
    # select earliest timepoint in each day
    t_day = [t[0]]
    x_day = [x[0]]
    for ii in range(1,len(t)):
        if math.floor(t[ii]) != math.floor(t[ii-1]):
            t_day.append(t[ii])
            x_day.append(x[ii])


    #store all Rt values appropriately
    if i ==0:
        time_series = [t_day]
        population_series = [x_day]
    else:
        time_series.append(t_day)
        population_series.append(x_day)

basefile = "simulations/network_sim_" + str(n) + "_" + str(final_timepoint) + "_" + str(total_pop) + "_"
dump_file(population_series, "population_series", basefile)
dump_file(time_series, "time_series", basefile)
