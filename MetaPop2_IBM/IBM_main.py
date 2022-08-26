#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from Gillespie_SEIR_function_effupdate import single_model_run_SEIR_eff
from reactionvector_function import build_reaction_vector_list
import pickle

# model parameters
# beta: probability of infection
b = 0.8
# gamma: probability of becoming infectious
g = 1
# mu: probability of recovery
m = 1/10
# lambda: probability of losing immunity
l = 2/365



# algorithm parameters
c = 4 # number of compartments
final_timepoint = 100 # final time point for simulations
number_of_runs = 20 # number of times we want to run the stochastic model

def dump_file(obj, obj_name, basefile):
    with open(basefile + obj_name + ".bin", "wb") as f:
        pickle.dump(obj, f)

#%% run simulations with overdispersed network
plt.figure()
for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
    # import network
    print('Import network...')
    # choose an overdispersed network
    with open("./example_network_overdispersed.pickle","rb") as input_file:
        temp = pickle.load(input_file)
    group_membership = temp[1]
    n = len(group_membership) #number of individuals in simulation
    adjacency_matrix = temp[0]
    # list of reactions
    print('Building list of reactions...')
    rxn = build_reaction_vector_list(n,c)
    print('Starting simulation...')
    t,x = single_model_run_SEIR_eff(final_timepoint,n,np.array([b,g,m,l]),adjacency_matrix,rxn)
    print("Finished simulation, plotting I over time...")


    number_of_I = np.zeros(len(t))
    for i in range(len(t)):
        number_of_I[i] = np.sum(x[i][2:-1:c])

    plt.plot(t,number_of_I,color='grey',lw=1)
    # process outputs of simulation to compute Rt values per simulation
    # select earliest timepoint in each day
    #t_day = [t[0]]
    #x_day = [x[0]]
    #for ii in range(1,len(t)):
    #    if math.floor(t[ii]) != math.floor(t[ii-1]):
    #        t_day.append(t[ii])
    #        x_day.append(x[ii])

    #if i ==0:
    #    time_series = [t_day]
    #    population_series = [x_day]
    #else:
    #    time_series.append(t_day)
    #    population_series.append(x_day)

#basefile = "simulations/network_sim_" + str(n) + "_" + str(final_timepoint) + "_" + "_"
#dump_file(population_series, "population_series", basefile)
#dump_file(time_series, "time_series", basefile)
#dump_file(adjacency_matrix, "adjacency_matrix", basefile)
    


plt.tick_params(direction='in',size=6) 
plt.title("Number of infected individuals over time: overdispersed network")

#%% run simulations with poisson network
plt.figure()
for i in range(number_of_runs):
    print("Starting instance %d" %(i+1))
    # import network
    print('Import network...')
    # choose a poisson network
    with open("./example_network_poisson.pickle","rb") as input_file:
        #with open("./simulations/network_sim_100_10_1000_global_Rts.bin","rb") as input_file:
        temp = pickle.load(input_file)
    group_membership = temp[1]
    n = len(group_membership) #number of individuals in simulation
    adjacency_matrix = temp[0]
    # list of reactions
    print('Building list of reactions...')
    rxn = build_reaction_vector_list(n,c)
    print('Starting simulation...')
    t_poisson,x_poisson = single_model_run_SEIR_eff(final_timepoint,n,np.array([b,g,m,l]),adjacency_matrix,rxn)
    print('Final timepoint')
    print(t_poisson[-1])
    print("Finished simulation, plotting I over time...")


    number_of_I = np.zeros(len(t_poisson))
    for i in range(len(t_poisson)):
        number_of_I[i] = np.sum(x_poisson[i][2:-1:c])

    plt.plot(t_poisson,number_of_I,color='grey',lw=1)
    

plt.tick_params(direction='in',size=6) 
plt.title("Number of infected individuals over time: Poisson network")

