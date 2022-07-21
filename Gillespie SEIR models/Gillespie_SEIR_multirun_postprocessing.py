#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 2021

@author: vandervegt
"""

import numpy as np
import matplotlib.pylab as plt
from local_Rt_function import compute_local_Rt
from regional_Rt_function import compute_regional_Rt
from global_Rt_function import compute_global_Rt
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
final_timepoint = 730 # final time point for simulations
total_pop = 10**4 # approximate total population over all nodes
number_of_runs = 1 # number of times we want to run the stochastic model


#%%

with open("./simulations/network_sim_100_730_10000_time_series.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_global_Rts.bin","rb") as input_file:
    time_series = pickle.load(input_file)

with open("./simulations/network_sim_100_730_10000_population_series.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_global_Rts.bin","rb") as input_file:
    population_series = pickle.load(input_file)

#%%
# NOTE: functions compute_X_Rt compute Rts for ONE run of the mode.
# in case of multiple runs, need to edit the below
if number_of_runs > 1:
    print('All code written below this statement is for one run of the model. You must adapt to allow for processing of multiple runs.')
print('Starting computation of local Rt values...')
local_Rts = compute_local_Rt(time_series[0],population_series[0],n,b,m)
#print('Starting computation of regional Rt values...')
#regional_Rts = compute_regional_Rt(time_series[0],population_series[0],n,b,m,adjacency_matrix)
print('Starting computation of global Rt values...')
global_Rts = compute_global_Rt(time_series[0],population_series[0],n,b,m,c)
#%%
print('Start plotting...')
# plot Rt outputs over time
color = plt.cm.rainbow(np.linspace(0,1,number_of_runs))
plt.figure(figsize = (9,3))
plt.plot(time_series[0],global_Rts,color = color[0],linewidth=0.5)
plt.title('Global Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

#%%
plotname = "SEIRmodel_Rt_test_global.eps"
plt.savefig(plotname)

plt.figure(figsize = (9,3))
for jj in range(n):
    plt.plot(time_series[0],local_Rts[jj],color = color[0],linewidth=0.5)
plt.title('Local Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_test_local.eps"
plt.savefig(plotname)

#%%
plt.figure(figsize = (9,3))
for jj in range(n):
    plt.plot(time_series[0],regional_Rts[jj],color = color[0],linewidth=0.5)
plt.title('Regional Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_test_regional.eps"
plt.savefig(plotname)

#%%
plt.figure(figsize = (9,9))
for jj in range(n):
    plt.plot(global_Rts,local_Rts[jj],color = color[0],linewidth=0.5)
plt.plot([x for x in range(5)],[x for x in range(5)],'k--',linewidth=0.5)
plt.title('Global vs Local Rt')
plt.xlabel('Global Rt')
plt.ylabel('Local Rt')
plt.xlim([0,4.1])
plt.ylim([0,4.1])
plt.show

plotname = "SEIRmodel_Rt_test.eps"
plt.savefig(plotname)

#%%
S_temp = np.zeros((100,731))
for ii in range(len(time_series[0])):
    for jj in range(n):
        S_temp[jj][ii] = population_series[0][ii][jj*4]

S_temp_sum = np.sum(S_temp,axis=0)
#%%
plt.figure()
for ii in range(n):
    plt.plot(time_series[0],S_temp[ii])

plt.figure()
plt.plot(time_series[0],S_temp_sum)