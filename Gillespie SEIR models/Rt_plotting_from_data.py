# plot Rt outputs over time

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 2021

@author: vandervegt
"""
import numpy as np
import matplotlib.pylab as plt
import pickle

#%%
# Import data
# ensure working directory is Gillespie SEIR models folder
with open("./simulations/network_sim_1000_1825_10000_global_Rts.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_global_Rts.bin","rb") as input_file:
    e = pickle.load(input_file)
    number_of_runs = len(e)-1
    global_Rts = e[0]

with open("./simulations/network_sim_1000_1825_10000_local_Rts.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_local_Rts.bin","rb") as input_file:
    e = pickle.load(input_file)
    local_Rts = e[0]
    
with open("./simulations/network_sim_1000_1825_10000_regional_Rts.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_regional_Rts.bin","rb") as input_file:
    e = pickle.load(input_file)
    regional_Rts = e[0]

with open("./simulations/network_sim_1000_1825_10000_time_series.bin","rb") as input_file:
#with open("./simulations/network_sim_100_10_1000_time_series.bin","rb") as input_file:
    e = pickle.load(input_file) #this outputs the time points at which the above data is collected
    time_series = e[0]

# number of nodes
n = 10

color = plt.cm.rainbow(np.linspace(0,1,number_of_runs+1))
#%%

plt.figure(figsize = (9,3))
if number_of_runs == 0:
    plt.plot(time_series,global_Rts,color = color[number_of_runs],linewidth=0.5)
else:
    for ii in range(number_of_runs):
        plt.plot(time_series[ii],global_Rts[ii],color = color[ii],linewidth=0.5)
plt.title('Global Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_testlargesim_global.eps"
plt.savefig(plotname)

plt.figure(figsize = (9,3))
if number_of_runs == 0:
    for jj in range(n):
        plt.plot(time_series,local_Rts[jj],color = color[number_of_runs],linewidth=0.5)
else:
    for ii in range(number_of_runs):
        for jj in range(n):
            plt.plot(time_series[ii],local_Rts[ii][jj],color = color[ii],linewidth=0.5)
plt.title('Local Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_testlargesim_local.eps"
plt.savefig(plotname)

plt.figure(figsize = (9,3))
if number_of_runs == 0:
    for jj in range(n):
        plt.plot(time_series,regional_Rts[jj],color = color[number_of_runs],linewidth=0.5)
else:
    for ii in range(number_of_runs):
        for jj in range(n):
            plt.plot(time_series[ii],regional_Rts[ii][jj],color = color[ii],linewidth=0.5)
plt.title('Regional Rt')
plt.xlabel('time')
plt.ylabel('Rt')
plt.show

plotname = "SEIRmodel_Rt_testlargesim_regional.eps"
plt.savefig(plotname)

#%%
plt.figure(figsize = (9,9))
if number_of_runs == 0:
    for jj in range(n):
        plt.plot(global_Rts,local_Rts[jj],color = color[number_of_runs],linewidth=0.5)
else:
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
#plt.savefig(plotname)