import numpy as np
import random
from propensities_function import compute_propensities
def single_model_run_SEIR_eff(t_final,n,b,g,m,l,adjacency_matrix,rxn,group_membership,init_infec):
    """
    Function which runs a single instance of the stochastic SEIR model.
    t_final: the final time point of the simulation.
    total_pop: approximate total number of individuals over all nodes.
    n: number of nodes.
    model parameters, being
        b: probability of infection
        g: probability of becoming infectious
        m: probability of recovery
        l: probability of loss of immunity
    adjacency_matrix: matrix defining the network graph. If there is an edge between nodes
    i and j, then adjacency_matrix[i][j] = adjacency_matrix[j][i] = 1. All other entries 0.
    rxn: list of reaction vectors which can be added to current population to affect
    the changes in population sizes associated with each reaction or diffusion process.
    group_membership: age class for each individual.
    init_infec: number of individuals infected at t=0.
    
    returns: a list of timesteps, and a 2D numpy array of population size per node,
    per compartment, and the edge matrix.
    """
    
    c = 4 # number of compartments
    x = np.repeat(0,n*c) # initialized compartments for all nodes
    x[0:-1:c] = 1 # initialize all individuals as susceptible
    # infect random node
    for ii in range(init_infec):
        init_node = random.randint(0,n-1)*c # pick one individual to be infectious
        x[init_node:init_node+c] = [x[init_node]-1, 0, 1, 0]
    x_store = [x]
    t_store = [0]
    t = 0
    
    # draw random numbers now
    rand_to_generate = 10**6
    random_numbers = np.random.uniform(0,1,rand_to_generate)
    random_number_count = 0
    
    current_rxn = []
    last_prop = []
            
    # start loop over time
    tracker = 0
    while t <= t_final:
        if tracker > 50:
            print(t)
            tracker = 0
        # compute propensities
        prop = compute_propensities(t,n,c,b,g,m,l,x,adjacency_matrix,current_rxn,last_prop,group_membership)
        prop_cumsummed = np.cumsum(prop)
        last_prop = prop.copy()
        
        # check if enough random numbers left
        if random_number_count + 1 >= rand_to_generate:
            random_numbers = np.random.uniform(0,1,rand_to_generate)
            random_number_count = 0
        
        # compute time step
        dt = - np.log(random_numbers[random_number_count])*(1/prop_cumsummed[-1])
        random_number_count += 1
        
        t = t + dt
        tracker += dt
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
        
    return t_store, x_store #np.array(x_store)