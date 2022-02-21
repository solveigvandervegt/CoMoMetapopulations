import numpy as np
import random
import math
def single_model_run_SEIR_eff(t_final,total_pop,n,par,edges_mat,rxn):
    """
    Function which runs a single instance of the stochastic SEIR model.
    t_final: the final time point of the simulation.
    total_pop: approximate total number of individuals over all nodes.
    n: number of nodes.
    par: model parameters, being
        b: probability of infection
        g: probability of becoming infectious
        m: probability of recovery
        l: probability of loss of immunity
        p: probability of movement between nodes
    edges_mat: matrix defining the network graph. If there is an edge between nodes
    i and j, then edges_mat[i][j] = edges_mat[j][i] = 1. All other entries 0.
    rxn: list of reaction vectors which can be added to current population to affect
    the changes in population sizes associated with each reaction or diffusion process.
    
    returns: a list of timesteps, and a 2D numpy array of population size per node,
    per compartment, and the edge matrix.
    """
    # basic parameter
    # beta: probability of infection
    b = par[0]
    # gamma: probability of becoming infectious
    g = par[1]
    # mu: probability of recovery
    m = par[2]
    # lambda: probability of losing immunity
    l = par[3]
    # p: probability of moving
    p = par[4]
    
    
    # try for n nodes, fully connected
    # initial conditions of the form S1, I1, R1, S2, I2, R2, etc.
    c = 4 # number of compartments
    approx_total_number_individuals = total_pop
    avg_pop_size = approx_total_number_individuals/n
    x = np.repeat(0,n*c) # initialized compartments for all nodes
    x[0:-1:c] = np.rint(np.random.normal(1,0.5,n)*avg_pop_size) # initial all susceptible populations
    # infect random node
    init_node = random.randint(0,n-1)*c
    x[init_node:init_node+c] = [x[init_node]*0.99, 0, x[init_node]*0.01, 0] # initial infections in first node
    x[x < 0] = 0
    x_store = [x]
    t_store = [0]
    t = 0
    
    # draw random numbers now
    rand_to_generate = 10**6
    random_numbers = np.random.uniform(0,1,rand_to_generate)
    random_number_count = 0
    
    
    n_edges_per_node = [np.sum(x) for x in edges_mat]
            
        
    # start loop over time
    tracker = 0
    while t <= t_final:
        #if tracker > 50:
            #print(t)
            #tracker = 0
        # compute propensities if first time point
        if t == 0:
            prop = np.array([])
            for i in range(n):
                # patch 1: S to E, E to I, I to R, R to S, move 1 to 2 (for each comp), move 1 to 3 (for each comp)
                N = np.sum(x[i*c:(i+1)*c])
                if N == 0:
                    pinfect = np.zeros(c)
                    pmove = np.zeros(n*c)
                else:
                    pinfect = np.array([(1-(1-b/N)**x[i*c+2])*x[i*c], g*x[i*c+1], m*x[i*c+2], l*x[i*c+3]])
                    pmove = np.zeros(n*c)
                    for j in range(n):
                        if edges_mat[i][j] == 1:
                            first_index = j*c
                            pmove[first_index:first_index+c] = p/n_edges_per_node[i] *np.array([x[i*c], x[i*c+1], x[i*c+2], x[i*c+3]])
                    pmove = np.array(pmove)
                p_temp = np.concatenate((pinfect, pmove))
                # total propensities: concatenate lists
                prop = np.concatenate((prop,p_temp))
        else:
            # we only update those propensities which have changed after the last reaction
            # get number of populations that have been changed
            pos_change = np.where(current_rxn == 1)
            neg_change = np.where(current_rxn == -1)
            for i in [math.floor(pos_change[0]/4), math.floor(neg_change[0]/4)]:
                N = np.sum(x[i*c:(i+1)*c])
                if N == 0:
                    pinfect = np.zeros(c)
                    pmove = np.zeros(n*c)
                else:
                    pinfect = np.array([(1-(1-b/N)**x[i*c+2])*x[i*c], g*x[i*c+1], m*x[i*c+2], l*x[i*c+3]])
                    pmove = np.zeros(n*c)
                    for j in range(n):
                        if edges_mat[i][j] == 1:
                            first_index = j*c
                            pmove[first_index:first_index+c] = p/n_edges_per_node[i] *np.array([x[i*c], x[i*c+1], x[i*c+2], x[i*c+3]])
                    pmove = np.array(pmove)
                p_temp = np.concatenate((pinfect, pmove))
                prop[i*(n+1)*c:(i+1)*(n+1)*c] = p_temp
        
        prop_cumsummed = np.cumsum(prop)
        
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
        
    return t_store, np.array(x_store)