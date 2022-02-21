import numpy as np
import random
import math
def single_model_run_SEIR(t_final,total_pop,n,par,tau,edges_mat):
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
    tau: time step for tau-leap algorithm
    edges_mat: matrix defining the network
    
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
    n_rxn = 4 # number of reaction within compartment
    approx_total_number_individuals = total_pop
    avg_pop_size = approx_total_number_individuals/n
    x = np.repeat(0,n*c) # initialized compartments for all nodes
    x[0:-1:c] = np.rint(np.random.normal(1,0.5,n)*avg_pop_size) # initial all susceptible populations
    init_node = random.randint(0,n-1)*c
    x[init_node:init_node+c] = [x[init_node]*0.95, 0, x[init_node]*0.05, 0] # initial infections in first node
    print(init_node)
    x[x < 0] = 0
    x_store = [x]
    t_store = [0]
    t = 0
    
    # draw random numbers now
    rand_to_generate = 10**6
    random_numbers = np.random.uniform(0,1,rand_to_generate)
    random_number_count = 0
    
    # comupte degree of each node
    n_edges_per_node = [np.sum(x) for x in edges_mat]

    # list of reactions
    rxn = np.zeros((n*(n_rxn+n*c), n*c))
    for i in range(n):
        # compartment reactions
        StoE = np.repeat(0,n*c)
        StoE[i*c] = -1
        StoE[i*c+1] = 1
        rxn[i*(n_rxn+(n*c))] = StoE
        EtoI = np.repeat(0,n*c)
        EtoI[i*c+1] = -1
        EtoI[i*c+2] = 1
        rxn[i*(n_rxn+(n*c))+1] = EtoI
        ItoR = np.repeat(0,n*c)
        ItoR[i*c+2] = -1
        ItoR[i*c+3] = 1
        rxn[i*(n_rxn+(n*c))+2] = ItoR
        RtoS = np.repeat(0,n*c)
        RtoS[i*c+3] = -1
        RtoS[i*c] = 1
        rxn[i*(n_rxn+(n*c))+3] = RtoS
        # movement reactions
        #count = 0
        for j in range(n):
            if edges_mat[i][j] == 1:
                Sitoj = np.repeat(0,n*c)
                Sitoj[i*c] = -1
                Sitoj[j*c] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c] = Sitoj
                Eitoj = np.repeat(0,n*c)
                Eitoj[i*c + 1] = -1
                Eitoj[j*c + 1] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 1] = Eitoj
                Iitoj = np.repeat(0,n*c)
                Iitoj[i*c + 2] = -1
                Iitoj[j*c + 2] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 2] = Iitoj
                Ritoj = np.repeat(0,n*c)
                Ritoj[i*c + 3] = -1
                Ritoj[j*c + 3] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 3] = Ritoj
                #count += 1
  
    # start loop over time
    while t <= t_final:
        if t%1 == 0:
            print(t)
        # compute propensities
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
            
        prop_cumsummed = np.cumsum(prop)
        prop_cumsummed = np.divide(prop_cumsummed,prop_cumsummed[-1])
        
        if random_number_count + len(prop) >= rand_to_generate:
            random_numbers = np.random.uniform(0,1,rand_to_generate)
            random_number_count = 0
        
        t = t + tau
        t_store.append(t)
        # compute number of occurences of each reaction
        store_events = np.zeros(len(prop))
        for jj in range(len(prop)):
            rn = random_numbers[random_number_count]
            random_number_count += 1
            prob = 0
            sum_p = ((prop_cumsummed[jj]*tau)**prob)/math.factorial(prob)
            while rn > np.exp(-prop_cumsummed[jj]*tau) * sum_p:
                prob += 1
                sum_p = sum_p + ((prop_cumsummed[jj]*tau)**prob)/math.factorial(prob)
            store_number_of_events = prob
            store_events = [sum(temp) for temp in zip(store_events, store_number_of_events*rxn[jj])]
        
        
        # add reaction to population vector
        x = [sum(temp) for temp in zip(x,store_events)]
        x = [max(temp,0) for temp in x]
        x_store.append(x)
        
    return t_store, np.array(x_store), edges_mat