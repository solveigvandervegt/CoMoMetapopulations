import numpy as np
import random
import math
def single_model_run_SEIR(t_final,total_pop,n,p_edge,par,tau):
    """
    Function which runs a single instance of the stochastic SEIR model.
    t_final: the final time point of the simulation.
    total_pop: approximate total number of individuals over all nodes.
    n: number of nodes.
    p_edge: probability of an edge existing between two nodes
    par: model parameters, being
        b: probability of infection
        g: probability of becoming infectious
        m: probability of recovery
        l: probability of loss of immunity
        p: probability of movement between nodes
    tau: time step for tau-leap algorithm
    
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
    x[0:c] = [x[0]*0.95, 0, x[0]*0.05, 0] # initial infections in first node
    x[x < 0] = 0
    x_store = [x]
    t_store = [0]
    t = 0
    
    # draw random numbers now
    rand_to_generate = 10**6
    random_numbers = np.random.uniform(0,1,rand_to_generate)
    random_number_count = 0
    
    # matrix of egdes between nodes
    min_edges = 2 # minimum number of edges for a node
    max_edges = np.floor(np.sqrt(n)) # maximum number of edges for a node
    max_edges = max_edges.astype(int)
    pdf = [x**(-3) for x in range(min_edges,max_edges+1)]
    cdf = np.cumsum(pdf)
    n_edges_per_node = np.random.uniform(0,cdf[-1],n)
    list_of_stubs = []
    for i in range(n):
        c_edges = 0
        # find minimum number of degrees for which U(0,1) is more than poisson prob.
        while n_edges_per_node[i] > cdf[c_edges]:
            c_edges += 1
        n_edges_per_node[i] = c_edges + min_edges # save number of edges
        list_of_stubs.extend([i]*(c_edges+min_edges)) # save number of stubs for network building
    if len(list_of_stubs)%2 != 0: #if the number of edges is not even, we need to fix that
        r_int = random.randint(0,n-1)
        while list_of_stubs.count(list_of_stubs[r_int]) <= min_edges: # check that we don't decrease degree belwo minimum
            r_int = random.randint(0,n-1)
        list_of_stubs.remove(r_int)
    edges_mat = np.zeros((n,n)) # initiate edges matrix
    list_of_edges = [] # initiate list of all edges in network
    while len(list_of_stubs)>0:
        if (len(list_of_stubs) == 2) & (list_of_stubs[0] == list_of_stubs[1]): # cannot connect to own node
            # break up previously made edge, and make two new ones
            edges_mat[last_edge[0]][last_edge[1]] = 0
            edges_mat[last_edge[1]][last_edge[0]] = 0
            edges_mat[last_edge[0]][list_of_stubs[0]] = 1
            edges_mat[last_edge[1]][list_of_stubs[1]] = 1
            edges_mat[list_of_stubs[0]][last_edge[0]]= 1
            edges_mat[list_of_stubs[1]][last_edge[1]] = 1
            break
        edge = random.sample(list_of_stubs,2) # create edge from two stubs
        if (edge[0] != edge[1]) & ~(edge in list_of_edges): # check if not connecting to self and edge doesn't already exist
            edges_mat[edge[0]][edge[1]] = 1 # connect nodes in edges matrix
            edges_mat[edge[1]][edge[0]] = 1
            list_of_stubs.remove(edge[0]) # remove stubs from list
            list_of_stubs.remove(edge[1])
            last_edge = edge
            list_of_edges.append(edge) 

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
        x_store.append(x)
        
    return t_store, np.array(x_store), edges_mat