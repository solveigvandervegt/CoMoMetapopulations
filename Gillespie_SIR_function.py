import numpy as np
def single_model_run(t_final,total_pop,nodes,p_edge,par):
    """
    Function which runs a single instance of the stochastic SIR model.
    t_final: the final time point of the simulation.
    total_pop: approximate total number of individuals over all nodes.
    nodes: number of nodes.
    p_edge: probability of an edge existing between two nodes
    par: model parameters, being
        b: probability of infection
        m: probability of recovery
        l: probability of loss of immunity
        p: probability of movement between nodes
    
    returns: a list of timesteps, and a 2D numpy array of population size per node,
    per compartment, and the edge matrix.
    """
    # basic parameter
    # beta: probability of infection
    b = par[0]
    # mu: probability of recovery
    m = par[1]
    # lambda: probability of losing immunity
    l = par[2]
    # p: probability of moving
    p = par[3]
    
    
    # try for n nodes, fully connected
    # initial conditions of the form S1, I1, R1, S2, I2, R2, etc.
    n = nodes # number of nodes
    c = 3 # number of compartments
    approx_total_number_individuals = total_pop
    avg_pop_size = approx_total_number_individuals/n
    x = np.repeat(0,n*c) # initialized compartments for all nodes
    x[0:-1:c] = np.rint(np.random.normal(1,0.5,10)*avg_pop_size) # initial all susceptible populations
    x[0:c] = [x[0]*0.95, x[0]*0.05, 0] # initial infections in first node
    x[x < 0] = 0
    x_store = [x]
    t_store = [0]
    t = 0
    
    # draw random numbers now
    rand_to_generate = 10**6
    random_numbers = np.random.uniform(0,1,rand_to_generate)
    random_number_count = 0
    
    # matrix of egdes between nodes
    random_mat = np.random.random((n,n))
    edges_mat = np.zeros((n,n))
    for i in range(n):
        for j in range(i-1):
            if (random_mat[i][j] <= p_edge and i != j):
                edges_mat[i][j] = 1
                edges_mat[j][i] = 1 # assuming the connections go both directions
    
    
    # list of reactions
    rxn = np.zeros((n*(3+(n-1)*c), n*c))
    for i in range(n):
        # compartment reactions
        StoI = np.repeat(0,n*c)
        StoI[i*c] = -1
        StoI[i*c+1] = 1
        rxn[i*(3+((n-1)*c))] = StoI
        ItoR = np.repeat(0,n*c)
        ItoR[i*c+1] = -1
        ItoR[i*c+2] = 1
        rxn[i*(3+((n-1)*c))+1] = ItoR
        RtoS = np.repeat(0,n*c)
        RtoS[i*c+2] = -1
        RtoS[i*c] = 1
        rxn[i*(3+((n-1)*c))+2] = RtoS
        # movement reactions
        count = 0
        for j in range(n):
            if edges_mat[i][j] == 1:
                Sitoj = np.repeat(0,n*c)
                Sitoj[i*c] = -1
                Sitoj[j*c] = 1
                rxn[i*(3+((n-1)*c)) + c + count*c] = Sitoj
                Iitoj = np.repeat(0,n*c)
                Iitoj[i*c + 1] = -1
                Iitoj[j*c + 1] = 1
                rxn[i*(3+((n-1)*c)) + c + count*c + 1] = Iitoj
                Ritoj = np.repeat(0,n*c)
                Ritoj[i*c+2] = -1
                Ritoj[j*c+2] = 1
                rxn[i*(3+((n-1)*c)) + c + count*c + 2] = Ritoj
                count += 1
        
    
    
    # start loop over time
    while t <= t_final:
        # compute propensities
        prop = np.array([])
        for i in range(n):
            # patch 1: S to I, I to R, R to S, move 1 to 2 (for each comp), move 1 to 3 (for each comp)
            N = np.sum(x[i*c:(i+1)*c])
            if N == 0:
                pinfect = np.zeros(c)
                pmove = np.zeros((n-1)*c)
            else:
                pinfect = np.array([(1-(1-b/N)**x[i*c+1])*x[i*c], m*x[i*c+1], l*x[i*c+2]])
                pmove = (n-1)*[(p/2)*(x[i*c]/N), (p/2)*(x[i*c+1]/N), (p/2)*(x[i*c+2]/N)]
                pmove = np.array(pmove)
            p_temp = np.concatenate((pinfect, pmove))
            # total propensities: concatenate lists
            prop = np.concatenate((prop,p_temp))
            
        prop_cumsummed = np.cumsum(prop)
        
        if random_number_count + 1 >= rand_to_generate:
            random_numbers = np.random.uniform(0,1,rand_to_generate)
            random_number_count = 0
        
        # compute time step
        dt = - np.log(random_numbers[random_number_count])*(1/prop_cumsummed[-1])
        random_number_count += 1
        
        t = t + dt
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
        
    return t_store, np.array(x_store), edges_mat