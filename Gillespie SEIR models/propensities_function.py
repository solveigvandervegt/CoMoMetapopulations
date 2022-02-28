import numpy as np
import math
def compute_propensities(t,n,c,params,x,adjacency_matrix,node_degrees,current_rxn,current_prop):
    """
    Function to compute the propensities of the reaction and diffusion processes
    at a given time step.
    
    t: time point
    n: number of nodes
    c: number of compartments
    params: model parameters, vector of the form [b,g,m,l,p]
        b: probability of infection
        g: probability of becoming infectious when exposed
        m: probability of recovering when infectious
        l: probability of losing immunity
        p: probability of moving to a neighbouring node
    x: current population vector
    adjacency_matrix: matrix noting which nodes are connected by edges
    node_degrees: vector indicating the degree of each node
    current_rxn: change-in-state vector of previous reaction
    current_prop: vector of propensities to be update if t>0
    
    returns: vector of propensities
    """
    
    # basic parameter
    # beta: probability of infection
    b = params[0]
    # gamma: probability of becoming infectious
    g = params[1]
    # mu: probability of recovery
    m = params[2]
    # lambda: probability of losing immunity
    l = params[3]
    # p: probability of moving
    p = params[4]
    
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
                    if adjacency_matrix[i][j] == 1:
                        first_index = j*c
                        pmove[first_index:first_index+c] = p/node_degrees[i] *np.array([x[i*c], x[i*c+1], x[i*c+2], x[i*c+3]])
                pmove = np.array(pmove)
            p_temp = np.concatenate((pinfect, pmove))
            # total propensities: concatenate lists
            prop = np.concatenate((prop,p_temp))
    else:
        # we only update those propensities which have changed after the last reaction
        # get number of populations that have been changed
        prop = current_prop.copy()
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
                    if adjacency_matrix[i][j] == 1:
                        first_index = j*c
                        pmove[first_index:first_index+c] = p/node_degrees[i] *np.array([x[i*c], x[i*c+1], x[i*c+2], x[i*c+3]])
                pmove = np.array(pmove)
            p_temp = np.concatenate((pinfect, pmove))
            prop[i*(n+1)*c:(i+1)*(n+1)*c] = p_temp
    return prop