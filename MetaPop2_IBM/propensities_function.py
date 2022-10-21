import numpy as np
import math
def compute_propensities(t,n,c,b,g,m,l,x,adjacency_matrix,current_rxn,current_prop,group_membership):
    """
    Function to compute the propensities of the reaction and diffusion processes
    at a given time step.
    
    t: time point
    n: number of nodes
    c: number of compartments
    model parameters
        b: probability of infection
        g: probability of becoming infectious when exposed
        m: probability of recovering when infectious
        l: probability of losing immunity
    x: current population vector
    adjacency_matrix: matrix noting which nodes are connected by edges
    node_degrees: vector indicating the degree of each node
    current_rxn: change-in-state vector of previous reaction
    current_prop: vector of propensities to be update if t>0
    group_membership: which age class does the individual belong to
    
    returns: vector of propensities
    """
    
    if t == 0:
        prop = np.array([])
        for i in range(n): #loop through each individual
            # for each path: S to E, E to I, I to R, R to S
            # compute the number of people that an individual is in touch with that are infected
            I_neighbors = 0
            for j in range(n): #loop through to find all neighbors
                if (adjacency_matrix[i,j] == 1 and x[j*c+2] == 1):
                    I_neighbors += 1
            pinfect = np.array([(1-(1-b[group_membership[i]])**I_neighbors)*x[i*c], g[group_membership[i]]*x[i*c+1], m[group_membership[i]]*x[i*c+2], l*x[i*c+3]])
            # total propensities: concatenate lists
            prop = np.concatenate((prop,pinfect))
    else:
        # we only update those propensities which have changed after the last reaction
        # get number of populations that have been changed
        prop = current_prop.copy()
        pos_change = np.where(current_rxn == 1)
        i = math.floor(pos_change[0]/4) # this is the individual we're changing the propensities for
        I_neighbors = 0
        for j in range(n): #looping through to find neighbors
            if (adjacency_matrix[i,j] == 1 and x[j*c+2] == 1):
                I_neighbors += 1
        pinfect = np.array([(1-(1-b[group_membership[i]])**I_neighbors)*x[i*c], g[group_membership[i]]*x[i*c+1], m[group_membership[i]]*x[i*c+2], l*x[i*c+3]])
        prop[i*c:(i+1)*c] = pinfect
    return prop