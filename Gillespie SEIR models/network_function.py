import numpy as np
import random
def build_powerlaw_network(n):
    """
    Function which builds a single-component network where the degrees of the nodes
    are distributed according to the power law P(k) = k^-3.
    
    n: number of nodes
    
    returns: adjecency matrix of the network
    """
    # new degree generation
    min_degrees = 2
    max_degrees = np.sqrt(n).astype(int)
    range_degrees = [x for x in range(min_degrees,max_degrees+1)]
    distribution_probabilities = [x**(-3) for x in range_degrees]
    degree_probabilities = [x/sum(distribution_probabilities) for x in distribution_probabilities]
    list_of_degrees = np.random.choice(range_degrees,size = n,p = degree_probabilities)
    
    # check if sum of degrees is even
    if np.sum(list_of_degrees)%2 != 0:
        r_int = random.randint(0,n-1)
        list_of_degrees[r_int] += 1
    
    # create list of stubs
    list_of_stubs =[]
    for i in range(n):
        list_of_stubs.extend(list_of_degrees[i]*[i])
    
    #build adjacency matrix
    list_of_stubs_store = list_of_stubs.copy()
    adjacency_matrix = np.zeros((n,n)) # initiate edges matrix
    list_of_edges = [] # initiate list of all edges in network
    while len(list_of_stubs)>0:
        if (len(list_of_stubs) == 2) & ((list_of_stubs[0] == list_of_stubs[1]) or (list_of_stubs in list_of_edges)): # cannot connect to own node
            # if fault in last edge, just start over
            list_of_stubs = list_of_stubs_store.copy()
            adjacency_matrix = np.zeros((n,n)) # initiate edges matrix
            list_of_edges = [] # initiate list of all edges in network
        else :
            edge = random.sample(list_of_stubs,2) # create edge from two stubs
            if (edge[0] != edge[1]) & ~(edge in list_of_edges): # check if not connecting to self and edge doesn't already exist
                adjacency_matrix[edge[0]][edge[1]] = 1 # connect nodes in edges matrix
                adjacency_matrix[edge[1]][edge[0]] = 1
                list_of_stubs.remove(edge[0]) # remove stubs from list
                list_of_stubs.remove(edge[1])
                list_of_edges.append(edge) 
                list_of_edges.append([edge[1], edge[0]])
    
    return adjacency_matrix