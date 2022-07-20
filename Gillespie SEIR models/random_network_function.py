import numpy as np
def build_random_network(n,p):
    """
    Function which builds a single-component network where two nodes are connected
    with a given probability p.
    
    n: number of nodes
    p: probability of an edge between two given nodes
    
    returns: adjecency matrix of the network
    """
    # create nxn matrix of random numbers between 0 and 1
    edge_matrix = np.zeros([n,n])
    start = np.random.randint(0,n)
    while np.sum(np.sum(edge_matrix, axis = 1) >= 1) <  n:
        go_to_node = np.random.randint(0,n)
        while start == go_to_node:
            go_to_node = np.random.randint(0,n)
        edge_matrix[start][go_to_node] = 1
        edge_matrix[go_to_node][start] = 1
        start = go_to_node
    
    return edge_matrix