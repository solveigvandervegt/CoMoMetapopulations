def compute_regional_Rt(t,x,n,b,m,adjacency_matrix,neighbors = 1):
    """ 
    compute Rt value for a neighborhood around each node in the network for one
    run of the model.
    
    t: time series with one entry per day
    x: 3D array with, for each run, the populations on each day
    n: number of nodes
    b: probability of infection
    m: probability of recovery
    adjacency_matrix: nxn matrix denoting which nodes are connected by edges
    
    return: Rt value at the each node at each time point
    """
    from Rt_function import Rt
    store_regional_Rt = [[0]*len(t) for _ in range(n)]#np.zeros((n,len(t_day)))
    for ii in range(n):
        for jj in range(len(t)):
            St = x[jj][ii*4]
            regional_total_pop = sum(x[jj][ii*4:ii*4+4])
            for kk in range(len(adjacency_matrix[ii])):
                if adjacency_matrix[ii][kk] == 1:
                    St = St + x[jj][int(kk)*4]
                    regional_total_pop = regional_total_pop + sum(x[jj][int(kk)*4:int(kk)*4+4])
            store_regional_Rt[ii][jj] = Rt(St/regional_total_pop,b,m)

    
    return store_regional_Rt
