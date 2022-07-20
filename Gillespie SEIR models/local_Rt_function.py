def compute_local_Rt(t,x,n,b,m):
    """ 
    compute Rt value for each node in the network.
    
    t: time series with one entry per day
    x: 3D array with, for each run, the populations on each day
    n: number of nodes
    b: probability of infection
    m: probability of recovery
    
    return: Rt value at the each node at each time point
    """
    from Rt_function import Rt
        
    store_local_Rt = [[0]*len(t) for _ in range(n)]#np.zeros((n,len(t_day)))
    for ii in range(n):
        for jj in range(len(t)):
            store_local_Rt[ii][jj] = Rt(x[jj][ii*4]/sum(x[jj][ii*4:ii*4+4]),b,m)
    
    return store_local_Rt
