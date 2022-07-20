def compute_global_Rt(t,x,n,b,m,c):
    """ 
    compute Rt value across the network for one run of the model.
    
    t: time series with one entry per day
    x: 3D array with, for each run, the populations on each day
    n: number of nodes
    b: probability of infection
    m: probability of recovery
    c: number of compartments
    
    return: Rt value at the each node at each time point
    """
    from Rt_function import Rt
    import numpy as np
    
    global_total_pop = sum(x[0])
    store_global_Rt = [0]*len(t) #np.zeros(len(t_day))
    for jj in range(len(t)):
        temp = x[jj][0:-1:c]
        St = np.sum(temp)
        store_global_Rt[jj] = Rt(St/global_total_pop,b,m)
    
    return store_global_Rt
