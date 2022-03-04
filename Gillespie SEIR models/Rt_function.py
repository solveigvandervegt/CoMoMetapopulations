
def Rt(St,b,g,d=0):
    """
    Function to compute R_t for a given number of susceptibles and parameter
    values.
    
    St: fraction of the population that is susceptible.
    b: beta parameter, rate of infection.
    g: gamma parameter, recovery rate.
    d: death rate, default zero as our model has no deaths
    
    returns: value of time-varying reproduction number, Rt.
    """
    return St*b/(g+d)