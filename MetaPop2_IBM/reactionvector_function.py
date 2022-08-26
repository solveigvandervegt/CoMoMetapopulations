import numpy as np
def build_reaction_vector_list(n,c):
    """
    Function builds the list of vectors for the reactions, i.e. infections etc.,
    in the model. This function considers the reactions of an SEIR models. 
    Each vector notes from which node/compartment
    an individual leaves and where it ends up, so when a reaction
    is selected, the appropriate vector can be added to the population vector
    to update.
    
    n: number of nodes
    c: number of compartments
    
    return: numpy array with population update vectors.
    """
    n_rxn = 4 # number of reaction within compartment
    rxn = np.zeros((n*n_rxn, n*c))
    for i in range(n):
        # compartment reactions
        StoE = np.repeat(0,n*c)
        StoE[i*c] = -1
        StoE[i*c+1] = 1
        rxn[i*n_rxn] = StoE
        EtoI = np.repeat(0,n*c)
        EtoI[i*c+1] = -1
        EtoI[i*c+2] = 1
        rxn[i*n_rxn+1] = EtoI
        ItoR = np.repeat(0,n*c)
        ItoR[i*c+2] = -1
        ItoR[i*c+3] = 1
        rxn[i*n_rxn+2] = ItoR
        RtoS = np.repeat(0,n*c)
        RtoS[i*c+3] = -1
        RtoS[i*c] = 1
        rxn[i*n_rxn+3] = RtoS

                
    return rxn