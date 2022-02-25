import numpy as np
def build_reaction_vector_list(n,c,adjacency_matrix):
    """
    Function builds the list of vectors for the reaction and diffusion actions
    in the model. This function considers the reactions of an SEIR model, with
    movement between connected nodes. Each vector notes from which node/compartment
    an individual leaves and where it ends up, so when a reaction or diffusion
    is selected, the appropriate vector can be added to the population vector
    to update.
    
    n: number of nodes
    c: number of compartments
    adjacency_matrix: matrix noting which nodes are connected.
    
    return: numpy array with population update vectors.
    """
    n_rxn = 4 # number of reaction within compartment
    rxn = np.zeros((n*(n_rxn+n*c), n*c))
    for i in range(n):
        # compartment reactions
        StoE = np.repeat(0,n*c)
        StoE[i*c] = -1
        StoE[i*c+1] = 1
        rxn[i*(n_rxn+(n*c))] = StoE
        EtoI = np.repeat(0,n*c)
        EtoI[i*c+1] = -1
        EtoI[i*c+2] = 1
        rxn[i*(n_rxn+(n*c))+1] = EtoI
        ItoR = np.repeat(0,n*c)
        ItoR[i*c+2] = -1
        ItoR[i*c+3] = 1
        rxn[i*(n_rxn+(n*c))+2] = ItoR
        RtoS = np.repeat(0,n*c)
        RtoS[i*c+3] = -1
        RtoS[i*c] = 1
        rxn[i*(n_rxn+(n*c))+3] = RtoS
        # movement reactions
        #count = 0
        for j in range(n):
            if adjacency_matrix[i][j] == 1:
                Sitoj = np.repeat(0,n*c)
                Sitoj[i*c] = -1
                Sitoj[j*c] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c] = Sitoj
                Eitoj = np.repeat(0,n*c)
                Eitoj[i*c + 1] = -1
                Eitoj[j*c + 1] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 1] = Eitoj
                Iitoj = np.repeat(0,n*c)
                Iitoj[i*c + 2] = -1
                Iitoj[j*c + 2] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 2] = Iitoj
                Ritoj = np.repeat(0,n*c)
                Ritoj[i*c + 3] = -1
                Ritoj[j*c + 3] = 1
                rxn[i*(n_rxn+(n*c)) + c + j*c + 3] = Ritoj
                
    return rxn