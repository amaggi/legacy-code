'''
Simple metropolis sampler class
'''

import numpy as np
np.random.seed(0)

def Metropolis(n_samples,fun_calcLLK,fun_verify,data,m_ini,prior_bounds,prop_cov):

    ''' 
    Metropolis algorithm 
    Args:
        * n_samples: Number of samples
        * fun_calcLLK: Function to calculate Log-likelihood
        * fun_verify: Function to verify model
        * data: Data
        * m_ini: Initial model
        * prior_bounds: Parameter bounds (uniform priors)
        * prop_cov: Covariance of the (gaussian) proposal PDF
    '''

    # Initiate model chain
    M   = np.zeros((n_samples,m_ini.size))
    LLK = np.zeros((n_samples,))
    M[0,:] = m_ini # Initialize the chain
    LLK[0] = fun_calcLLK(m_ini,data) 

    # Proposal mean
    prop_mean = np.zeros(m_ini.shape)

    # Iteration loop
    count = 0 # Number of accepted models
    for i in range(1,n_samples): 
    
        # Random walk
        Mnew = M[i-1,:].copy()
        Mnew += np.random.multivariate_normal(prop_mean,prop_cov)

        # Check if model sample is within prior
        valid = fun_verify(Mnew,prior_bounds)
        if not valid:
            M[i,:] = M[i-1,:].copy()
            LLK[i] = LLK[i-1]
            continue

        # LLK of new model        
        logLLKn = fun_calcLLK(Mnew, data)
        dLLK    = logLLKn - LLK[i-1]
        
        # Metropolis acceptance/rejection
        accept = 0
        U = np.log(np.random.rand())
        if U < dLLK:
            accept=1        
        else:
            accept=0

        # Update model
        if accept:
            M[i,:]  = Mnew.copy()
            LLK[i]  = logLLKn
            count += 1
        else:
            M[i,:] = M[i-1,:].copy()
            LLK[i] = LLK[i-1]

    return M,LLK,count
