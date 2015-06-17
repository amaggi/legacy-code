#!/usr/bin/env python
#
# Infinite fault (strike-slip, interseismic)
#
# Invert for locking depth and slip rate - Burn-in
#


# Externals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time

# Personals
from Metropolis import *


# -------------------------------------------
def verify(M,prior_bounds):
    '''
    Verify Model
    '''

    # Check stuff
    for i in range(M.size):
        if M[i]<prior_bounds[i,0] or M[i]>prior_bounds[i,1]:
            return False

    # All done
    return True


# -------------------------------------------
def calcLLK(M,data_dict):
    '''
    Compute Log likelihood
    Args:
    * M: Model
    * data_dict: Input data
    '''

    # Parse input
    x    = data_dict['x']
    data = data_dict['data']
    sig  = data_dict['sigma']

    # Normalization factor
    Norm = -(np.log(sig)+0.5*np.log(2*np.pi))*data.size

    # Residual vector
    res  = (data - calcPred(M,x))/sig # residuals

    # Log-Likelihood
    logLLK = Norm-0.5*(res*res).sum() # log of model likelihood

    # All done
    return logLLK


# -------------------------------------------
def calcPred(M,x1):
    '''
    Compute Green's functions
    Args:
    * M: Model (M[0] = slip, M[1] = locking depth)
    * x1: receiver locations
    '''
    GM = M[0]/np.pi*np.arctan(x1/M[1])

    assert not np.isnan(GM).any(), 'NaN value in Green''s functions'

    # All done
    return GM


# -------------------------------------------
def plot_result(M,mtarget,show_ini=True,title = None):
    '''
    Plot results
    Args:
    * M: sample set
    * show_ini: Show initial sample (black star)
    * mtarget
    '''
    xbounds = (0.55,1.1)
    ybounds = (0.65,1.35)    

    ybins = np.linspace(ybounds[0],ybounds[1],250)
    xbins = np.linspace(xbounds[0],xbounds[1],250)

    ### Slip and slip depth histograms
    H, x, y = np.histogram2d(M[:,0],M[:,1],bins=(xbins,ybins),normed = True)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    Hmask = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

    # Figure and axis objects
    figM = plt.figure(facecolor='w')

    gs1 = gridspec.GridSpec(2, 2,width_ratios=[4,1],height_ratios=[4,1])
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])
    ax3 = plt.subplot(gs1[2])

    # 2D Hist
    im = ax1.pcolormesh(x,y,Hmask,zorder=0,vmin=0)#,vmax=0.75*Cm)
    if show_ini:
        ax1.plot(M[0,0],M[0,1],'k*',markersize=10,zorder=2)
    ax1.plot(mtarget[0],mtarget[1],'r*',markersize=10,zorder=2)
    ax1.set_xlabel('Slip rate')
    ax1.xaxis.set_label_coords(0.5, -0.06)
    ax1.set_ylabel('Locking depth')    
    ax1.set_ylim(ybounds[0],ybounds[1])
    ax1.set_xlim(xbounds[0],xbounds[1])
    ax1.invert_yaxis()

    ## Slip depth marginal distribution
    bins = np.linspace(ybounds[0],ybounds[1],50)
    ax2.xaxis.set_ticklabels('')
    ax2.yaxis.set_ticklabels('')
    ax2.hist(M[:,1],bins=bins,range=(ybounds[0],ybounds[1]),normed=True,orientation='horizontal')
    ax2.plot([ax2.get_xlim()[0],ax2.get_xlim()[1]],[mtarget[1],mtarget[1]],'r',lw=4)
    ax2.set_ylim(ybounds[0],ybounds[1])
    ax2.invert_yaxis()        

    ## Slip marginal distribution
    bins = np.linspace(xbounds[0],xbounds[1],50)
    ax3.xaxis.set_ticklabels('')
    ax3.yaxis.set_ticklabels('')
    ax3.hist(M[:,0],bins=bins,range=(xbounds[0],xbounds[1]),normed=True)
    ax3.plot([mtarget[0],mtarget[0]],[ax3.get_ylim()[0],ax3.get_ylim()[1]],'r',lw=4)
    ax3.set_xlim(xbounds[0],xbounds[1])
    ax3.invert_yaxis()

    # Title
    if title is not None:
        plt.suptitle(title)    



### Define input parameters ###
##############################


## Number of samples 
n_samples = 50000


## Target model (Actual "true" solution) 

# Slip rate
slip     = 1.

# Locking depth
lockdepth = 1.0

# Data error 
sigdata = 0.03

# Target model vector
mtarget = np.array([slip,lockdepth])


## Prior Bounds (assuming uniform prior) 
prior_bounds = np.array([[0.0  ,4.0],  # Slip
                         [0.1  ,4.0]]) # Slip depth


## First sample in the markov chain 
m_ini = np.array((0.57,1.32))


## Covariance of the proposal PDF
prop_sigma = np.array([0.035,0.035])
prop_cov   = np.diag(prop_sigma*prop_sigma)


## Calculate synthetic data

# Receiver locations
x = np.linspace(-8*lockdepth,8*lockdepth,100)/lockdepth 
i = np.where(x==0)[0]
if i.size>0:
    x = np.delete(x,i)
    
# Creation of noise free synthetic data
pred = calcPred(mtarget,x)       # Unnormalized noise free data

# Noisy data
data = pred + sigdata*np.random.randn(pred.size); # Noisy data
data_dict = {'x':x,'data':data,'sigma':sigdata}

## Run Metropolis

start_time = time.time()
M,LLK,accepted = Metropolis(n_samples,calcLLK,verify,data_dict,m_ini,prior_bounds,prop_cov)
run_time = time.time() - start_time

# Burn in
iburn = np.where(LLK>=LLK[n_samples/2:].mean())[0][0]

# Mean/STD
M_mean = M[iburn:,:].mean(axis=0)
M_std  = M[iburn:,:].std(axis=0)


## Output display & figures

# Print information
print("--- %s seconds ---" % (run_time))
print("Acceptance rate : %f" %(float(accepted)/n_samples))
print("Posterior mean : %f , %f" %(M_mean[0] ,M_mean[1]))
print("2-sigma error  : %f , %f" %(2*M_std[0],2*M_std[1] ))

# Plot data
pred = calcPred(M_mean[:2],x)
plt.figure(facecolor='w')
ax=plt.subplot(111)
ax.plot(x,data)
ax.plot(x,pred,'r-',lw=2)
ax.legend(('Data','Predictions'),loc='lower right')
ax.set_xlabel('Distance / Locking depth')
ax.set_ylabel('Velocity / Slip rate')

# Burn-in period
plt.figure(facecolor='w')
ax=plt.subplot(111)
xgen = np.arange(n_samples)
ax.semilogx(xgen[:iburn+1],LLK[:iburn+1],'r-')
ax.semilogx(xgen[iburn:],LLK[iburn:],'b-')
ax.legend(('Burn-in period','Post burn-in Samples'),loc='lower right')
ax.set_xlabel('Generations')
ax.set_ylabel('Log(Posterior PDF)')

# Plot results without burn-in
plot_result(M,mtarget,title='Locked fault model - No Burn in')

# Plot results with burn-in
plot_result(M[iburn:,:],mtarget,False,'Locked fault model - Burn in: %d samples'%(iburn))

plt.show()
