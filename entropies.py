''' Entropy functions from LES: entropies.h '''

from parameters import *
import numpy as np


# dry entropy
def sd_c(pd,T):
    return sd_tilde + cpd*np.log(T/T_tilde) - Rv*np.log(pd/p_tilde)

# moist entropy
def sv_c(pv,T):
    if pv == 0:
        print('in sv_c: pv=0')
        return 0.0
    else:
        return sv_tilde + cpv*np.log(T/T_tilde) - Rd*np.log(pv/p_tilde)

# liquid water entropy
def sc_c(L, T):
    return -L/T



