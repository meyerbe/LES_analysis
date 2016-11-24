import numpy as np
import sys
from parameters import *
from scipy.integrate import odeint

from thermo_aux import latent_heat


def main():
    p, T = ClausiusClapeyron(200,400)
    return

# ------------------------------------------------
''' Integrate the ClausiusClapeyron equation '''
#def initialize(self, namelist, LatentHeat):
#def ClausiusClapeyron():
def ClausiusClapeyron(Tmin, Tmax):
#    try:
#        Tmin = namelist['ClausiusClapeyron']['temperature_min']
#    except:
#        print('Clasius-Clayperon lookup table temperature_min not '
#                       'given in name list taking default of 100.15 K')
#        Tmin = 100.15
#    
#    try:
#        Tmax = namelist['ClausiusClapeyron']['temperature_max']
#    except:
#        print('Clasius-Clayperon lookup table temperature_max not '
#                       'given in name list taking default of 380 K')
#        Tmax = 380.0
#    
#    try:
#        n_lookup = namelist['ClausiusClapeyron']['n_lookup']
#    except:
#        print('Clasius-Clayperon lookup table n_lookup not '
#                       'given in name list taking default of 128')
#        n_lookup = 512

    n_lookup = 512
    
    #Generate array of equally space temperatures
    T = np.linspace(Tmin, Tmax, n_lookup)
    
    #Find the maximum index of T where T < T_tilde
    tp_close_index = np.max(np.where(T<=Tt))
    
    #Check to make sure that T_tilde is not in T
    if T[tp_close_index] == Tt:
        print('Array of temperatures for ClasiusClapyeron lookup table contains Tt  \n')
        print('Pick different values for ClasiusClapyeron Tmin and Tmax in lookup table \n')
        print('Killing Simulation now!')
        sys.exit()
    else:
        print('Tt not in T-array')
    
    #Now prepare integration
    T_above_Tt= np.append([Tt],T[tp_close_index+1:])
    T_below_Tt= np.append(T[:tp_close_index+1],[Tt])[::-1]

    #set the initial condition
    pv0 = np.log(pv_star_t)

    #Integrate
    pv_above_Tt = np.exp(odeint(rhs,pv0,T_above_Tt,hmax=0.1)[1:])
    pv_below_Tt = np.exp(odeint(rhs,pv0,T_below_Tt,hmax=0.1)[1:])[::-1]
    pv = np.append(pv_below_Tt,pv_above_Tt )
#    self.LT.initialize(T,pv)

    print(pv.shape)
    print(T.shape)


    return pv, T


# Set up the RHS
def rhs(z,T_):
    # lam = LH.Lambda(T_)
    lam = 1.0
    # L = LH.L(T_,lam)
    L = latent_heat(T_)
    return L/(Rv * T_ * T_)


# ------------------------------------------------
if __name__ == '__main__':
    main()
