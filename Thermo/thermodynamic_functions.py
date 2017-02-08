''' Auxiliary thermodynamic functions from LES '''
import sys

from math import *
import numpy as np

sys.path.append("..")
from parameters import *



# ---------------------------------------------------------------------------
'''Pressure '''
# partial vapour pressure
def pv_c(p0,qt,qv):
    return p0 * eps_vi * qv / (1.0-qt+eps_vi*qv)

# ---------------------------------------------------------------------------
''' Entropy '''
# inline double sd_c(double pd, double T){
#     return sd_tilde + cpd*log(T/T_tilde) -Rd*log(pd/p_tilde);
# }
#
# inline double sv_c(double pv, double T){
#     return sv_tilde + cpv*log(T/T_tilde) - Rv * log(pv/p_tilde);
# }
#
# inline double sc_c(double L, double T){
#     return -L/T;
# }
def s_dry(pd, T):
    return sd_tilde

def s_vap(pv, T):
    return sv_tilde

def s_cond(L, T):
    return -L/T
# ---------------------------------------------------------------------------
''' Temperature '''
# T(pd,pv,s,qt)
def temperature_no_ql(pd,pv,s,qt):
    cp = ((1.0 - qt) * cpd + qt * cpv)
    if pv > 0:
        temp = T_tilde * exp( ( s - (1.0-qt)*(sd_tilde-Rd*log(pd/p_tilde)) - qt*(sv_tilde-Rv*log(pv/p_tilde)) ) / cp )
    elif pv == 0:
        temp = T_tilde * exp( (s-(1.0-qt)*(sd_tilde-Rd*log(pd/p_tilde)))/cp ) * exp( -qt*sv_tilde/cp ) * (pv/cp)**(pv/p_tilde)
    else:
        print('pv < 0: not physical!!')
    return temp

# Entropy Potential Temperature
def thetas_c(s, qt):
    return T_tilde*np.exp((s-(1.0-qt)*sd_tilde - qt*sv_tilde)/cpm_c(qt))

# ---------------------------------------------------------------------------
''' Latent Heat '''
# from Microphysics.pxd - L_fp(T,lambda)
def latent_heat(T):
    TC = T - 273.15
    return (2500.8 - 2.36 * TC + 0.0016 * TC *
            TC - 0.00006 * TC * TC * TC) * 1000.0

def latent_heat_constant(T, Lambda):
    return 2.501e6

def cpm_c(qt):
    return (1.0-qt) * cpd + qt * cpv

# ---------------------------------------------------------------------------
''' Condensation / Clausius Clapeyron '''
#def CC_const_Lv(T):
#    pv_star = 0.0
#    return pv_star

def CC_Magnus(T):
    T = T - 273.15
    pv_star = 6.1094 * np.exp((17.625*T)/(T+243.04))*100
    return pv_star

def qv_star_c(p0,qt,pv):
    pd = p0-pv
        #if pd < 0:
        #print('in qv_star_c pd<0 (pv='+str(pv)+')')
    return eps_v * (1.0 - qt) * pv / pd


# inversion of qv_star_c()
def pv_star_from_qv_star(p0,qt,qv):
    a = qv*eps_vi/(1.0-qt)
    pv = a*p0 / (1.0+a)
    return pv

