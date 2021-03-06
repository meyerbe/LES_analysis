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
    return p0 * eps_vi * qv / (1.0 - qt + eps_vi * qv)

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
    # print('s dry: ', type(sd_tilde), type(cpd), type(T), type(T_tilde), type(Rd), type(pd), type(p_tilde))
    temp = sd_tilde + cpd*log(T/T_tilde) - Rd*log(pd/p_tilde)
    # print(type(temp))
    return temp
    # return sd_tilde + cpd*log(T/T_tilde) - Rd*log(pd/p_tilde)

def s_vap(pv, T):
    temp = sv_tilde + cpv*log(T/T_tilde) - Rv * log(pv/p_tilde)
    # print('s vap: ', type(sv_tilde), type(cpv), type(T), type(T_tilde), type(Rv), type(pv), type(p_tilde))
    # print(type(temp))
    return temp
    # return sv_tilde + cpv*log(T/T_tilde) - Rv * log(pv/p_tilde)

def s_cond(L, T):
    # print('s cond: ', type(L), type(T), type(L/T))
    return -L/T
# ---------------------------------------------------------------------------
''' Temperature '''
# T(pd,pv,s,qt)
def temperature_no_ql_from_entropy(pd,pv,s,qt):
    cp = ((1.0 - qt) * np.double(cpd) + qt * cpv)
    # temp = T_tilde * np.exp(
    #     (s - (1.0 - qt) * (sd_tilde - Rd * log(pd / p_tilde)) - qt * (sv_tilde - Rv * log(pv / p_tilde))) / cp)
    if pv > 0:
        temp = T_tilde * np.exp( ( s - (1.0-qt)*(sd_tilde-Rd*log(pd/p_tilde)) - qt*(sv_tilde-Rv*log(pv/p_tilde)) ) / cp )
    elif pv == 0:
        temp = T_tilde * np.exp( (s-(1.0-qt)*(sd_tilde-Rd*log(pd/p_tilde)))/cp ) * np.exp( -qt*sv_tilde/cp ) * (pv/cp)**(pv/p_tilde)
    else:
        print('pv < 0: not physical!!')
        temp = 1
    # print('types temp no ql from entropy:', type(cp), type(qt), type(cpd), type(cpv))
    # print(type(T_tilde), type(s), type(pd), type(pv), type(p_tilde))
    # print(type(sd_tilde), type(sv_tilde), type(Rd), type(Rv))
    # print(type(temp))
    return temp

def temperature_no_ql_from_thetali(pd, pv, thl, qt):
    p0 = pd + pv
    # cp = ((1.0 - qt) * cpd + qt * cpv)
    # return thl * exner_c_moist(p0, Rd/cp)
    return thl * exner_c(p0)


# Exner Function
def exner_c(p0):
    return np.power((p0/p_tilde),kappa)

def exner_c_moist(p0, kap):
    return np.power((p0/p_tilde),kap)

# Dry potential temperature
def theta_c(p, T):
    return T / exner_c(p)

# Entropy Potential Temperature
def thetas_c(s, qt):
    return T_tilde*np.exp((s-(1.0-qt)*sd_tilde - qt*sv_tilde)/cpm_c(qt))

# Liquid ice potential temperature consistent with Triopoli and Cotton (1981)
# def thetali(p, T, qt, ql, qi, L):
def theta_li(p, T, qt, ql, qi):
    L = latent_heat(T)
    return theta_c(p, T) * np.exp(-latent_heat(T)*(ql/(1.0 - qt) + qi/(1.0 -qt))/(T*cpd))

# ---------------------------------------------------------------------------
''' Latent Heat '''
# from Microphysics.pxd - L_fp(T,lambda)
def latent_heat(T):
    TC = T - 273.15
    temp = (2500.8 - 2.36 * TC + 0.0016 * TC *
            TC - 0.00006 * TC * TC * TC) * 1000.0
    # print('latent heat: ', type(T), type(TC), type(temp))       # ok
    return temp
    # return (2500.8 - 2.36 * TC + 0.0016 * TC *
    #         TC - 0.00006 * TC * TC * TC) * 1000.0
    # return latent_heat_constant(T)

def latent_heat_constant(T):
    return 2.501e6

def cpm_c(qt):
    return (1.0-qt) * cpd + qt * cpv

# ---------------------------------------------------------------------------
''' Condensation / Clausius Clapeyron '''
#def CC_const_Lv(T):
#    pv_star = 0.0
#    return pv_star

def CC_Magnus(T):
    TC = T - 273.15
    pv_star = 6.1094 * np.exp((17.625*TC)/(TC+243.04))*100
    # print('CC Magnus: ', type(TC), type(pv_star))       # ok
    return pv_star

def qv_star_c(p0,qt,pv):
    pd = (p0-pv)
        #if pd < 0:
        #print('in qv_star_c pd<0 (pv='+str(pv)+')')
    temp = eps_v * (1.0 - qt) * pv / pd
    # print('qv star: ', type(eps_v), type(qt), type(pv), type(pd), type(temp))   # ok
    return temp
    # return eps_v * (1.0 - qt) * pv / pd


# inversion of qv_star_c()
def pv_star_from_qv_star(p0,qt,qv):
    a = qv*eps_vi/(1.0-qt)
    pv = a*p0 / (1.0+a)
    return pv

