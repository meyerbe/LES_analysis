''' Studying Clausiu Clapeyron Relationship: Magnus Form '''

import numpy as np
import pylab as plt
label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

from parameters import *
from entropies import sd_c, sv_c, sc_c
from thermo_aux import CC_Magnus
from thermo_aux import latent_heat
from CC_PyCLES import ClausiusClapeyron


def main():
    global p0
    global s, qt
    p0 = 1e5

#    ''' only CC Magnus: T --> pv_star '''
#    Tmin = 250
#    Tmax = 400
#    T = np.linspace(Tmin,Tmax,101)
#    print(T)
#    pv_star = CC_Magnus(T)
#    global pv_star_ref, T_ref
#    pv_star_ref, T_ref = ClausiusClapeyron(Tmin,Tmax)
#    plot_CC(T,pv_star)


    ''' (S,qt)->T & CC Magnus: (S,qt) --> T --> pv_star '''
    from thermo_aux import temperature_no_ql
    from thermo_aux import pv_c, qv_star_c
    ns = 21
    nqt = 11
    s = np.linspace(6500,7500,ns)
    qt = np.linspace(0.0,0.05,nqt)
    qv = np.zeros((nqt))
    pd = np.zeros((nqt))
    pv = np.zeros((nqt))
    T = np.zeros((ns,nqt))
    pv_star_1 = np.zeros((ns,nqt))
    qv_star_1 = np.zeros((ns,nqt))
    
    '''
        (1) starting point (first guess): pv, pd, T
            - assumption: no liquid, i.e. qv=qt
    '''
    qv = qt
    pv = pv_c(p0,qt,qv)
    pd = p0*np.ones((nqt))-pv
    plot_pd(T,pd)
    # --> pv, pd well behaved
    for i in xrange(ns):
        for j in xrange(nqt):
            #        T = temperature_no_ql(pd[i,j],pv[i,j],s[i],qt[j])
            T[i,j] = temperature_no_ql(pd[j],pv[j],s[i],qt[j])
            pv_star_1[i,j] = CC_Magnus(T[i,j])
            qv_star_1[i,j] = qv_star_c(p0,qt[j],pv_star_1[i,j])

    plot_temp(T)
    plot_CC_all(T,pv_star_1,qv_star_1)
    # T well defined for all s, qt
    # pv_star<p0 & negative qv_star_1 occur for high s (and low qt)
    
    

    '''
        (2) Check if saturated or not
    '''
    sigma_1 = np.ones((ns,nqt))
    sigma_1 = qt - qv_star_1
    s_1 = np.zeros((ns,nqt))
    pv_star_2 = np.zeros((ns,nqt))
    qv_star_2 = np.zeros((ns,nqt))
    T_2 = np.zeros((ns,nqt))
    delta_T12 = np.zeros((ns,nqt))
    # test variables
    sd = np.zeros((ns,nqt))
    sc = np.zeros((ns,nqt))
    for i in xrange(ns):
        for j in xrange(nqt):
            if (qt[j]<=qv_star_1[i,j]):
                print("not saturated: (s,qt)=", round(s[i],2), round(qt[j],4))
                ql = 0.0
                qi = 0.0
            else:
                print("saturated: (s,qt)=", round(s[i],2), round(qt[j],4))
#                print("saturated: (s,qt)=", round(s[i],2), round(qt[i,j],4), "pv=", pv[i,j])
#                T_1 = T[i,j]
#                pd_1 = pd[i,j]
#                pv_1 = pv[i,j]
#                lam_1 = 1.0
#                L_1 = latent_heat(T_1)
#
#                s_1[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j]) + sv_c(pv_1,T_1)*qt[i,j] + sc_c(L_1,T_1)*sigma_1[i,j]
#                # test
#                sd[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j])
#                sc[i,j] = sc_c(L_1,T_1)*sigma_1[i,j]
#                # --
    print(np.isnan(sigma_1).any())
#    print(sigma_1)
#    plot_sigma(sigma_1)





# ------------------------------------------------


def plot_1D(data):
    return
def plot_2D(field):
    return
def plot_sigma(sigma):
    plt.figure()
#    plt.contour(sigma,levels=np.linspace(-2e-1,1,10))
    plt.contourf(sigma)
    plt.colorbar()
#    plt.xlabel(r'$q_t$'+' moisture',fontsize=12)
#    plt.ylabel(r'$s$'+' entropy',fontsize=12)
#    plt.title('sigma 1',fontsize=15)
#    ax = plt.gca()
#    labels = ax.get_xticks()
#    for i in range(labels.shape[0]):
#        labels[i] = qt[int(labels[i])]
#    ax.set_xticklabels(labels)
#    labels = ax.get_yticks().astype(int)
#    for i in range(labels.shape[0]):
#        labels[i] = np.round(s[labels[i]], 2)
#    ax.set_yticklabels(labels)
    plt.savefig('figures/4_sigma1.png')
#    plt.close()
    return

def plot_CC_all(T,pv_star_1,qv_star_1):
    a = pv_star_1-p0
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    ax1 = plt.contourf(T)
    ax2 = plt.contour(pv_star_1-p0,levels=[0.0])
#    plt.text(qt[3],s[10],'jo',fontsize=20)
    plt.xlabel(r'$q_t$'+' moisture',fontsize=12)
    plt.ylabel(r'$s$'+' entropy',fontsize=12)
    plt.title('temperature (no ql)',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels = ax.get_xticks()
    for i in range(labels.shape[0]):
        labels[i] = qt[int(labels[i])]
    ax.set_xticklabels(labels)
    labels = ax.get_yticks().astype(int)
    for i in range(labels.shape[0]):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)
    plt.subplot(1,4,2)
    plt.contourf(pv_star_1)
    #plt.contour(1e5)
    plt.xlabel(r'$q_t$'+' moisture',fontsize=12)
    #plt.ylabel(r'$s$'+' entropy',fontsize=12)
    plt.title('Magnus Formula: pv_star',fontsize=15)
    plt.colorbar()
    ax = plt.gca()
    labels = ax.get_xticks()
    for i in range(labels.shape[0]):
        labels[i] = qt[int(labels[i])]
    ax.set_xticklabels(labels)
    labels = ax.get_yticks().astype(int)
    for i in range(labels.shape[0]):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)
    plt.subplot(1,4,3)
    ax1 = plt.contourf(pv_star_1-p0)
    cont = np.linspace(0.0,1)
    ax2 = plt.contour(pv_star_1-p0, levels = cont)
    #plt.contour(1e5)
    plt.xlabel(r'$q_t$'+' moisture',fontsize=12)
    #plt.ylabel(r'$s$'+' entropy',fontsize=12)
    plt.title('pv_star-p0',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels = ax.get_xticks()
    for i in range(labels.shape[0]):
        labels[i] = qt[int(labels[i])]
    ax.set_xticklabels(labels)
    labels = ax.get_yticks().astype(int)
    for i in range(labels.shape[0]):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)
    plt.subplot(1,4,4)
    ax1 = plt.contourf(qv_star_1)
    ax2 = plt.contour(qv_star_1,levels=[0.0],linewidth=2)
    plt.xlabel(r'$q_t$'+' moisture',fontsize=12)
    #plt.ylabel(r'$s$'+' entropy',fontsize=12)
    plt.title('qv_star_1',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels = ax.get_xticks()
    for i in range(labels.shape[0]):
        labels[i] = qt[int(labels[i])]
    ax.set_xticklabels(labels)
    labels = ax.get_yticks().astype(int)
    for i in range(labels.shape[0]):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)
    #plt.show()
    plt.savefig('figures/3_CC_all_Magnus.png')
#    plt.close()
    return
def plot_pd(T,pd):
    plt.figure()
    plt.plot(pd,linewidth=2)
    plt.ylabel(r'$p_d$',fontsize=21)
    plt.savefig('figures/1_pd.png')
#    plt.close()
    return
def plot_temp(T):
    plt.figure()
    plt.contourf(T.T)
    plt.xlabel(r'$s$'+' entropy',fontsize=18)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=18)
    plt.title('temperature (no ql)',fontsize=21)
    plt.colorbar()
    plt.savefig('figures/2_temp.png')
#    plt.close()
    return
def plot_CC(T,pv_star):
    plt.figure()
    plt.plot(T,pv_star,linewidth=2)
    plt.plot([373,373],[0.1,max(pv_star)])
    plt.plot([min(T),max(T)],[1e5,1e5],linewidth=2)
    plt.grid(b=True, which='major', color='grey', linestyle='-')
    plt.grid(b=True, which='minor', color='grey', linestyle='--')
    plt.text(270, 1e2,'Magnus Form',color='blue')
    plt.text(np.min(T)+30, 1.2e5,'p0=10^5',color='red')
    plt.text(373, 1e2,'T=373 K',color='green')
    plt.yscale('log')
    plt.ylabel('pv_star',fontsize=18)
    plt.xlabel('temperature [K]',fontsize=18)
    plt.title('Magnus Formula',fontsize=24)
    #plt.show()
    plt.savefig('figures/CC_Magnus.png')
    plt.close()
    return

# ------------------------------------------------
if __name__ == '__main__':
    main()

