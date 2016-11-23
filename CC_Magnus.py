''' Studying Clausiu Clapeyron Relationship: Magnus Form '''

import numpy as np
import pylab as plt

from thermo_aux import CC_Magnus


def main():
    global p0
    global s, qt
    p0 = 1e5

    ''' only CC Magnus: T --> pv_star '''
    T = np.linspace(0,500,501)
    T = np.linspace(200,400,101)
    print(T)
    pv_star = CC_Magnus(T)
    plot_CC(T,pv_star)


    ''' (S,qt)->T & CC Magnus: (S,qt) --> T --> pv_star '''
    from thermo_aux import temperature_no_ql
    from thermo_aux import pv_c, qv_star_c
    ns = 101
    nqt = 11
    s = np.linspace(6000,7500,ns)
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
    print(pd.shape)
    print(pv.shape)
    for i in xrange(ns):
        for j in xrange(nqt):
            #        T = temperature_no_ql(pd[i,j],pv[i,j],s[i],qt[j])
            T[i,j] = temperature_no_ql(pd[j],pv[j],s[i],qt[j])
            pv_star_1[i,j] = CC_Magnus(T[i,j])
            qv_star_1[i,j] = qv_star_c(p0,qt[j],pv_star_1[i,j])

    plot_temp(T)
    plot_CC_all(T,pv_star_1,qv_star_1)
    # T well defined for all s, qt
    
    
    
    

    '''
        (2) 
    '''




# ------------------------------------------------


def plot_1D(data):
    return
def plot_2D(field):
    return
def plot_CC_all(T,pv_star_1,qv_star_1):
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.contourf(T.T)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.title('temperature (no ql)',fontsize=18)
    plt.colorbar()
    plt.subplot(1,3,2)
    plt.contourf(pv_star_1.T)
    plt.contour(1e5)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.title('Magnus Formula: pv_star',fontsize=18)
    plt.colorbar()
    plt.subplot(1,3,3)
    plt.contourf(qv_star_1.T)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.title('qv_star_1',fontsize=18)
    plt.colorbar()
    #plt.show()
    plt.savefig('figures/3_CC_all_Magnus.png')
    return
def plot_pd(T,pd):
    plt.figure()
    plt.plot(pd,linewidth=2)
    plt.ylabel(r'$p_d$',fontsize=21)
    plt.savefig('figures/1_pd.png')
    plt.close()
    return
def plot_temp(T):
    plt.figure()
    plt.contourf(T.T)
    plt.xlabel(r'$s$'+' entropy',fontsize=18)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=18)
    plt.title('temperature (no ql)',fontsize=21)
    plt.colorbar()
    plt.savefig('figures/2_temp.png')
    plt.close()
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
    return

# ------------------------------------------------
if __name__ == '__main__':
    main()

