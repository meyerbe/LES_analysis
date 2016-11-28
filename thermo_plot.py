''' Auxiliary thermodynamic functions from LES '''
import pylab as plt
import numpy as np
from parameters import *
from math import *

from thermo_aux import pv_c
from thermo_aux import CC_Magnus
from CC_PyCLES import ClausiusClapeyron

def main():
    plot_CC()
    return

''' Clausius Clapeyron: Magnus vs. PyCLES: T --> pv_star '''
def plot_CC():
    Tmin = 250
    Tmax = 400
    T = np.linspace(Tmin,Tmax,101)
    print(T)
    pv_star = CC_Magnus(T)
    global pv_star_ref, T_ref
    pv_star_ref, T_ref = ClausiusClapeyron(Tmin,Tmax)

    plt.figure()
    plt.plot(T,pv_star,linewidth=2,label='Magnus')
    plt.plot(T_ref,pv_star_ref,'--',linewidth=2,label='Ref (PyCLES)')
    plt.plot([373,373],[min(pv_star),max(pv_star)],'k')
    plt.plot([min(T),max(T)],[1e5,1e5],'r',linewidth=2)
    plt.legend(loc=3)
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
    plt.savefig('figures/CC_Magnus_log.png')
    plt.close()
    
    plt.figure()
    plt.plot(T,pv_star,linewidth=2,label='Magnus')
    plt.plot(T_ref,pv_star_ref,'--',linewidth=2,label='Ref (PyCLES)')
    plt.plot([373,373],[min(pv_star),max(pv_star)],'k')
    plt.plot([min(T),max(T)],[1e5,1e5],'r',linewidth=2)
    plt.legend(loc=3)
    plt.grid(b=True, which='major', color='grey', linestyle='-')
    plt.grid(b=True, which='minor', color='grey', linestyle='--')
    plt.text(345, 1.2e5,'Magnus Form',color='blue')
    plt.text(np.min(T)+30, 1.1e5,'p0=10^5',color='red')
    plt.text(373, 1e4,'T=373 K',color='green')
    plt.ylabel('pv_star',fontsize=18)
    plt.xlabel('temperature [K]',fontsize=18)
    plt.title('Magnus Formula',fontsize=24)
    #plt.show()
    plt.savefig('figures/CC_Magnus.png')
    plt.close()
    return


''' Vapor pressure '''
def plot_pvsat(T,pv_star):
    qv = np.linspace(0,0.5,11)
    qt_ = np.linspace(0.0,0.8,5)
    #qt = 0.3
    p0 = 1e5
    plt.figure()
    for qt in qt_:
        plt.plot(qv,pv_c(p0,qt,qv),label='qt='+str(qt))
    plt.legend()
    plt.xlabel(r'$q_v$')
    plt.ylabel('vapor pressure '+r'$p_v$')
    plt.title('vapor pressure')
    plt.savefig('figures/vapor_pressure.png')

''' Liquid potential temperature '''


# ------------------------------------------------
if __name__ == '__main__':
    main()