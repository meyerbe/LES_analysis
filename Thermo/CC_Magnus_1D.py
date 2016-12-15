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
    nqt = 1
    s = np.linspace(6800,7200,ns)
    qt = np.zeros((ns,nqt))
    qv = np.zeros((ns,nqt))
    pd = np.zeros((ns,nqt))
    pv = np.zeros((ns,nqt))
    T = np.zeros((ns,nqt))
    pv_star_1 = np.zeros((ns,nqt))
    qv_star_1 = np.zeros((ns,nqt))
    
    '''
    (1) starting point (first guess): pv, pd, T
        - assumption: no liquid, i.e. qv=qt
    '''
    global nan_index, nan_min_s
    nan_index = np.zeros((ns,nqt))
    nan_min_s = 9999
    
    qv = qt
    pv = pv_c(p0,qt,qv)
    pd = p0*np.ones((nqt))-pv
    plot_pd(T,pd)
    # --> pv, pd well behaved
    for i in xrange(ns):
        for j in xrange(nqt):
            T[i,j] = temperature_no_ql(pd[i,j],pv[i,j],s[i],qt[i,j])
            pv_star_1[i,j] = CC_Magnus(T[i,j])
            qv_star_1[i,j] = qv_star_c(p0,qt[i,j],pv_star_1[i,j])
            if qv_star_1[i,j]<0:        # meaning problem
                nan_index[i,j] = -1
                if i<nan_min_s:
                    nan_min_s = i
            else:                       # meaning all ok
                nan_index[i,j] = 1

    plot_temp(T)
    plot_CC_all_1(T,pv_star_1,qv_star_1,'3_CC_all_Magnus_firstguess')
    # T well defined for all s, qt
    # pv_star<p0 & negative qv_star_1 occur for high s (and low qt)
    
    
    '''
        (2) Check if saturated or not
    '''
    sigma_1 = qt - qv_star_1        # = - qv_star_1 > 0
    s_1 = np.zeros((ns,nqt))
    pv_star_2 = np.zeros((ns,nqt))
    qv_star_2 = np.zeros((ns,nqt))
    pv_star_3 = np.zeros((ns,nqt))
    qv_star_3 = np.zeros((ns,nqt))
    T_2 = np.zeros((ns,nqt))
    delta_T12 = np.zeros((ns,nqt))
    delta_T23 = np.zeros((ns,nqt))
    # test variables
    sd = np.zeros((ns,nqt))
    sc = np.zeros((ns,nqt))
    global sat_index
    sat_index = np.zeros((ns,nqt))
    for i in xrange(ns):
        for j in xrange(nqt):
            if (qt[i,j]<=qv_star_1[i,j] or qv_star_1[i,j]<0.0):
                print("not saturated: (s,qt)=", round(s[i],2), round(qt[i,j],4))
                sat_index[i,j] = -1
                ql = 0.0
                qi = 0.0
            else:
                print("saturated: (s,qt)=", round(s[i],2), round(qt[i,j],4), "pv=", pv[i,j])
                ''' ---- T_2 ---- '''
                sat_index[i,j] = 1
                T_1 = T[i,j]
                pd_1 = pd[i,j]
                pv_1 = pv[i,j]
                lam_1 = 1.0
                L_1 = latent_heat(T_1)
                if pv[i,j] == 0:
                    s_1[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j])
#                    T_2[i,j] = T_1 + sigma_1[i,j]*L_1 / ((1.0-qt[i,j])*cpd + qv_star_1[i,j]*cpv)
                    # continue
                else:
                    s_1[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j]) + sv_c(pv_1,T_1)*qt[i,j] + sc_c(L_1,T_1)*sigma_1[i,j]
#                    T_2[i,j] = T_1 + sigma_1[i,j]*L_1 / ((1.0-qt[i,j])*cpd + qv_star_1[i,j]*cpv)
                # s_1[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j]) + sv_c(pv_1,T_1)*qt[i,j] + sc_c(L_1,T_1)*sigma_1[i,j]
                print('qt',qt[i,j])
                # test
                sd[i,j] = sd_c(pd_1,T_1)*(1.0-qt[i,j])
                sc[i,j] = sc_c(L_1,T_1)*sigma_1[i,j]
                # --
                f_1 = s[i] - s_1[i,j]
                T_2[i,j] = T_1 + sigma_1[i,j]*L_1 / ((1.0-qt[i,j])*cpd + qv_star_1[i,j]*cpv)
                delta_T12[i,j] = np.fabs(T_2[i,j]-T_1)
                
                pv_star_2[i,j] = CC_Magnus(T_2[i,j])
                qv_star_2[i,j] = qv_star_c(p0, qt[i,j], pv_star_2[i,j])
                ql_2 = qt[i,j] - qv_star_2[i,j]
                lam_2 = 1.0

                count = 0
                if (delta_T12[i,j]>= 1.0e-3 or ql_2<0.0):
                    print("start while loop")
                    ''' ---- T_3 ---- '''
                    pv_star_3[i,j] = CC_Magnus(T_2[i,j])
                    qv_star_3[i,j] = qv_star_c(p0,qt[i,j],pv_star_3[i,j])
                    


    plot_sigma(sigma_1)
    plot_s(s_1,sd,sc)
    plot_CC_all_2(T,T_2,delta_T12,pv_star_1,qv_star_1,pv_star_2,qv_star_2,'6_CC_all_Magnus_secondguess')







# ------------------------------------------------


def plot_1D(data):
    return
def plot_2D(field):
    return
#def plot_temp(T):
#    plt.figure()
#    plt.plot(s,T[:,0])
#    plt.xlabel(r'$s$'+' entropy',fontsize=18)
#    plt.ylabel(r'$T$'+' temperature',fontsize=18)
#    plt.title('temperature (no ql)',fontsize=21)
#    plt.savefig('figures/1D_2_temp.png')
#    #    plt.close()
#    return
def plot_s(s_,sd_,sc_):
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(s,s,label='s init')
    plt.plot(s,s_,label='s new')
    plt.plot(s,sd_,'--',label='s dry')
    plt.plot(s,sc_,'--',label='s cond')
    plt.legend(loc=3,fontsize=9)
    plt.title('initial s, new s')
    plt.subplot(1,2,2)
    plt.plot(s,s_[:,0],'x-',label='s_{new}')
    plt.plot(s,s,'-',label='s_{init}')
    plt.plot(s,1e3*sat_index[:,0],'g',label='sat=1,unsat=-1')
    plt.plot(s,1e3*nan_index[:,0],'-o',label='qv* pos=1,qv* neg=-1')
    plt.plot([s[nan_min_s-1],s[nan_min_s-1]],[min(s_),max(s_)],'k--')
    plt.legend(loc=3,fontsize=8)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$s_{new}$',fontsize=12)
    plt.title(r'$s_{new}$',fontsize=15)
    plt.savefig('figures/1D_5_s1.png')
    plt.close()
def plot_sigma(sigma):
    global sat_index, nan_index
    plt.figure()
    plt.plot(s,sigma[:,0],label='sigma_1')
    plt.plot([s[0],s[-1]],[0,0],'k--')
    plt.plot([s[14],s[14]],[min(sigma),max(sigma)],'k--')
    plt.plot(s,sat_index[:,0],'g',label='sat=1,unsat=-1')
    plt.plot(s,nan_index[:,0],'-o',label='qv* pos=1,qv* neg=-1')
    plt.legend(loc=3)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$\sigma_1$',fontsize=12)
    plt.title(r'$\sigma_1 = q_t - q_{v,1}^*$',fontsize=15)
    plt.savefig('figures/1D_4_sigma1.png')
    plt.close()
    return
def plot_CC_all_1(T,pv_star_1,qv_star_1,message):
    a = pv_star_1-p0
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    ax1 = plt.plot(s,T[:])
    plt.plot([s[14],s[14]],[T[0],T[-1]],'k--')
    plt.plot([s[0],s[-1]],[T[nan_min_s-1],T[nan_min_s-1]],'k--')
    #    plt.text(qt[3],s[10],'jo',fontsize=20)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel('temperature',fontsize=12)
    plt.title('temperature (no ql, qt='+str(qt[0,0])+')',fontsize=15)
    plt.subplot(1,4,2)
    plt.plot(s,pv_star_1[:,0])
    plt.plot([s[0],s[-1]],[1e5,1e5],'k--')
    plt.plot([s[14],s[14]],[pv_star_1[0],pv_star_1[-1]],'k-')
    plt.plot([s[nan_min_s-1],s[nan_min_s-1]],[pv_star_1[0],pv_star_1[-1]],'g-')
    plt.plot([s[nan_min_s],s[nan_min_s]],[pv_star_1[0],pv_star_1[-1]],'g-')
    plt.text(s[15],pv_star_1[0],'s[14]='+str(s[14]),fontsize=8)
    plt.yscale('log')
    plt.ylabel(r'$p_v^*$')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('Magnus Formula: '+r'$p_v^*$',fontsize=15)
    plt.subplot(1,4,3)
    ax1 = plt.plot(s,pv_star_1[:,0]-p0)
    plt.plot([s[0],s[-1]],[0,0],'k--')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$p_{v,1}^*-p_0$')
    plt.title(r'$p_{v,1}^*-p_0$',fontsize=15)
    plt.subplot(1,4,4)
    ax1 = plt.plot(s,qv_star_1[:,0],'-x')
    plt.plot([s[14],s[14]],[min(qv_star_1),max(qv_star_1)],'k--')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_{v,1}^*$'+' saturation moisture',fontsize=12)
    plt.title(r'$q_{v,1}^*$',fontsize=15)
    plt.savefig('figures/1D_'+message+'.png')
    plt.close()
    return
def plot_CC_all_2(T,T_new,deltaT,pv_star_1,qv_star_1,pv_star_2,qv_star_2,message):
    a = pv_star_1-p0
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    plt.plot(s,T[:],label='T init')
    plt.plot(s,T_new[:],label='T_2',linewidth=3)
    plt.plot(s,deltaT[:],'-x',label=r'$\Delta T$')
    plt.plot([s[nan_min_s-1],s[nan_min_s-1]],[min(T_new),max(deltaT)],'k--')
    plt.plot([s[0],s[-1]],[T[nan_min_s-1],T[nan_min_s-1]],'k--')
    plt.legend(loc=3,fontsize=9)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel('temperature',fontsize=12)
    plt.title('temperature (no ql, qt='+str(qt[0,0])+')',fontsize=15)
    plt.subplot(1,4,2)
    plt.plot(s,pv_star_1[:,0],label=r'$p_{v,1}^*$')
    plt.plot(s,pv_star_2[:,0],label=r'$p_{v,2}^*$')
    plt.plot([s[0],s[-1]],[1e5,1e5],'k--')
    plt.plot([s[14],s[14]],[pv_star_2[0],max(pv_star_2)],'k--')
    plt.text(s[15],pv_star_1[0],'s[14]='+str(s[14]),fontsize=8)
    plt.legend(loc=2,fontsize=8)
    plt.yscale('log')
    plt.ylabel(r'$p_v^*$')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('Magnus Formula: '+r'$p_v^*$',fontsize=15)
    plt.subplot(1,4,3)
    plt.plot(s,pv_star_2[:,0]-p0,label=r'$p_{v,2}^*-p0')
    plt.plot(s,pv_star_1[:,0]-p0,label=r'$p_{v,1}^*-p0')
    plt.plot([s[0],s[-1]],[0,0],'k--')
    plt.plot([s[14],s[14]],[pv_star_2[0],max(pv_star_2)],'k--')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$p_{v,1}^*-p_0$')
    plt.title(r'$p_{v,1}^*-p_0$',fontsize=15)
    plt.subplot(1,4,4)
    ax1 = plt.plot(s,qv_star_1[:,0],'-x')
    plt.plot([s[14],s[14]],[min(qv_star_1),max(qv_star_1)],'k--')
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_{v,1}^*$'+' saturation moisture',fontsize=12)
    plt.title(r'$q_{v,1}^*$',fontsize=15)
    plt.savefig('figures/1D_'+message+'.png')
    plt.close()
    return
def plot_pd(T,pd):
    plt.figure()
    plt.plot(s,pd[:,0],linewidth=2)
    plt.xlabel('entropy',fontsize=21)
    plt.ylabel(r'$p_d$',fontsize=21)
    plt.savefig('figures/1D_1_pd.png')
    #    plt.close()
    return
def plot_temp(T):
    plt.figure()
    plt.plot(s,T[:,0])
    plt.xlabel(r'$s$'+' entropy',fontsize=18)
    plt.ylabel(r'$T$'+' temperature',fontsize=18)
    plt.title('temperature (no ql)',fontsize=21)
    plt.savefig('figures/1D_2_temp.png')
    #    plt.close()
    return
def plot_CC(T,pv_star):
    plt.figure()
    plt.plot(T,pv_star,linewidth=2,label='Magnus')
    plt.plot(T_ref,pv_star_ref,'--',linewidth=2,label='Ref (PyCLES)')
    plt.plot([373,373],[0.1,max(pv_star)],'k')
    plt.plot([min(T),max(T)],[1e5,1e5],'r',linewidth=2)
    plt.legend(loc=3)
    plt.xlim([max(T[0],T_ref[0]),max(T[-1],T_ref[-1])])
    plt.ylim([5e1,1e6])
    plt.grid(b=True, which='major', color='grey', linestyle='-')
    plt.grid(b=True, which='minor', color='grey', linestyle='--')
    #plt.text(270, 1e2,'Magnus Form',color='blue')
    plt.text(np.min(T)+30, 1.2e5,'p0=10^5',color='red')
    plt.text(373, 1e2,'T=373 K',color='green')
    plt.yscale('log')
    plt.ylabel('pv_star',fontsize=18)
    plt.xlabel('temperature [K]',fontsize=18)
    plt.title('Magnus Formula',fontsize=24)
    plt.savefig('figures/CC_Magnus.png')
    plt.close()
    return

# ------------------------------------------------
if __name__ == '__main__':
    main()

