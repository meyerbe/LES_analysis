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
    global s, qt, ns, nqt
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
    ns = 150
    nqt = 100
    s = np.linspace(6800,7300,ns)
    qt = np.linspace(0.0,0.05,nqt)
    qv = np.zeros((ns,nqt))
    pd_1 = np.zeros((ns,nqt))
    pv_1 = np.zeros((nqt))
    T = np.zeros((ns,nqt))

    
    s_1 = np.zeros((ns,nqt))
    s_2 = np.zeros((ns,nqt))
    s_3 = np.zeros((ns,nqt))
    pv_star_1 = np.zeros((ns,nqt))
    qv_star_1 = np.zeros((ns,nqt))
    pv_star_2 = np.zeros((ns,nqt))
    qv_star_2 = np.zeros((ns,nqt))
    pv_star_3 = np.zeros((ns,nqt))
    qv_star_3 = np.zeros((ns,nqt))
    pv_star_4 = np.zeros((ns,nqt))
    qv_star_4 = np.zeros((ns,nqt))
    sigma_1 = np.zeros((ns,nqt))
    T_1 = np.zeros((ns,nqt))
    T_2 = np.zeros((ns,nqt))
    T_3 = np.zeros((ns,nqt))
    T_4 = np.zeros((ns,nqt))
    delta_T21 = np.zeros((ns,nqt))
    delta_T32 = np.zeros((ns,nqt))
    delta_T43 = np.zeros((ns,nqt))
    pv_2 = np.zeros((ns,nqt))
    pd_2 = np.zeros((ns,nqt))
    ql_2 = np.zeros((ns,nqt))
    pv_3 = np.zeros((ns,nqt))
    pd_3 = np.zeros((ns,nqt))
    ql_3 = np.zeros((ns,nqt))
    
    # -- test variables
    sd = np.zeros((ns,nqt))
    sv = np.zeros((ns,nqt))
    sc = np.zeros((ns,nqt))
    cpd_ = np.zeros((ns,nqt))
    # --
    global nan_index, sat_index, nan_index2
    global plt_count
    plt_count = 1
    nan_index = np.zeros((ns,nqt))
    sat_index = np.zeros((ns,nqt))
    nan_index2 = np.zeros((ns,nqt))
    nan_index3 = np.zeros((ns,nqt))
    count = np.zeros((ns,nqt))
    #
            
    
    
    '''
        (1) starting point (first guess): pv, pd, T
            - assumption: no liquid, i.e. qv=qt
    '''
    for i in xrange(ns):
        qv[i,:] = qt[:]
    pv_1 = pv_c(p0,qt,qv)
    pd_1 = p0*np.ones((nqt))-pv_1
    # --> pv, pd well behaved
    for i in xrange(ns):
        for j in xrange(nqt):
            T_1[i,j] = temperature_no_ql(pd_1[i,j],pv_1[i,j],s[i],qt[j])
            pv_star_1[i,j] = CC_Magnus(T_1[i,j])
            qv_star_1[i,j] = qv_star_c(p0,qt[j],pv_star_1[i,j])
            if qv_star_1[i,j]<0:        # meaning problem
                nan_index[i,j] = 1
            else:                       # meaning all ok
                nan_index[i,j] = 0

    # RESULTS:
    # T well defined for all s, qt
    # pv_star<p0 & negative qv_star_1 occur for high s (and low qt)

    

    '''
        (2) Check if saturated or not
    '''
    for i in xrange(ns):
        for j in xrange(nqt):
            if (qt[j]<=qv_star_1[i,j]):
            # if (qt[j]<=qv_star_1[i,j] or qv_star_1[i,j]<0.0):
                # print("not saturated: (s,qt)=", round(s[i],2), round(qt[j],4))
                sat_index[i,j] = 0
                ql = 0.0
                qi = 0.0
            else:
                '''
                    (3) if saturated: calculate second starting point T_2
                '''
                print("saturated: (s,qt)=", round(s[i],2), round(qt[j],4), "pv=", pv_1[i,j])
                sat_index[i,j] = 1.0
                ''' ---- s_1 ---- '''
                lam_1 = 1.0
                L_1 = latent_heat(T_1[i,j])
                sigma_1[i,j] = qt[j] - qv_star_1[i,j]
                s_1[i,j] = sd_c(pd_1[i,j],T_1[i,j])*(1.0-qt[j]) + sv_c(pv_1[i,j],T_1[i,j])*qt[j] + sc_c(L_1,T_1[i,j])*sigma_1[i,j]
                # test
                sd[i,j] = sd_c(pd_1[i,j],T_1[i,j])*(1.0-qt[j])
                sv[i,j] = sv_c(pv_1[i,j],T_1[i,j])*qt[j]
                sc[i,j] = sc_c(L_1,T_1[i,j])*sigma_1[i,j]
                cpd_[i,j] = (1.0-qt[j])*cpd + qv_star_1[i,j]*cpv
                # --
                f_1 = s[i] - s_1[i,j]
                ''' ---- T_2 ---- '''
                T_2[i,j] = T_1[i,j] + sigma_1[i,j]*L_1 / ((1.0-qt[j])*cpd + qv_star_1[i,j]*cpv)
                delta_T21[i,j] = np.fabs(T_2[i,j]-T_1[i,j])

                pv_star_2[i,j] = CC_Magnus(T_2[i,j])
                qv_star_2[i,j] = qv_star_c(p0, qt[j], pv_star_2[i,j])
                ql_2[i,j] = qt[j] - qv_star_2[i,j]
                lam_2 = 1.0

                count[i,j] = 0
                if (delta_T21[i,j]>= 1.0e-3 or ql_2[i,j]<0.0):
                    print("start while loop (ql_2="+ str(ql_2[i,j])+")" )
                    pv_star_2[i,j] = CC_Magnus(T_2[i,j])
                    qv_star_2[i,j] = qv_star_c(p0,qt[j],pv_star_2[i,j])
                    ''' ---- s_2 ---- '''
                    pv_2[i,j] = pv_c(p0,qt[j],qv_star_2[i,j])
                    pd_2[i,j] = p0 - pv_2[i,j]
                    # --> qv_star_2 and hence pd_2 can be negative !!!!
                    lam_2 = 1.0
                    L_2 = latent_heat(T_2[i,j])
                    if 1 == 1:
                    # if qv_star_2[i,j]>0:    # equivalent to "if pd_2<0"; problem: sd(pd) not well-defined for pd<0 since sd~log(pd)
                        s_2[i,j] = sd_c(pd_2[i,j],T_2[i,j])*(1.0 - qt[j]) + sv_c(pv_2[i,j],T_2[i,j])*qt[j] + sc_c(L_2,T_2[i,j])*ql_2[i,j]
                        f_2 = s[i] - s_2[i,j]
                        ''' ---- T_3 ---- '''
                        T_3[i,j] = T_2[i,j] - f_2*(T_2[i,j] - T_1[i,j])/(f_2 - f_1)
                        delta_T32[i, j] = np.abs(T_3[i, j] - T_2[i, j])
                        count[i, j] += 1
                    ''''''
                    if qv_star_2[i,j] < 0:
                        nan_index2[i,j] = 1
                    ''''''
                    # else:
                    #     nan_index2[i,j] = 1
                    #     ''' !!!!!!!?????????!!!!!!!!!!!!! '''
                    #     qv_star_2 = 0.0
                    #     a = qv_star_2*eps_vi/(1.0-qt)
                    #     pv_star_2 = a*p0 / (1.0+a)
                    #     T3 = T_2 + sigma_1 * L_1/((1.0 - qt)*cpd + qv_star_2 * cpv)
                    #     delta_T = np.fabs(T_2 - T_1)
                    #     ql_2 = qt-qv_star_2
                    #                    delta_T = 1e-5
                    #     ''' !!!!!!!?????????!!!!!!!!!!!!! '''
                    #
                    #     ''' ???? --> need to define T_3!! what to do?? --> new definition of qv_star_2??? '''
                    #     print('!! qv_star_2<0')
                    #     T_3[i,j] = T_2[i,j]
                    #     delta_T32[i,j]  = np.abs(T_3[i,j] - T_2[i,j])       # will be zero
            
                if (delta_T32[i,j] >= 1.0e-3 or ql_2[i,j]<0.0):
                    pv_star_3[i,j] = CC_Magnus(T_3[i,j])
                    qv_star_3[i,j] = qv_star_c(p0,qt[j],pv_star_3[i,j])
                    ''' ---- s_3 ---- '''
                    pv_3[i,j] = pv_c(p0,qt[j],qv_star_3[i,j])
                    pd_3[i,j] = p0 - pv_3[i,j]
                    
                    lam_2 = 1.0
                    L_3 = latent_heat(T_2[i,j])
                    if 1 == 1:
                    # if qv_star_3[i,j]>0:
                        s_3[i,j] = sd_c(pd_3[i,j],T_3[i,j])*(1.0 - qt[j]) + sv_c(pv_3[i,j],T_3[i,j])*qt[j] + sc_c(L_3,T_3[i,j])*ql_3[i,j]
                        f_3 = s[i] - s_3[i,j]
                        ''' ---- T_4 ---- '''
                        T_4[i,j] = T_3[i,j] - f_3*(T_3[i,j] - T_2[i,j])/(f_3 - f_2)
                        delta_T43[i,j]  = np.abs(T_4[i,j] - T_3[i,j])
                        count[i,j] += 1
                    # else:
                    #     nan_index3 [i,j] = 1
                    #     print('!! qv_star_3<0')
                    ''''''
                    if qv_star_3[i, j] < 0:
                        nan_index3[i, j] = 1
                    ''''''



    plot_CC_all(T_1,pv_star_1,qv_star_1)
    plot_s(s_1,sd,sv,sc,'s1')
    plot_sigma(sigma_1)
    plot_cp(cpd_,'cpd_tilde')
    plot_1D(L_1*sigma_1/cpd_,'sigma_1_cpd')
    plot_CC_all_2(T_2,pv_star_2,qv_star_2,pd_2)
    plot_temp(delta_T21,'deltaT_12')
    plot_sat(sat_index,nan_index,pd_1,'index1')
    plot_CC_all_3(T_3,pv_star_3,qv_star_3,'T3')
    plot_temp(delta_T32,'deltaT_23')
    plot_sat(sat_index,nan_index2,pd_2,'index2')
#    plot_temp(-f_2*(T_2-T_1)/(f_2-f_1))
    plot_s(s_3, sd, sv, sc,'s3')
    plot_CC_all_3(T_4,pv_star_4,qv_star_4,'T4')
    plot_sat(sat_index,nan_index3,pd_3,'index3')
    plot_1D(count,'iterations')




# ------------------------------------------------


def plot_1D(data,title):
    global plt_count
    plt.figure()
    plt.contourf(data.T)
    plt.colorbar()
    plt.title(title)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.savefig('figures/'+str(plt_count)+'_'+title+'.png')
    plt.close()
    plt_count += 1
    return
def plot_2D(field):
    return
def plot_cp(data,title):
    global plt_count
    plt.contourf(data.T,levels=np.append([0.0,1.0],np.linspace(900,1050,16)))
    plt.colorbar()
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.savefig('figures/'+str(plt_count)+'_'+title+'.png')
    plt.close()
    plt_count += 1
def plot_sat(sat_index,nan_index,pd_,name):
    global ns, plt_count
    plt.figure(figsize=(12,4))
    plt.subplot(1,3,1)
    ax1 = plt.contourf(pd_.T)
    ax2 = plt.contour(pd_.T,levels=[0.0],linewidths=2)
    plt.clabel(ax2,inline=1)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.title('pd_'+name[-1])
    plt.subplot(1,3,2)
    ax1 = plt.contourf(nan_index.T)
    ax2 = plt.contour(pd_.T,levels=[0.0])
    plt.clabel(ax2,inline=1)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.title('nan-'+name+' (1=qt*<0,0=qt*>0')
    plt.subplot(1,3,3)
    ax1 = plt.contourf(sat_index.T)
    ax2 = plt.contour(nan_index.T,'k',levels=[-1.0,1.0],linewidth=5)
#    ax3 = plt.contourf(nan_index, hatch="/")
#    ax3 = plt.fill(nan_index, fill=False, hatch='\\')
    plt.clabel(ax2,inline=1,fontsize=10)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.title('sat-index1 (0=unsat,1=sat)')
    plt.savefig('figures/'+str(plt_count)+'_sat_index.png')
    plt.close()
    plt_count += 1
    return

def plot_s(s_,sd_,sv_,sc_,name):
    global nan_index, ns, plt_count
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    ax1 = plt.contourf(s_.T)
    ax2 = plt.contour(sd_.T+sc_.T,levels=[0.0],linewidths=2)
    plt.clabel(ax2, inline=1, fontsize=10)#, text=r'$\alpha$')
    ax3 = plt.contour(nan_index.T,levels=[-1.0,1.0],linewidths=1)
    plt.clabel(ax3, inline=1, fontsize=10)#, text=r'$\alpha$')
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.title(name +' (cont=sd+sc, nan_index)')
    plt.subplot(1,4,2)
    ax1 = plt.contourf(sd_.T)
    ax2 = plt.contour(nan_index.T,levels=[-1.0,1.0],linewidth=2)
    plt.clabel(ax2, inline=1, fontsize=10)#, text=r'$\alpha$')
    plt.colorbar(ax1)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.title('s dry (cont=nan-ind)')
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    plt.subplot(1,4,3)
    ax1 = plt.contourf(sv_.T)
    ax2 = plt.contour(s_.T,levels=[-4e5,0.0],colors='k',linewidths=2)
    plt.clabel(ax2, inline=1, fontsize=10)#, text=r'$\alpha$')
    plt.colorbar(ax1)
    plt.title('s vap (cont=s1)',fontsize=15)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    plt.subplot(1,4,4)
    ax1 = plt.contourf(sc_.T)
    ax2 = plt.contour(s_.T,levels=[-4e5,0.0],colors='k',linewidths=2)
    plt.clabel(ax2, inline=1, fontsize=10)#, text=r'$\alpha$')
    plt.colorbar(ax1)
    plt.title('s liquid (cont=s1)',fontsize=15)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    plt.savefig('figures/'+str(plt_count)+'_s1.png')
    plt.close()
    plt_count += 1
    return

def plot_sigma(sigma):
    global ns, plt_count
    plt.figure()
#    plt.contour(sigma,levels=np.linspace(-2e-1,1,10))
    plt.contourf(sigma.T)
    plt.xlabel(r'$s_{init}$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$ moisture',fontsize=12)
    plt.colorbar()
    plt.savefig('figures/'+str(plt_count)+'_sigma1.png')
    plt.close()
    plt_count += 1
    return

def plot_CC_all_3(T,pv_star_,qv_star_,name):
    global ns, nan_index2, plt_count
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    # ax1 = plt.contourf(T.T,levels=np.linspace(np.amin(T),np.amax(T),250))
    ax1 = plt.contourf(T.T)
    ax2 = plt.contour(p0-pv_star_.T,levels=[0.0])
    # ax2 = plt.contour(nan_index2.T,levels=[0.0,1.0])
    ax3 = plt.contour(T.T,levels=[int(250),int(300)],colors='w')
    plt.colorbar(ax1)
    plt.clabel(ax2)
    plt.clabel(ax3,inline=1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.title('temperature '+name+' (no ql)',fontsize=15)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,2)
    ax1 = plt.contourf(pv_star_.T)
    ax2 = plt.contour(pv_star_.T , levels = [1e5], text='1e5')
    plt.clabel(ax2,inline=1)
    plt.colorbar(ax1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('Magnus Formula: '+r'$p_{v}^*$',fontsize=15)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,3)
    ax1 = plt.contourf(p0-pv_star_.T)
    ax2 = plt.contour(p0-pv_star_.T , levels = [0.0])
    plt.clabel(ax2,inline=1)
    #plt.contour(1e5)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('p0-pv_star',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,4)
    ax1 = plt.contourf(qv_star_.T)
    ax2 = plt.contour(p0-pv_star_.T,levels=[0.0],linewidth=2)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('qv_star',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.savefig('figures/'+str(plt_count)+'_CC_all_Magnus_'+name+'_guess.png')
    plt.close()
    plt_count += 1
    return

def plot_CC_all_2(T,pv_star_,qv_star_,pd_):
    global ns, plt_count
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    ax1 = plt.contourf(T.T,levels=np.linspace(np.amin(T),np.amax(T),500))
    print('!!! max(T)', np.amax(T))
    ax2 = plt.contour(p0-pv_star_.T,levels=[0.0],linewidths=3)
    ax3 = plt.contour(T.T,levels=[0.0,373],colors='w',linewidths=2)
    plt.clabel(ax2)
    plt.clabel(ax3,inline=1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.title(r'temperature $T_2$',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,2)
    # ax1 = plt.contourf(pv_star_.T)
    # ax2 = plt.contour(pv_star_.T , levels = [1e5], text='1e5')
    # plt.clabel(ax2,inline=1)
    # plt.colorbar(ax1)
    # plt.xlabel(r'$s$'+' entropy',fontsize=12)
    # plt.title('Magnus Formula: '+r'$p_{v,2}^*$',fontsize=15)
    # ax = plt.gca()
    # labels_x = ax.get_xticks().astype(int)
    # labels_y = ax.get_yticks()
    # lx, ly = set_ticks(labels_x,labels_y)
    # ax.set_xticklabels(lx)
    # ax.set_yticklabels(ly)
    # plt.subplot(1,4,3)
    ax1 = plt.contourf(p0-pv_star_.T)
    ax2 = plt.contour(p0-pv_star_.T , levels = [0.0])
    plt.clabel(ax2,inline=1)
    #plt.contour(1e5)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title(r'$p_0-p_{v,2}^*$ (Magnus Formula)',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    # plt.subplot(1,4,4)
    plt.subplot(1,4,3)
    ax1 = plt.contourf(qv_star_.T,levels=np.linspace(-10,10,250))
    ax2 = plt.contour(p0-pv_star_.T,levels=[0.0], linewidths=3)
    plt.clabel(ax2, inline=1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title(r'$q_{v,2}^*$',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,4)
    ax1 = plt.contourf(pd_.T)
    ax2 = plt.contour(p0 - pv_star_.T, levels=[0.0], linewidths=3)
    ax3 = plt.contour(pd_.T, levels=[0.0], colors='w', linewidths=1)
    plt.clabel(ax2, inline=1)
    plt.clabel(ax3, inline=1)
    plt.xlabel(r'$s$' + ' entropy', fontsize=12)
    plt.title(r'partial pressure: $p_{d,2}^*$', fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x, labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.savefig('figures/'+str(plt_count)+'_CC_all_Magnus_secondguess.png')
    plt.close()
    plt_count += 1
    return

def plot_CC_all(T,pv_star_1,qv_star_1):
    global ns, plt_count
    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1)
    ax1 = plt.contourf(T.T,levels=np.linspace(0,550,250))
    ax2 = plt.contour(p0-pv_star_1.T,levels=[0.0],colors='k',linewidths=2)
    ax3 = plt.contour(T.T,levels=np.linspace(0,500,6),colors='k',linestyle='--')
    plt.clabel(ax3,inline=1)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=12)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('temperature T1 (no ql)',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,2)
    plt.contourf(pv_star_1.T)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('Magnus Formula: '+r'$p_{v,1}^*$',fontsize=15)
    plt.colorbar()
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,3)
    ax1 = plt.contourf(p0-pv_star_1.T)#,levels=np.linspace(-9e6,1e6))
    ax2 = plt.contour(p0-pv_star_1.T, levels = [0.0])
    plt.clabel(ax2,inline=1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('p0-pv_star',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.subplot(1,4,4)
    ax1 = plt.contourf(qv_star_1.T,levels=np.linspace(-10,10,250))
    ax2 = plt.contour(qv_star_1.T,levels=[-1,0,1],linewidth=2,fontsize=6)
    plt.clabel(ax2,inline=1)
    plt.xlabel(r'$s$'+' entropy',fontsize=12)
    plt.title('qv_star_1',fontsize=15)
    plt.colorbar(ax1)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
    plt.savefig('figures/'+str(plt_count)+'_CC_all_Magnus_firstguess.png')
    plt.close()
    plt_count += 1
    return

def plot_pd(T,pd):
    global plt_count
    plt.figure()
    plt.plot(pd,linewidth=2)
    plt.ylabel(r'$p_d$',fontsize=21)
    plt.savefig('figures/1_pd.png')
    plt.close()
    plt_count += 1
    return

def plot_temp(T,name):
    global plt_count
    global ns, nqt
    plt.figure()
    plt.contourf(T.T)
    ax = plt.gca()
    labels_x = ax.get_xticks().astype(int)
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx)
    ax.set_yticklabels(ly)
#    lab_x = [s[0]]
#    lab_y = [qt[0]]
#    for i in range(1,labels_x.shape[0]):
#        if labels_x[i] < ns:
#            lab_x= np.append(lab_x,int(s[int(labels_x[i])]))
#    for i in range(1,labels_y.shape[0]):
#        if labels_y[i] < nqt:
#            lab_y = np.append(lab_y,np.round(qt[int(labels_y[i])],2))
#    ax.set_xticklabels(lab_x)
#    ax.set_yticklabels(lab_y)
    plt.xlabel(r'$s$'+' entropy',fontsize=18)
    plt.ylabel(r'$q_t$'+' moisture',fontsize=18)
    plt.title(name,fontsize=21)
    plt.colorbar()
    plt.savefig('figures/'+str(plt_count)+'_temp.png')
    plt_count += 1
    return


def plot_CC(T,pv_star):
    global plt_count
    plt.figure()
    plt.plot(T,pv_star,linewidth=2)
    plt.plot(T_ref,pv_star_ref,'--',linewidth=2)
    plt.plot([373,373],[min(pv_star),max(pv_star)])
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
def set_ticks(labels_x,labels_y):
    global plt_count
    lab_x = [s[0]]
    lab_y = [qt[0]]
    for i in range(1,labels_x.shape[0]):
        if labels_x[i] < ns:
            lab_x= np.append(lab_x,int(s[int(labels_x[i])]))
    lab_x = lab_x.astype(int)
    for i in range(1,labels_y.shape[0]):
        if labels_y[i] < nqt:
            lab_y = np.append(lab_y,np.round(qt[int(labels_y[i])],2))
    return lab_x, lab_y


# ------------------------------------------------
if __name__ == '__main__':
    main()

