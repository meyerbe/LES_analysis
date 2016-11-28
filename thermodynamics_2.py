
''' offline thermodynamic computations as in LES '''


from parameters import *
import numpy as np
import pylab as plt


def main():
    # sat_adj modified for qt = 0:
    #   --> pv = 0 for qt=qv=0 (from pv_c(p0,qt,qv)
    #       --> temperature_no_ql(pv=0) ill-defined, due to T~log(pv/p_tilde)
    #           --> redefined in terms of (pv/p_tilde)**a

    # def eos_c(lam_fp, L_fp, p0, s, qt):
    p0 = 1e5

    global ns, nqt, s, qt
    ns = 250
    nqt = 250

    s = np.linspace(6800,7200,ns)
    qt = np.linspace(0,0.045,nqt)
    T = np.zeros(shape=(ns, nqt))
    ql = np.zeros(shape=(ns, nqt))

    name = 'sat_adj_'
    for i in range(ns):
        for j in range(nqt):
            T[i,j], ql[i,j], qi = sat_adj(p0, s[i], qt[j])
            if np.isnan(T[i,j]):
                print('T is nan for s,qt =', s[i], qt[j])
            if np.isnan(ql[i,j]):
                print('ql is nan for s,qt =', s[i], qt[j])
    # print('ns', ns, 'nqt', nqt)
    plot(T,ql,s,qt,name+'corr_qvstar')

    return


# ---------------------------------------------------------------------------

# def eos_c(lam_fp, L_fp, p0, s, qt):
def sat_adj(p0, s, qt):
    '''
    Use saturation adjustment scheme to compute temperature and ql given s and qt.
    :param p0: pressure [Pa]
    :param s: entropy  [K]
    :param qt:  total water specific humidity
    :return: T, ql, qi

    Functions from Microphyiscs.pxd:
        liquid fraction:
            lam_fp(T) = 1.0 (lambda_constant)
        Latent Heat:
            L_fp(T, lambda(T)) = (2500.8 - 2.36 * TC + 0.0016 * TC**2 - 0.00006 * TC**3) * 1000.0
                with: TC = T - 273.15

    Definitions (c.f. Pressel et al., 2016):
        saturation vapor pressure: pv_star(T)
            --> in LES from Clausius Clapeyron (Lookup table for integration), here use Magnus formula
        saturation specific humidity: qv_star(p,qt,pv_star)
            --> ideal gas law
        saturation excess: sigma = qt - qv_star
        partial entropies: sd_c(pd,T), sv_c(pv,T), sc_c(L,T)
            --> defined in Csrc/entropies.h
    '''

    # Compute temperature
    pv_1 = pv_c(p0,qt,qt)                       # eos_c: same
    pd_1 = p0 - pv_1                            # eos_c: same
    print('p0, pv_1, pd_1', p0, pv_1, pd_1)
    T_1 = temperature_no_ql(pd_1,pv_1,s,qt)     # eos_c: same
    #Compute saturation vapor pressure
    pv_star_1 = get_pv_star(T_1)                # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    #Compute saturation mixing ratio
    qv_star_1 = qv_star_c(p0, qt, pv_star_1)    # eos_c: same

     # If not saturated
    if(qt <= qv_star_1 or qv_star_1 < 0.0):
        T = T_1
        ql = 0.0
        qi = 0.0
        if np.isnan(T):
            print('T is nan')
        print("not saturated: (s,qt)=", round(s,2), round(qt,4))
        return T, ql, qi
    else:
        # print("saturated: start iterations");
        sigma_1 = qt - qv_star_1           # eos_c: same; thermo.py = ql_1
        lam_1 = lam_fp(T_1);                # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid); here lam_fp = 1.0
        L_1 = latent_heat(T_1)              # eos_c: L_1 = L_fp(T_1, lam_1)
        s_1 = sd_c(pd_1,T_1)*(1.0 - qt) + sv_c(pv_1,T_1)*qt + sc_c(L_1,T_1)*sigma_1
        f_1 = s - s_1
        T_2 = T_1 + sigma_1 * L_1 /((1.0 - qt)*cpd + qv_star_1 * cpv)
        delta_T = np.fabs(T_2 - T_1)

        count = 0

        pv_star_2 = get_pv_star(T_2)  # pv_star_2 = lookup(LT, T_2)
        qv_star_2 = qv_star_c(p0, qt, pv_star_2)
        ql_2 = qt - qv_star_2
        lam_2 = lam_fp(T_2)
        # while((delta_T >= 1.0e-3 or ql_2 < 0.0) and count < 2):
        while (delta_T >= 1.0e-3 or ql_2 < 0.0):
            pv_star_2 = get_pv_star(T_2)    # pv_star_2 = lookup(LT, T_2)
            qv_star_2 = qv_star_c(p0, qt, pv_star_2)
            ql_2 = qt - qv_star_2

            pv_2 = pv_c(p0,qt,qv_star_2)
            pd_2 = p0 - pv_2
            # if qv_star_2 < 0:
            #     print('')
            #     print('!!! qv_star < 0')
            # if pd_2 < 0:
            #     print('p0,pv_2', p0, pv_2, qt, qv_star_2, pv_star_2)
            #     print('pv_2, pv_star_2', pv_2, pv_star_2)
            #     print('T_2', T_2)

            lam_2 = lam_fp(T_2)         # eos_c
            L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
            # print('loop: T_2, pd_2', T_2, pd_2)
            # if np.isnan(T_2):
            #     print('')
            #     print('T_2 is nan')
            s_2 = sd_c(pd_2,T_2) * (1.0 - qt) + sv_c(pv_2,T_2) * qt + sc_c(L_2,T_2)*ql_2
            # if np.isnan(s_2):
            #     print('')
            #     print('s_2 is nan')
            #     print('pv_2', pv_2)     # pv_2 > p0 --> pd_2 < 0
            #     print('pd_2', pd_2)     # pd_2 < 0
            #     print('qv_star_2', qv_star_2)
            #     print('pv_star_2', pv_star_2)

            f_2 = s - s_2
            T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
            T_1 = T_2
            T_2 = T_n
            f_1 = f_2
            delta_T  = np.abs(T_2 - T_1)
            count += 1
        T  = T_2
        qv = qv_star_2
        ql = lam_2 * ql_2           # ??? lam_2
        qi = (1.0 - lam_2) * ql_2   # ??? lam_2
        print("saturated: (s,qt)=", round(s,2), round(qt,4)," iterations = ",count)
        # if np.isnan(T):
        #     print('T is nan', pd_1, pd_2)
        # if np.isnan(ql):
        #     print('ql is nan')
        return T, ql, qi

# ---------------------------------------------------------------------------

    # def eos_c(lam_fp, L_fp, p0, s, qt):
def sat_adj_lin(p0, s, qt):
    pv_1 = pv_c(p0, qt, qt)  # eos_c: same
    pd_1 = p0 - pv_1  # eos_c: same
    T_1 = temperature_no_ql(pd_1, pv_1, s, qt)  # eos_c: same
    pv_star_1 = get_pv_star(T_1)  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    qv_star_1 = qv_star_c(p0, qt, pv_star_1)  # eos_c: same

    if (qt <= qv_star_1):
        T = T_1
        ql = 0.0
        qi = 0.0
        # print("not saturated: no iteration");
        if np.isnan(T):
            print('T is nan')
        return T, ql, qi
    else:
        # print("saturated: start iterations");
        sigma_1 = qt - qv_star_1  # eos_c: same; thermo.py = ql_1
        lam_1 = lam_fp(
            T_1);  # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid)
        L_1 = latent_heat(T_1)  # eos_c: L_1 = L_fp(T_1, lam_1)
        # print('T_1, pd_1', T_1, pd_1)
        s_1 = sd_c(pd_1, T_1) * (1.0 - qt) + sv_c(pv_1, T_1) * qt + sc_c(L_1, T_1) * sigma_1
        f_1 = s - s_1
        T_2 = T_1 + sigma_1 * L_1 / ((1.0 - qt) * cpd + qv_star_1 * cpv)
        # thermo.py: t_2 = t_1 + ql_1 * lv / cpd
        # T_2 = T_1 + sigma_1 * L(T_1) / cp
        # cp = (1.0 - qt)*cpd + qv_star_1 * cpv
        delta_T = np.fabs(T_2 - T_1)

        count = 0

        pv_star_2 = get_pv_star(T_2)  # pv_star_2 = lookup(LT, T_2)
        qv_star_2 = qv_star_c(p0, qt, pv_star_2)
        ql_2 = qt - qv_star_2
        lam_2 = lam_fp(T_2)

        T = T_2
        qv = qv_star_2
        ql = lam_2 * ql_2  # ??? lam_2
        qi = (1.0 - lam_2) * ql_2  # ??? lam_2
        # print("saturated: iterations = %d\n",count)
        # if np.isnan(T):
        #     print('T is nan', pd_1, pd_2)
        # if np.isnan(ql):
        #     print('ql is nan')
        return T, ql, qi

# ---------------------------------------------------------------------------
def sat_adj_firstguess(p0, s, qt):
    pv_1 = pv_c(p0, qt, qt)  # eos_c: same
    pd_1 = p0 - pv_1  # eos_c: same
    T_1 = temperature_no_ql(pd_1, pv_1, s, qt)  # eos_c: same
    pv_star_1 = get_pv_star(T_1)  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    qv_star_1 = qv_star_c(p0, qt, pv_star_1)  # eos_c: same

    if (qt <= qv_star_1):
        T = T_1
        ql = 0.0
        qi = 0.0
        print("not saturated: no iteration");
        if np.isnan(T):
            print('T is nan')
        return T, ql, qi
    else:
        sigma_1 = qt - qv_star_1  # eos_c: same; thermo.py = ql_1
        lam_1 = lam_fp(T_1);  # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid)
        T = T_1
        ql = lam_1 * sigma_1
        qi = (1.0-lam_1) * sigma_1

        return T, ql, qi

# ---------------------------------------------------------------------------
''' entropies.h '''
def sd_c(pd, T):
    print('pd', pd, 'T', T)
    if T<0 or pd<0:
        print('sd: negative pd or T', pd, T)
    return sd_tilde + cpd*np.log(T/T_tilde) - Rd*np.log(pd/p_tilde)

def sv_c(pv, T):
    return sv_tilde + cpv*np.log(T/T_tilde) - Rv * np.log(pv/p_tilde)

def sc_c(L, T):
    return -L/T


def pv_c(p0, qt, qv):
    return p0 * eps_vi * qv /(1.0 - qt + eps_vi * qv)

def qv_star_c(p0, qt, pv):
    return eps_v * (1.0 - qt) * pv / (p0 - pv)

''' Microphysics.pxd - lambda_constant(T) '''
def lam_fp(T):
    return 1.0

''' Microphysics.pxd - L_fp(T,lambda) '''
def latent_heat(T):
    TC = T - 273.15
    return (2500.8 - 2.36 * TC + 0.0016 * TC *
            TC - 0.00006 * TC * TC * TC) * 1000.0

'''  Magnus formula '''
def get_pv_star(T):
    T = T - 273.15
    pv_star = 6.1094*np.exp((17.625*T)/float(T+243.04))*100
    return pv_star

''' thermodynamics_sa.h '''
def temperature_no_ql(pd, pv, s, qt):
    # print('temperature_no_ql')
    if pv > 0:
        temp = T_tilde * np.exp(  (s - (1.0-qt)*(sd_tilde - Rd * np.log(pd/p_tilde))
                            - qt * (sv_tilde - Rv * np.log(pv/p_tilde)))
                            /((1.0-qt)*cpd + qt * cpv)  )
    elif pv == 0:
        cp = ((1.0 - qt) * cpd + qt * cpv)
        temp = T_tilde * np.exp(  (s - (1.0-qt)*(sd_tilde - Rd * np.log(pd/p_tilde)) ) / cp) \
               * np.exp( -qt*sv_tilde / cp ) \
               * (pv/p_tilde)**(-qt*Rv/cp)
    else:
        print('pv < 0: not physical!!')
    return temp


# ---------------------------------------------------------------------------

def plot(T,ql,s,qt,name):
    levels_t = np.linspace(270,420,250)
    plt.figure(figsize=(25,8))
    plt.subplot(1,3,1)
    i_ = 0
    for i in xrange(ns):
        if s[i] < 7060:
            i_ = i
    plt.plot([i_,i_],[0,nqt],'k',linewidth=2)
    plt.contourf(T.T,levels=levels_t)

    plt.title('temperature T (min/max: '+np.str(np.round(np.nanmin(T),1))+', '+np.str(np.round(np.nanmax(T),1))+')',fontsize=18)
    plt.xlabel('entropy s',fontsize=15)
    plt.ylabel('qt',fontsize=15)
    plt.colorbar()
    
    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels_x = ax.get_xticks()
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx,fontsize=9)
    ax.set_yticklabels(ly,fontsize=9)
    
    plt.subplot(1, 3, 2)
    levels_ql = np.linspace(0, np.amax(qt), 250)
    plt.contourf(ql.T,levels=levels_ql)
    ax2 = plt.contour(ql.T,levels=np.append([0.0,0.01,0.02,0.03],np.linspace(0.2,1.0,5)),colors='w')
    plt.clabel(ax2,inline=1)
    plt.title('ql',fontsize=24)
    plt.xlabel('entropy s',fontsize=18)
    plt.ylabel('qt',fontsize=18)
    plt.title('ql (min/max: '+np.str(np.round(np.nanmin(ql),1))+', '+np.str(np.round(np.nanmax(ql),1))+')',fontsize=18)
    plt.xlabel('entropy s',fontsize=18)
    plt.ylabel('qt',fontsize=18)
    plt.colorbar()
    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels_x = ax.get_xticks()
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx,fontsize=9)
    ax.set_yticklabels(ly,fontsize=9)

    plt.subplot(1, 3, 3)
    levels_ql = np.linspace(0, 50, 250)
    ax1 = plt.contourf(ql.T, levels=levels_ql)
    ax2 = plt.contour(ql.T,levels=np.append([0.0,0.01,0.02,0.03],np.linspace(0.2,1.0,5)),colors='w')
    plt.clabel(ax2,inline=1)
    plt.title('ql',fontsize=24)
    plt.xlabel('entropy s',fontsize=18)
    plt.ylabel('qt',fontsize=18)
    plt.colorbar(ax1)
    
    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels_x = ax.get_xticks()
    labels_y = ax.get_yticks()
    lx, ly = set_ticks(labels_x,labels_y)
    ax.set_xticklabels(lx,fontsize=9)
    ax.set_yticklabels(ly,fontsize=9)
    
    plt.savefig('./figures/' + name + '.png')
    # plt.show()
    
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
            lab_y = np.append(lab_y,np.round(qt[int(labels_y[i])],3))
    return lab_x, lab_y

if __name__ == '__main__':
    main()
