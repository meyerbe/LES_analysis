from parameters import *
import numpy as np
import pylab as plt


def main():
    # sat_adj modified for qt = 0:
    #   --> pv = 0 for qt=qv=0 (from pv_c(p0,qt,qv)
    #       --> temperature_no_ql(pv=0) ill-defined, due to T~log(pv/p_tilde)
    #           --> redefined in terms of (pv/p_tilde)**a

    ijk_ = 13841

    case_name = 'Bomex'
    path = './../Output.' + case_name + '.e0f33/'
    print(path)

    # def eos_c(lam_fp, L_fp, p0, s, qt):
    p0 = 1e5

    ns = 30
    nqt = 20


    T = np.zeros(shape=(ns, nqt))
    ql = np.zeros(shape=(ns, nqt))

    # QL
    s = np.linspace(6900,7060,ns)
    qt = np.linspace(0,0.045,nqt)
    name = 'sat_adj_ql'
    for i in range(ns):
        for j in range(nqt):
            T[i,j], ql[i,j], qi = sat_adj(p0, s[i], qt[j])
            if np.isnan(T[i,j]):
                print('T is nan for s,qt =', s[i], qt[j])
            if np.isnan(ql[i,j]):
                print('ql is nan for s,qt =', s[i], qt[j])
    # print('ns', ns, 'nqt', nqt)
    plot_ql(T,ql,s,qt,name)

    # name = 'sat_adj_T1_ql'
    # print(name)
    # for i in range(ns):
    #     for j in range(nqt):
    #         T[i, j], ql[i, j], qi = sat_adj_firstguess(p0, s[i], qt[j])
    #         if np.isnan(T[i, j]):
    #             print('T is nan for s,qt =', s[i], qt[j])
    #         if np.isnan(ql[i, j]):
    #             print('ql is nan for s,qt =', s[i], qt[j])
    # plot_ql(T, ql, s, qt, name)

    # full
    s = np.linspace(6900, 6970, ns)
    qt = np.linspace(0, 0.02, nqt)
    name = 'sat_adj_full'
    for i in range(ns):
        for j in range(nqt):
            T[i, j], ql[i, j], qi = sat_adj(p0, s[i], qt[j])
    plot_full(T, ql, s, qt, name)
    #
    # name = 'sat_adj_T1_full'
    # print(name)
    # for i in range(ns):
    #     for j in range(nqt):
    #         T[i, j], ql[i, j], qi = sat_adj_firstguess(p0, s[i], qt[j])
    #         if np.isnan(T[i, j]):
    #             print('T is nan for s,qt =', s[i], qt[j])
    #         if np.isnan(ql[i, j]):
    #             print('ql is nan for s,qt =', s[i], qt[j])
    # plot_full(T, ql, s, qt, name)

    # nqt = 2
    # qt = np.linspace(0.04,0.045,2)
    # s = 6980
    # for i in range(nqt):
    #     T[0,i], ql[0,i], qi = sat_adj(p0,s,qt[i])
    #     # print('qt: ', qt[i], 'ql: ', ql[0,i], 'T:', T[0,i])




    return



def plot_ql(T,ql,s,qt,name):
    levels_t = np.linspace(283,365,250)
    levels_ql = np.linspace(-300, 700, 250)
    plt.figure(figsize=(30,10))
    plt.subplot(1,3,1)
    plt.contourf(T,levels=levels_t)
    #plt.contourf(T)
    plt.title('temperature T (min/max: '+np.str(np.round(np.nanmin(T),1))+', '+np.str(np.round(np.nanmax(T),1))+')',fontsize=21)
    plt.xlabel('qt',fontsize=18)
    plt.ylabel('entropy s',fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    for i in range(labels.shape[0]-1):
        labels[i] = np.round(qt[labels[i]], 4)
    ax.set_xticklabels(labels)

    labels = ax.get_yticks()
    i = int(0)
    for i in range(0,labels.shape[0]-1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.subplot(1, 3, 2)
    plt.contourf(ql,levels=levels_ql)
    plt.title('ql',fontsize=24)
    plt.xlabel('qt',fontsize=18)
    plt.ylabel('entropy s',fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    for i in range(labels.shape[0]- 1):
        labels[i] = np.round(qt[labels[i]],4)
    ax.set_xticklabels(labels)
    labels = ax.get_yticks()
    for i in range(labels.shape[0]-1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.subplot(1, 3, 3)
    levels_ql = np.linspace(0, 700, 250)
    plt.contourf(ql, levels=levels_ql)
    plt.title('ql',fontsize=24)
    plt.xlabel('qt',fontsize=18)
    plt.ylabel('entropy s',fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    for i in range(labels.shape[0]-1):
        labels[i] = np.round(qt[labels[i]], 4)
    ax.set_xticklabels(labels)
    labels = ax.get_yticks()
    for i in range(labels.shape[0]-1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.savefig('./' + name + '.png')
    # plt.show()

    plt.close()

    return

def plot_full(T,ql,s,qt,name):
    #levels_t = np.linspace(270,330,250)
    levels_ql = np.linspace(0, 0.018, 250)
    plt.figure(figsize=(30,10))
    plt.subplot(1,3,1)
    #plt.contourf(T,levels=levels_t)
    plt.contourf(T)
    plt.title('temperature T',fontsize=24)
    plt.title('temperature T (min/max: '+np.str(np.round(np.nanmin(T),1))+', '+np.str(np.round(np.nanmax(T),1))+')',fontsize=21)
    plt.xlabel('qt',fontsize=18)
    plt.ylabel('entropy s',fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    i = int(0)
    for i in range(labels.shape[0]- 1):
        labels[i] = np.round(qt[labels[i]], 4)
    ax.set_xticklabels(labels)
    labels = ax.get_yticks()
    for i in range(labels.shape[0]- 1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.subplot(1, 3, 2)
    plt.contourf(ql,levels=levels_ql)
    plt.title('ql',fontsize=24)
    plt.xlabel('qt',fontsize=18)
    plt.ylabel('entropy s',fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    for i in range(labels.shape[0] - 1):
        labels[i] = np.round(qt[labels[i]], 4)
    ax.set_xticklabels(labels)
    labels = ax.get_yticks()
    for i in range(labels.shape[0] - 1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.subplot(1, 3, 3)
    levels_ql = np.linspace(0, 0.008, 250)
    plt.contourf(ql, levels=levels_ql)
    plt.title('ql', fontsize=24)
    plt.xlabel('qt', fontsize=18)
    plt.ylabel('entropy s', fontsize=18)
    plt.colorbar()

    ax = plt.gca()
    ax.tick_params(direction='out', pad=0)
    labels = ax.get_xticks()
    i = int(0)
    for i in range(labels.shape[0] - 1):
        labels[i] = np.round(qt[labels[i]], 4)
    ax.set_xticklabels(labels)
    labels = ax.get_yticks()
    for i in range(labels.shape[0] - 1):
        labels[i] = np.round(s[labels[i]], 2)
    ax.set_yticklabels(labels)

    plt.savefig('./' + name + '.png')
    # plt.show()

    plt.close()
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
    T_1 = temperature_no_ql(pd_1,pv_1,s,qt)     # eos_c: same
    #Compute saturation vapor pressure
    pv_star_1 = get_pv_star(T_1)                # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    #Compute saturation mixing ratio
    qv_star_1 = qv_star_c(p0, qt, pv_star_1)    # eos_c: same

     # If not saturated
    if(qt <= qv_star_1):
        T = T_1
        ql = 0.0
        qi = 0.0
        # print("not saturated: no iteration");
        if np.isnan(T):
            print('T is nan')
        return T, ql, qi
    else:
        # print("saturated: start iterations");
        sigma_1 = qt - qv_star_1           # eos_c: same; thermo.py = ql_1
        lam_1 = lam_fp(T_1);                # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid)
        L_1 = latent_heat(T_1)              # eos_c: L_1 = L_fp(T_1, lam_1)
        # print('T_1, pd_1', T_1, pd_1)
        s_1 = sd_c(pd_1,T_1)*(1.0 - qt) + sv_c(pv_1,T_1)*qt + sc_c(L_1,T_1)*sigma_1
        f_1 = s - s_1
        T_2 = T_1 + sigma_1 * L_1 /((1.0 - qt)*cpd + qv_star_1 * cpv)
        # thermo.py: t_2 = t_1 + ql_1 * lv / cpd
        # T_2 = T_1 + sigma_1 * L(T_1) / cp
        # cp = (1.0 - qt)*cpd + qv_star_1 * cpv
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
        print("saturated: iterations = ",count)
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
        # print("not saturated: no iteration");
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
# from entropies.h
def sd_c(pd, T):
    # print('pd', pd, 'T', T)
    return sd_tilde + cpd*np.log(T/T_tilde) - Rd*np.log(pd/p_tilde)

def sv_c(pv, T):
    return sv_tilde + cpv*np.log(T/T_tilde) - Rv * np.log(pv/p_tilde)

def sc_c(L, T):
    return -L/T


def pv_c(p0, qt, qv):
    return p0 * eps_vi * qv /(1.0 - qt + eps_vi * qv)

def qv_star_c(p0, qt, pv):
    return eps_v * (1.0 - qt) * pv / (p0 - pv)

# from Microphysics.pxd - lambda_constant(T)
def lam_fp(T):
    return 1.0

# from Microphysics.pxd - L_fp(T,lambda)
def latent_heat(T):
    TC = T - 273.15
    return (2500.8 - 2.36 * TC + 0.0016 * TC *
            TC - 0.00006 * TC * TC * TC) * 1000.0

def get_pv_star(T):
#    Magnus formula
    T = T - 273.15
    pv_star = 6.1094*np.exp((17.625*T)/float(T+243.04))*100
    return pv_star

#from thermodynamics_sa.h
def temperature_no_ql(pd, pv, s, qt):
    # print('temperature_no_ql')
    if pv > 0:
        # print('pv>0')
        temp = T_tilde * np.exp(  (s - (1.0-qt)*(sd_tilde - Rd * np.log(pd/p_tilde))
                            - qt * (sv_tilde - Rv * np.log(pv/p_tilde)))
                            /((1.0-qt)*cpd + qt * cpv)  )
    elif pv == 0:
        # print('pv=0')
        cp = ((1.0 - qt) * cpd + qt * cpv)
        temp = T_tilde * np.exp(  (s - (1.0-qt)*(sd_tilde - Rd * np.log(pd/p_tilde)) ) / cp) \
               * np.exp( -qt*sv_tilde / cp ) \
               * (pv/p_tilde)**(-qt*Rv/cp)
    else:
        print('pv < 0: not physical!!')
    return temp


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    main()