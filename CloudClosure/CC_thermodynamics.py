import numpy as np
import sys

sys.path.append("..")
from parameters import *
sys.path.append("../Thermo/")



# # ---------------------------------------------------------------------------
# # def sat_adj_theta(p0, thetali, qt):
# #     '''
# #         (1) first guess: pv, pd, T
# #             - assumption: no liquid, i.e. qv=qt
# #     '''
# #     pv_1 = pv_c(p0, qt, qt)
# #     pd_1 = p0 - pv_1
# #     # Temperature for ql = 0: theta_l(ql=0) = theta
# #     T_1 = thetali * exner(p0)
# #
# #
# #
# #     ql = np.array(qt, copy=True)
# #     T = np.array(qt, copy=True)
# #     return T, ql
# # ---------------------------------------------------------------------------
# # def reference_state(zrange):
# #     from scipy.integrate import odeint
# #     print('Reference pressure profile initialization')
# #
# #     # sg = Thermodynamics.entropy(self.Pg, self.Tg, self.qtg, 0.0, 0.0)
# #
# #     # Form a right hand side for integrating the hydrostatic equation to
# #     # determine the reference pressure
# #     def rhs(p, z):
# #         # T, ql, qi = Thermodynamics.eos(np.exp(p), self.sg, self.qtg)
# #         return -g / (Rd * T * (1.0 - self.qtg + eps_vi * (self.qtg - ql - qi)))
# #
# #     # Construct arrays for integration points
# #     z = zrange
# #     nz = zrange.size
# #     # z = np.array(Gr.z[Gr.dims.gw - 1:-Gr.dims.gw + 1])
# #     # z_half = np.append([0.0], np.array(Gr.z_half[Gr.dims.gw:-Gr.dims.gw]))
# #
# #     # We are integrating the log pressure so need to take the log of the surface pressure
# #     # p0 = np.log(self.Pg)
# #     p0 = np.log(1e5)
# #
# #     p = np.zeros(nz, dtype=np.double, order='c')
# #     p_half = np.zeros(nz, dtype=np.double, order='c')
# #
# #     # Perform the integration
# #     p[:] = odeint(rhs, p0, z, hmax=1.0)[:, 0]
# #
# #     p = np.exp(p)
# #
# #     return
#
#
# ---------------------------------------------------------------------------

#----------------------------------------------------------------------
def sat_adj_fromentropy(p, s, qt):
    print ('')
    print('saturation adjustment from entropy')
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

    sys.path.append("../Thermo/")
    from thermodynamic_functions import pv_c, temperature_no_ql_from_entropy, CC_Magnus, qv_star_c, latent_heat
    from thermodynamic_functions import s_dry, s_vap, s_cond
    from parameters import cpd, cpv

    qv = qt.astype(np.double)
    ql = np.double(0.0)

    ''' (1) starting point (first guess): pv, pd, T (assumption: no liquid, i.e. qv=qt) '''
    pv_1 = pv_c(p,qt,qt)
    pd_1 = p - pv_1
    T_1 = temperature_no_ql_from_entropy(pd_1, pv_1, s, qt)
    pv_star_1 = CC_Magnus(T_1)                  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    qv_star_1 = qv_star_c(p, qt, pv_star_1)    # eos_c: same

    ''' (2) Check if saturated or not '''
    if (qt <= qv_star_1):
        T = T_1
        ql = 0.0
        print("not saturated: (s,qt)=", round(s,2), round(qt,4), round(T,1))
        alpha = 0
    else:
        ''' (3) if saturated: calculate second starting point T_2 '''
        print("air saturated: (s,qt)=", round(s,2), round(qt,4))
        alpha = 1
        ql_1 = qt - qv_star_1               # eos_c: same; thermo.py = ql_1
        # lam_1 = lam_fp(T_1);              # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid)
        L_1 = latent_heat(T_1)              # eos_c: L_1 = L_fp(T_1, lam_1)
        s_1 = s_dry(pd_1, T_1)*(1.0-qt) + s_vap(pv_1, T_1)*qt + s_cond(L_1,T_1)*ql_1
        f_1 = s - s_1
        T_2 = T_1 + ql_1*L_1 / ((1.0-qt)*cpd + qv_star_1*cpv)
        delta_T = np.abs(T_2 - T_1)

        T_2ndestimate = T_2
        count = 0
        pv_star_2 = CC_Magnus(T_2)  # pv_star_2 = lookup(LT, T_2)
        qv_star_2 = qv_star_c(p, qt, pv_star_2)
        ql_2 = qt - qv_star_2
        while(delta_T >= 1.0e-3 or ql_2 < 0.0):
            print('do loop:',  'ql2='+str(ql_2) + ', deltaT='+str(delta_T))
            pv_star_2 = CC_Magnus(T_2)      # pv_star_2 = lookup(LT, T_2)
            qv_star_2 = qv_star_c(p, qt, pv_star_2)
            pv_2 = pv_c(p, qt, qv_star_2)
            pd_2 = p - pv_2
            ql_2 = qt - qv_star_2
            if ql_2 < 0:
                print('ql_2 negative (count='+str(count)+')!', ql_2, 'delta T:', delta_T)
            if ql_2 < 0:
                print('ql_2 negative (count='+str(count)+')!', ql_2, 'delta T:', delta_T)

            L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
            s_2 = s_dry(pd_2,T_2) * (1.0 - qt) + s_vap(pv_2,T_2) * qt + s_cond(L_2,T_2)*ql_2
            f_2 = s - s_2
            T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)

            # if f_2 != 0:
            #     T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
            # else:
            #     T_n = T_2
            #     print('f_2 = 0: ', T_n, T_2, 'ql: ', ql_2)
            # # if np.isnan(T_n):
            #     print('!!!! f1: ', f_1, 'f2: ', f_2)
            #     break
            T_1 = T_2
            T_2 = T_n
            f_1 = f_2
            delta_T = np.abs(T_2 - T_1)
            count += 1
        T = T_2
        qv = qv_star_2
        ql = ql_2
        print("saturated: (s,qt)=", round(s,2), round(qt,4)," iterations = ",count)
    return T, ql, alpha


def sat_adj_fromthetali(p, thl, qt):
    print ('')
    print('saturation adjustment from thetali')
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

    sys.path.append("../Thermo/")
    from thermodynamic_functions import pv_c, temperature_no_ql_from_thetali, CC_Magnus, qv_star_c, latent_heat
    # from thermodynamic_functions import s_dry, s_vap, s_cond
    from thermodynamic_functions import theta_li
    from parameters import cpd, cpv

    qv = qt.astype(np.double)
    ql = np.double(0.0)

    ''' (1) starting point (first guess): pv, pd, T (assumption: no liquid, i.e. qv=qt) '''
    pv_1 = pv_c(p,qt,qt)
    pd_1 = p - pv_1
    # !!!!!!!
    T_1 = temperature_no_ql_from_thetali(pd_1, pv_1, thl, qt)
    # !!!!!!!
    pv_star_1 = CC_Magnus(T_1)                  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
    qv_star_1 = qv_star_c(p, qt, pv_star_1)    # eos_c: same
    print('T1: ', T_1)

    ''' (2) Check if saturated or not '''
    if (qt <= qv_star_1):
        T = T_1
        ql = 0.0
        print("not saturated: (thl,qt)=", round(thl,2), round(qt,4), round(T,1))
        alpha = 0
    else:
        ''' (3) if saturated: calculate second starting point T_2 '''
        print("air saturated: (thl,qt)=", round(thl,2), round(qt,4))
        alpha = 1
        ql_1 = qt - qv_star_1               # eos_c: same; thermo.py = ql_1
        # !!!!
        L_1 = latent_heat(T_1)  # eos_c: L_1 = L_fp(T_1, lam_1)
        thl_1 = theta_li(p,T_1,qt,ql_1,0)
        f_1 = thl - thl_1
        # s_1 = s_dry(pd_1, T_1) * (1.0 - qt) + s_vap(pv_1, T_1) * qt + s_cond(L_1, T_1) * ql_1
        # f_1 = s - s_1
        # !!!!
        T_2 = T_1 + ql_1 * L_1 / ((1.0 - qt) * cpd + qv_star_1 * cpv)
        delta_T = np.abs(T_2 - T_1)

        count = 0
        # pv_star_2 = CC_Magnus(T_2)  # pv_star_2 = lookup(LT, T_2)
        # qv_star_2 = qv_star_c(p, qt, pv_star_2)
        # ql_2 = qt - qv_star_2
        ql_2 = np.double(0.0)
        while(delta_T >= 1.0e-3 or ql_2 < 0.0):
            print('do loop: T2=' + str(T_2)+ ', ql2='+str(ql_2) + ', deltaT='+str(delta_T))
            pv_star_2 = CC_Magnus(T_2)      # pv_star_2 = lookup(LT, T_2)
            qv_star_2 = qv_star_c(p, qt, pv_star_2)
            pv_2 = pv_c(p, qt, qv_star_2)
            pd_2 = p - pv_2
            ql_2 = qt - qv_star_2
            if ql_2 < 0:
                print('ql_2 negative in satadj_thl (count='+str(count)+')!', ql_2, 'delta T:', delta_T)

            # !!!!!
            L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
            thl_2 = theta_li(p, T_2, qt, ql_2, 0)
            f_2 = thl - thl_2
            # s_2 = s_dry(pd_2,T_2)*(1.0-qt) + s_vap(pv_2,T_2)*qt + s_cond(L_2,T_2)*ql_2
            # f_2 = s - s_2
            # !!!!!
            T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
            T_1 = T_2
            T_2 = T_n
            f_1 = f_2
            delta_T = np.abs(T_2 - T_1)
            count += 1
        T = T_2
        qv = qv_star_2
        ql = ql_2
        print('count = ', count)

    return T, ql, alpha


