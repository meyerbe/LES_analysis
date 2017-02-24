import numpy as np
cimport numpy as np
import sys
import cython

from scipy.integrate import odeint

sys.path.append("..")
# sys.path.append("../Thermo/")
include '../parameters.pxi'

# __________________________________________________________________
cdef class Lookup:
    def __init__(self):
        print('initialize cpdef Lookup Table')
        return

    cpdef initialize(self, double[:] x, double[:] y):
        cdef long n = np.shape(y)[0]
        print('LT_c.initialize')
        init_table( &self.LookupStructC, n, &x[0], &y[0])
        return

    cpdef finalize(self):
        free_table(&self.LookupStructC)
        return

    cpdef table_bounds(self):
        return self.LookupStructC.x_min, self.LookupStructC.x_max, self.LookupStructC.y_min, self.LookupStructC.y_max

    cpdef lookup(self, double x):
        return lookup(&self.LookupStructC, x)

    cdef inline double fast_lookup(self, double x) nogil:
        return lookup(&self.LookupStructC, x)
# __________________________________________________________________
cdef class ClausiusClapeyron_c:
    def __init__(self):
        print('CC c')
        return

    cpdef rhs(self, double z, double T_):
        # sys.path.append('../Thermo/')
        from thermodynamic_functions import latent_heat
        # lam = LH.Lambda(T_)
        # L = LH.L(T_, lam)
        lam = 1.0
        L = latent_heat(T_)
        return L / (Rv * T_ * T_)

    cpdef initialize(self):
        cdef:
            double Tmin, Tmax
            long n_lookup
            double [:] pv

        Tmin = 100.15
        Tmax = 380.0
        n_lookup = 512

        # Generate array of equally space temperatures
        T = np.linspace(Tmin, Tmax, n_lookup)
        # Find the maximum index of T where T < T_tilde
        tp_close_index = np.max(np.where(T <= Tt))

        # Check to make sure that T_tilde is not in T
        if T[tp_close_index] == Tt:
            print('Array of temperatures for ClasiusClapyeron lookup table contains Tt  \n')
            print('Pick different values for ClasiusClapyeron Tmin and Tmax in lookup table \n')
            print('Killing Simulation now!')
            sys.exit()

        # Now prepare integration
        T_above_Tt = np.append([Tt], T[tp_close_index + 1:])
        T_below_Tt = np.append(T[:tp_close_index + 1], [Tt])[::-1]

        # Set the initial condition
        pv0 = np.log(pv_star_t)

        # Integrate
        pv_above_Tt = np.exp(odeint(self.rhs, pv0, T_above_Tt, hmax=0.1)[1:])
        pv_below_Tt = np.exp(odeint(self.rhs, pv0, T_below_Tt, hmax=0.1)[1:])[::-1]
        pv = np.append(pv_below_Tt, pv_above_Tt)
        # self.LT.initialize(T, pv)

        self.LT = Lookup()
        self.LT.initialize(T,pv)
        print(self.LT.lookup(298.0))
        print(self.LT.lookup(296.0))
        # print('pv: ', pv)

        # LT_c = LookupTable_c()
        # LT_c.make_lookup_table(T,pv)
        # print('given value of T: ', LT_c.lookup_table[T[0]])        # only possible to lookup for values of T that are already calculated --> need interpolation
        # print('interpolated value of T: ', LT_c.lookup(T, 298.0))
        return
# __________________________________________________________________
cdef class LatentHeat:
    # def __init__(self,namelist):
    def __init__(self):
        return

    cpdef L(self,double T, double Lambda):
        '''
        Provide a python interface to the latent heat function pointer.
        :param T (Thermodynamic Temperature):
        :return L (Latent Heat):
        '''
        return self.L_fp(T, Lambda)

    cpdef Lambda(self, double T):
        return self.Lambda_fp(T)
# __________________________________________________________________
cdef extern from "thermodynamics_sa.h":
# cdef extern from "../Thermo/thermodynamics_sa.h":
    inline double alpha_c(double p0, double T, double qt, double qv) nogil
    # void eos_c(Lookup.LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double p0, double s, double qt, double *T, double *qv, double *ql, double *qi) nogil
    void eos_c(LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double p0, double s, double qt, double *T, double *qv, double *ql, double *qi) nogil
#     void eos_update(Grid.DimStruct *dims, Lookup.LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double *p0, double *s, double *qt, double *T,
#                     double * qv, double * ql, double * qi, double * alpha)
#     void buoyancy_update_sa(Grid.DimStruct *dims, double *alpha0, double *alpha, double *buoyancy, double *wt)
#     void bvf_sa(Grid.DimStruct * dims, Lookup.LookupStruct * LT, double(*lam_fp)(double), double(*L_fp)(double, double), double *p0, double *T, double *qt, double *qv, double *theta_rho, double *bvf)
#     void thetali_update(Grid.DimStruct *dims, double (*lam_fp)(double), double (*L_fp)(double, double), double *p0, double *T, double *qt, double *ql, double *qi, double *thetali)
#     void clip_qt(Grid.DimStruct *dims, double  *qt, double clip_value)

cdef extern from "thermodynamic_functions.h":
# cdef extern from "../Thermo/thermodynamic_functions.h":
    # Dry air partial pressure
    inline double pd_c(double p0, double qt, double qv) nogil
    # Water vapor partial pressure
    inline double pv_c(double p0, double qt, double qv) nogil

# cpdef eos(self, double p0, double s, double qt):
cpdef sat_adj_fromentropy_c(double p0, double s, double qt, microphysics):
    print('calling here')
    cdef:
        double T, qv, qc, ql, qi, lam
        int alpha = 0

    CC = ClausiusClapeyron_c()
    # CC.initialize(namelist, LH, Par)
    CC.initialize()

    LH = LatentHeat()
    if microphysics == 'dry':
        LH.Lambda_fp = lambda_constant
        LH.L_fp = latent_heat_constant
    elif microphysics == 'sa':
        LH.Lambda_fp = lambda_constant
        LH.L_fp = latent_heat_variable
    L_fp = LH.L_fp
    Lambda_fp = LH.Lambda_fp
    print('ok, lets try it')
    # eos_c(&self.CC.LT.LookupStructC, self.Lambda_fp, self.L_fp, p0, s, qt, &T, &qv, &ql, &qi)
    eos_c(&CC.LT.LookupStructC, Lambda_fp, L_fp, p0, s, qt, &T, &qv, &ql, &qi)
    print('done')
    if ql > 0.0:
        print('saturated')
        alpha = 1
    else:
        print('dry')
    # return T, ql, qi
    return T, ql, alpha


# cpdef sat_adj_fromentropy(p, s, qt):
#     print ('')
#     print('saturation adjustment from entropy')
#     '''
#     Use saturation adjustment scheme to compute temperature and ql given s and qt.
#     :param p0: pressure [Pa]
#     :param s: entropy  [K]
#     :param qt:  total water specific humidity
#     :return: T, ql, qi
#
#     Functions from Microphyiscs.pxd:
#         liquid fraction:
#             lam_fp(T) = 1.0 (lambda_constant)
#         Latent Heat:
#             L_fp(T, lambda(T)) = (2500.8 - 2.36 * TC + 0.0016 * TC**2 - 0.00006 * TC**3) * 1000.0
#                 with: TC = T - 273.15
#
#     Definitions (c.f. Pressel et al., 2016):
#         saturation vapor pressure: pv_star(T)
#             --> in LES from Clausius Clapeyron (Lookup table for integration), here use Magnus formula
#         saturation specific humidity: qv_star(p,qt,pv_star)
#             --> ideal gas law
#         saturation excess: sigma = qt - qv_star
#         partial entropies: sd_c(pd,T), sv_c(pv,T), sc_c(L,T)
#             --> defined in Csrc/entropies.h
#     '''
#
#     sys.path.append("../Thermo/")
#     from thermodynamic_functions import pv_c, temperature_no_ql_from_entropy, CC_Magnus, qv_star_c, latent_heat
#     from thermodynamic_functions import s_dry, s_vap, s_cond
#     from parameters import cpd, cpv
#
#     cdef:
#         double T, qv
#         double ql = 0.0
#         double pv_1, pd_1, pv_star_1, qv_star_1, T_1
#         double pv_2, pd_2, pv_star_2, qv_star_2, T_2
#         double T_n
#
#     ''' (1) starting point (first guess): pv, pd, T (assumption: no liquid, i.e. qv=qt) '''
#     pv_1 = pv_c(p,qt,qt)
#     pd_1 = p - pv_1
#     T_1 = temperature_no_ql_from_entropy(pd_1, pv_1, s, qt)
#     pv_star_1 = CC_Magnus(T_1)                  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
#     qv_star_1 = qv_star_c(p, qt, pv_star_1)    # eos_c: same
#
#     ''' (2) Check if saturated or not '''
#     if (qt <= qv_star_1):
#         T = T_1
#         ql = 0.0
#         print("not saturated: (s,qt)=", round(s,2), round(qt,4), round(T,1))
#         alpha = 0
#     else:
#         ''' (3) if saturated: calculate second starting point T_2 '''
#         print("air saturated: (s,qt)=", round(s,2), round(qt,4))
#         alpha = 1
#         ql_1 = qt - qv_star_1               # eos_c: same; thermo.py = ql_1
#         # lam_1 = lam_fp(T_1);              # eos_c: lam_fp gives the liquid fraction for mixed - phase clouds(fraction of supercooled liquid)
#         L_1 = latent_heat(T_1)              # eos_c: L_1 = L_fp(T_1, lam_1)
#         s_1 = s_dry(pd_1, T_1)*(1.0-qt) + s_vap(pv_1, T_1)*qt + s_cond(L_1,T_1)*ql_1
#         f_1 = s - s_1
#         T_2 = T_1 + ql_1*L_1 / ((1.0-qt)*cpd + qv_star_1*cpv)
#         delta_T = np.abs(T_2 - T_1)
#
#         T_2ndestimate = T_2
#         count = 0
#         pv_star_2 = CC_Magnus(T_2)  # pv_star_2 = lookup(LT, T_2)
#         qv_star_2 = qv_star_c(p, qt, pv_star_2)
#         ql_2 = qt - qv_star_2
#         while(delta_T >= 1.0e-3 or ql_2 < 0.0):
#             print('do loop:',  'ql2='+str(ql_2) + ', deltaT='+str(delta_T))
#             pv_star_2 = CC_Magnus(T_2)      # pv_star_2 = lookup(LT, T_2)
#             qv_star_2 = qv_star_c(p, qt, pv_star_2)
#             pv_2 = pv_c(p, qt, qv_star_2)
#             pd_2 = p - pv_2
#             ql_2 = qt - qv_star_2
#             if ql_2 < 0:
#                 print('ql_2 negative (count='+str(count)+')!', ql_2, 'delta T:', delta_T)
#             if ql_2 < 0:
#                 print('ql_2 negative (count='+str(count)+')!', ql_2, 'delta T:', delta_T)
#
#             L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
#             s_2 = s_dry(pd_2,T_2) * (1.0 - qt) + s_vap(pv_2,T_2) * qt + s_cond(L_2,T_2)*ql_2
#             f_2 = s - s_2
#             T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
#
#             # if f_2 != 0:
#             #     T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
#             # else:
#             #     T_n = T_2
#             #     print('f_2 = 0: ', T_n, T_2, 'ql: ', ql_2)
#             # # if np.isnan(T_n):
#             #     print('!!!! f1: ', f_1, 'f2: ', f_2)
#             #     break
#             T_1 = T_2
#             T_2 = T_n
#             f_1 = f_2
#             delta_T = np.abs(T_2 - T_1)
#             count += 1
#         T = T_2
#         qv = qv_star_2
#         ql = ql_2
#         print("saturated: (s,qt)=", round(s,2), round(qt,4)," iterations = ",count)
#     return T, ql, alpha


def sat_adj_fromthetali_c(p, thl, qt):
#     print ('')
#     print('saturation adjustment from thetali')
#     '''
#     Use saturation adjustment scheme to compute temperature and ql given s and qt.
#     :param p0: pressure [Pa]
#     :param s: entropy  [K]
#     :param qt:  total water specific humidity
#     :return: T, ql, qi
#
#     Functions from Microphyiscs.pxd:
#         liquid fraction:
#             lam_fp(T) = 1.0 (lambda_constant)
#         Latent Heat:
#             L_fp(T, lambda(T)) = (2500.8 - 2.36 * TC + 0.0016 * TC**2 - 0.00006 * TC**3) * 1000.0
#                 with: TC = T - 273.15
#
#     Definitions (c.f. Pressel et al., 2016):
#         saturation vapor pressure: pv_star(T)
#             --> in LES from Clausius Clapeyron (Lookup table for integration), here use Magnus formula
#         saturation specific humidity: qv_star(p,qt,pv_star)
#             --> ideal gas law
#         saturation excess: sigma = qt - qv_star
#         partial entropies: sd_c(pd,T), sv_c(pv,T), sc_c(L,T)
#             --> defined in Csrc/entropies.h
#     '''
#
#     sys.path.append("../Thermo/")
#     from thermodynamic_functions import pv_c, temperature_no_ql_from_thetali, CC_Magnus, qv_star_c, latent_heat
#     # from thermodynamic_functions import s_dry, s_vap, s_cond
#     from thermodynamic_functions import theta_li
#     from parameters import cpd, cpv
#
#     qv = qt.astype(np.double)
#     ql = np.double(0.0)
#
#     ''' (1) starting point (first guess): pv, pd, T (assumption: no liquid, i.e. qv=qt) '''
#     pv_1 = pv_c(p,qt,qt)
#     pd_1 = p - pv_1
#     # !!!!!!!
#     T_1 = temperature_no_ql_from_thetali(pd_1, pv_1, thl, qt)
#     # !!!!!!!
#     pv_star_1 = CC_Magnus(T_1)                  # eos_c: pv_star_1 = lookup(LT, T_1) # ???
#     qv_star_1 = qv_star_c(p, qt, pv_star_1)    # eos_c: same
#
#     ''' (2) Check if saturated or not '''
#     if (qt <= qv_star_1):
#         T = T_1
#         ql = 0.0
#         # print("not saturated: (thl,qt)=", round(thl,2), round(qt,4), round(T,1))
#         alpha = 0
#     else:
#         ''' (3) if saturated: calculate second starting point T_2 '''
#         # print("air saturated: (thl,qt)=", round(thl,2), round(qt,4))
#         alpha = 1
#         ql_1 = qt - qv_star_1               # eos_c: same; thermo.py = ql_1
#         # !!!!
#         L_1 = latent_heat(T_1)  # eos_c: L_1 = L_fp(T_1, lam_1)
#         thl_1 = theta_li(p,T_1,qt,ql_1,0)
#         f_1 = thl - thl_1
#         # s_1 = s_dry(pd_1, T_1) * (1.0 - qt) + s_vap(pv_1, T_1) * qt + s_cond(L_1, T_1) * ql_1
#         # f_1 = s - s_1
#         # !!!!
#         T_2 = T_1 + ql_1 * L_1 / ((1.0 - qt) * cpd + qv_star_1 * cpv)
#         delta_T = np.abs(T_2 - T_1)
#
#         count = 0
#         pv_star_2 = CC_Magnus(T_2)  # pv_star_2 = lookup(LT, T_2)
#         qv_star_2 = qv_star_c(p, qt, pv_star_2)
#         ql_2 = qt - qv_star_2
#         while(delta_T >= 1.0e-3 or ql_2 < 0.0):
#             # print('do loop: T2=' + str(T_2)+ ', ql2='+str(ql_2) + ', deltaT='+str(delta_T))
#             pv_star_2 = CC_Magnus(T_2)      # pv_star_2 = lookup(LT, T_2)
#             qv_star_2 = qv_star_c(p, qt, pv_star_2)
#             pv_2 = pv_c(p, qt, qv_star_2)
#             pd_2 = p - pv_2
#             ql_2 = qt - qv_star_2
#             if ql_2 < 0:
#                 print('ql_2 negative in satadj_thl (count='+str(count)+')!', ql_2, 'delta T:', delta_T)
#
#             # !!!!!
#             L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
#             thl_2 = theta_li(p, T_2, qt, ql_2, 0)
#             f_2 = thl - thl_2
#             # s_2 = s_dry(pd_2,T_2)*(1.0-qt) + s_vap(pv_2,T_2)*qt + s_cond(L_2,T_2)*ql_2
#             # f_2 = s - s_2
#             # !!!!!
#             T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
#             T_1 = T_2
#             T_2 = T_n
#             f_1 = f_2
#             delta_T = np.abs(T_2 - T_1)
#             count += 1
#         T = T_2
#         qv = qv_star_2
#         ql = ql_2
#         print('count = ', count)
#
    return T, ql, alpha