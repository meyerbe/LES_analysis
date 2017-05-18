import numpy as np
cimport numpy as np
import sys
import cython
import time

from scipy.integrate import odeint

sys.path.append("..")
# sys.path.append("../Thermo/")
include '../parameters.pxi'

# __________________________________________________________________
''' Microphysics Schemes:
 No_Microphysics_Dry:
        self.thermodynamics_type = 'dry'
        LatentHeat.Lambda_fp = lambda_constant
        LatentHeat.L_fp = latent_heat_constant
 No_Microphysics_SA:
        self.thermodynamics_type = 'SA'
        LatentHeat.Lambda_fp = lambda_constant
        LatentHeat.L_fp = latent_heat_variable
Microphysics_SB_Liquid:
        self.thermodynamics_type = 'SA'
        LatentHeat.Lambda_fp = lambda_constant
        LatentHeat.L_fp = latent_heat_variable

- Bomex: namelist['microphysics']['scheme'] = 'None_SA' --> 'SA', LH=variable
- ZGILS6: namelist['microphysics']['scheme'] = 'SB_Liquid' --> 'SA', LH=variable
'''




# __________________________________________________________________
cdef class Lookup:
    def __init__(self):
        return

    cpdef initialize(self, double[:] x, double[:] y):
        cdef long n = np.shape(y)[0]
        print('initialize cpdef Lookup Table')
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
cdef class ClausiusClapeyron:
    def __init__(self):
        print('Initialize CC')
        return

    def initialize(self, namelist, LatentHeat LH):
        self.LT = Lookup()

        #Now integrate the ClausiusClapeyron equation
        cdef:
            double Tmin
            double Tmax
            long n_lookup
            double [:] pv

        try:
            Tmin = namelist['ClausiusClapeyron']['temperature_min']
        except:
            print('Clasius-Clayperon lookup table temperature_min not '
                           'given in name list taking default of 180 K')
            Tmin = 100.15

        try:
            Tmax = namelist['ClausiusClapeyron']['temperature_max']
        except:
            print('Clasius-Clayperon lookup table temperature_max not '
                           'given in name list taking default of 340 K')
            Tmax = 380.0

        try:
            n_lookup = namelist['ClausiusClapeyron']['n_lookup']
        except:
            print('Clasius-Clayperon lookup table n_lookup not '
                           'given in name list taking default of 512')
            n_lookup = 512

        #Generate array of equally space temperatures
        T = np.linspace(Tmin, Tmax, n_lookup)
        #Find the maximum index of T where T < T_tilde
        tp_close_index = np.max(np.where(T<=Tt))

        #Check to make sure that T_tilde is not in T
        if T[tp_close_index] == Tt:
            print('Array of temperatures for ClasiusClapyeron lookup table contains Tt  \n')
            print('Pick different values for ClasiusClapyeron Tmin and Tmax in lookup table \n')
            print('Killing Simulation now!')
            sys.exit()

        #Now prepare integration
        T_above_Tt= np.append([Tt],T[tp_close_index+1:])
        T_below_Tt= np.append(T[:tp_close_index+1],[Tt])[::-1]

        #Now set up the RHS
        def rhs(z,T_):
            lam = LH.Lambda(T_)
            L = LH.L(T_,lam)
            return L/(Rv * T_ * T_)

        #set the initial condition
        pv0 = np.log(pv_star_t)

        #Integrate
        pv_above_Tt = np.exp(odeint(rhs,pv0,T_above_Tt,hmax=0.1)[1:])
        pv_below_Tt = np.exp(odeint(rhs,pv0,T_below_Tt,hmax=0.1)[1:])[::-1]
        pv = np.append(pv_below_Tt,pv_above_Tt )
        self.LT.initialize(T,pv)

        # print(self.LT.lookup(298.0))
        # print(self.LT.lookup(296.0))

        return

    cpdef finalize(self):
        self.LT.finalize()
        return
# __________________________________________________________________
# __________________________________________________________________
cdef class LatentHeat:
    def __init__(self,namelist):
        if(namelist['microphysics']['scheme'] == 'None_Dry'):
            print(namelist['microphysics']['scheme'])
            self.Lambda_fp = lambda_constant
            self.L_fp = latent_heat_constant
        elif(namelist['microphysics']['scheme'] == 'None_SA') or (namelist['microphysics']['scheme'] == 'SB_Liquid'):
            print(namelist['microphysics']['scheme'])
            self.Lambda_fp = lambda_constant
            self.L_fp = latent_heat_variable
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
    void eos_c(LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double p0, double s, double qt, double *T, double *qv, double *ql, double *qi) nogil

cdef extern from "thermodynamic_functions.h":
# cdef extern from "../Thermo/thermodynamic_functions.h":
    # Dry air partial pressure
    inline double pd_c(double p0, double qt, double qv) nogil
    # Water vapor partial pressure
    inline double pv_c(double p0, double qt, double qv) nogil


cpdef sat_adj_fromentropy(double p0, double s, double qt, ClausiusClapeyron CC, LatentHeat LH):
    cdef:
        double T, qv, qc, ql, qi, lam
        int alpha = 0
    # eos_c(&self.CC.LT.LookupStructC, self.Lambda_fp, self.L_fp, p0, s, qt, &T, &qv, &ql, &qi)
    eos_c(&CC.LT.LookupStructC, LH.Lambda_fp, LH.L_fp, p0, s, qt, &T, &qv, &ql, &qi)
    if ql > 0.0:
        # print('saturated')
        alpha = 1
    else:
        # print('dry')
        pass
    return T, ql, alpha



cpdef sat_adj_fromthetali(double p, double thl, double qt, ClausiusClapeyron CC, LatentHeat LH):
    # print ('')
    time_a = time.clock()
    # print('saturation adjustment from thetali')
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
    # cdef extern from "thermodynamics_sa.h":
    #     inline double temperature_no_ql(double pd, double pv, double s, double qt)
    cdef extern from "thermodynamic_functions.h":
        inline double pv_c(double p0, const double qt, const double qv)
        inline double qv_star_c(const double p0, const double qt, const double pv)
        inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)
    # from thermodynamic_functions import pv_c, temperature_no_ql_from_thetali, CC_Magnus, qv_star_c, latent_heat
    from thermodynamic_functions import temperature_no_ql_from_thetali

    cdef:
        double qv = qt
        double ql = 0.0
        double T
        double pv_1, pd_1, T_1, pv_star_1, qv_star_1
        double pv_2, pd_2, T_2, T_n, delta_T, pv_star_2, qv_star_2
        double thl_1, thl_2, f_1, f_2
        double ql_1, ql_2

    lookuptable = &CC.LT.LookupStructC

    ''' (1) starting point (first guess): pv, pd, T (assumption: no liquid, i.e. qv=qt) '''
    pv_1 = pv_c(p,qt,qt)
    pd_1 = p - pv_1
    # !!!!!!!
    T_1 = temperature_no_ql_from_thetali(pd_1, pv_1, thl, qt)
    # !!!!!!!
    pv_star_1 = lookup(lookuptable, T_1)
    # pv_star_1 = CC_Magnus(T_1)                    # eos_c: pv_star_1 = lookup(LT, T_1)
    qv_star_1 = qv_star_c(p, qt, pv_star_1)         # eos_c: same

    ''' (2) Check if saturated or not '''
    if (qt <= qv_star_1):
        T = T_1
        ql = 0.0
        # print("not saturated: (thl,qt)=", round(thl,2), round(qt,4), round(T,1))
        alpha = 0
    else:
        ''' (3) if saturated: calculate second starting point T_2 '''
        # print("air saturated: (thl,qt)=", round(thl,2), round(qt,4))
        alpha = 1
        ql_1 = qt - qv_star_1               # eos_c: same; thermo.py = ql_1
        # L_1 = latent_heat(T_1)            # eos_c: L_1 = L_fp(T_1, lam_1)
        L_1 = LH.L(T_1, 1.0)
        # L_1 = LH.L_fp(T_1, 1.0)
        # thl_1 = theta_li(p,T_1,qt,ql_1,0)
        thl_1 = thetali_c(p,T_1,qt,ql_1,0, L_1)
        f_1 = thl - thl_1
        T_2 = T_1 + ql_1 * L_1 / ((1.0 - qt) * cpd + qv_star_1 * cpv)
        delta_T = np.abs(T_2 - T_1)

        count = 0
        pv_star_2 = lookup(lookuptable, T_2)
        # pv_star_2 = CC_Magnus(T_2)  # pv_star_2 = lookup(LT, T_2)
        qv_star_2 = qv_star_c(p, qt, pv_star_2)
        ql_2 = qt - qv_star_2
        while(delta_T >= 1.0e-3 or ql_2 < 0.0):
            # print('do loop: T2=' + str(T_2)+ ', ql2='+str(ql_2) + ', deltaT='+str(delta_T))
            pv_star_2 = lookup(lookuptable, T_2)
            # pv_star_2 = CC_Magnus(T_2)      # pv_star_2 = lookup(LT, T_2)
            qv_star_2 = qv_star_c(p, qt, pv_star_2)
            pv_2 = pv_c(p, qt, qv_star_2)
            pd_2 = p - pv_2
            ql_2 = qt - qv_star_2
            # L_2 = latent_heat(T_2)      # eos_c: L_2 = L_fp(T_2,lam_2)
            L_2 = LH.L(T_2, 1.0)
            # L_2 = LH.L_fp(T_2, 1.0)
            # thl_2 = theta_li(p, T_2, qt, ql_2, 0)
            thl_2 = thetali_c(p, T_2, qt, ql_2, 0, L_2)
            f_2 = thl - thl_2
            T_n = T_2 - f_2*(T_2 - T_1)/(f_2 - f_1)
            T_1 = T_2
            T_2 = T_n
            f_1 = f_2
            delta_T = np.abs(T_2 - T_1)
            count += 1
        T = T_2
        qv = qv_star_2
        ql = ql_2
        # print('count = ', count)
    time_b = time.clock()
    return T, ql, alpha#, time_b-time_a
