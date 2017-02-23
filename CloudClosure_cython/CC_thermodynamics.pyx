import numpy as np
cimport numpy as np
import sys
import cython

from scipy.integrate import odeint

sys.path.append("..")
include '../parameters.pxi'

def aux():
    print('hoihoi')
    return

# cdef aux_c():
#     print('hoihoi c')
#     return

class LookupTable:
    def __init__(self):
        print('initialize Lookup Table')
        self.lookup_table = {}
        self.x_min = 9e16
        self.x_max = -9e16
        self.y_min = 9e16
        self.y_max = -9e16
        # self.dx =
        # self.dxi
        return

    def make_lookup_table(self, x,y):
        if x.size == y.size:
            n = x.size
            print('x.size = y.size: ', x.size)
        else:
            n = np.min(x.size, y.size)
            print('x.size != y.size')
        self.x_min = min(np.amin(x), self.x_min)
        self.x_max = max(np.amax(x), self.x_max)
        self.y_min = min(np.amin(y), self.y_min)
        self.y_max = max(np.amax(y), self.y_max)

        for i in range(n):
            self.lookup_table[x[i]] = y[i]

        self.dx = x[1] - x[0]
        self.dxi = 1.0/self.dx
        return

    def lookup(self, x, x_in):
        # indx = np.int((x_in-self.x_min)*self.dxi)
        indx = ((x_in-self.x_min)*self.dxi)
        print('indx: ', indx, type(indx), self.y_max, len(self.lookup_table))
        indx = np.int(indx)
        x0 = self.x_min + self.dx*indx
        print('x0', x0, self.dx, indx)
        y1 = self.lookup_table[x[indx]]
        y2 = self.lookup_table[x[indx+1]]
        x_out = y1 + (x_in-x0)*(y2-y1)*self.dxi
        return x_out


cdef class LookupTable_c:
    def __init__(self):
        print('initialize Lookup Table')
        self.lookup_table = {}
        self.x_min = 9e16
        self.x_max = -9e16
        self.y_min = 9e16
        self.y_max = -9e16
        # self.dx =
        # self.dxi
        return

    cdef make_lookup_table(self, double[:] x, double[:] y):
        cdef int n
        if x.size == y.size:
            n = x.size
            print('x.size = y.size: ', x.size)
        else:
            n = np.min(x.size, y.size)
            print('x.size != y.size')
        self.x_min = min(np.amin(x), self.x_min)
        self.x_max = max(np.amax(x), self.x_max)
        self.y_min = min(np.amin(y), self.y_min)
        self.y_max = max(np.amax(y), self.y_max)

        for i in range(n):
            self.lookup_table[x[i]] = y[i]

        self.dx = x[1] - x[0]
        self.dxi = 1.0/self.dx
        return

    cdef lookup(self, double[:] x, double x_in):
        cdef double x_out, y1, y2, x0
        # indx = np.int((x_in-self.x_min)*self.dxi)
        indx = ((x_in-self.x_min)*self.dxi)
        print('indx: ', indx, type(indx), self.y_max, len(self.lookup_table))
        indx = np.int(indx)
        x0 = self.x_min + self.dx*indx
        print('x0', x0, self.dx, indx)
        y1 = self.lookup_table[x[indx]]
        y2 = self.lookup_table[x[indx+1]]
        x_out = y1 + (x_in-x0)*(y2-y1)*self.dxi
        return x_out


cdef class LookupTable_h:
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




class ClausiusClapeyron:
    def __init__(self):
        print('CC')
        return

    def initialize(self):
        cdef:
            double Tmin, Tmax
            long n_lookup
            double [:] pv

        Tmin = 100.15
        Tmax = 380.0
        n_lookup = 512
        # try:
        #     Tmin = namelist['ClausiusClapeyron']['temperature_min']
        # except:
        #     Par.root_print('Clasius-Clayperon lookup table temperature_min not '
        #                    'given in name list taking default of 180 K')
        #     Tmin = 100.15
        # try:
        #     Tmax = namelist['ClausiusClapeyron']['temperature_max']
        # except:
        #     Par.root_print('Clasius-Clayperon lookup table temperature_max not '
        #                    'given in name list taking default of 340 K')
        #     Tmax = 380.0
        #
        # try:
        #     n_lookup = namelist['ClausiusClapeyron']['n_lookup']
        # except:
        #     Par.root_print('Clasius-Clayperon lookup table n_lookup not '
        #                    'given in name list taking default of 128')
        #     n_lookup = 512
#

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

        # Now set up the RHS
        def rhs(z, T_):
            sys.path.append('../Thermo/')
            from thermodynamic_functions import latent_heat
            # lam = LH.Lambda(T_)
            # L = LH.L(T_, lam)
            lam  = 1.0
            L = latent_heat(T_)
            return L / (Rv * T_ * T_)

        # set the initial condition
        pv0 = np.log(pv_star_t)

        # Integrate
        pv_above_Tt = np.exp(odeint(rhs, pv0, T_above_Tt, hmax=0.1)[1:])
        pv_below_Tt = np.exp(odeint(rhs, pv0, T_below_Tt, hmax=0.1)[1:])[::-1]
        pv = np.append(pv_below_Tt, pv_above_Tt)
        # self.LT.initialize(T, pv)

        print('T:')
        # print(np.array(T))
        print('pv: ', pv.shape)
        # print(np.array(pv))

        LT = LookupTable()
        print(LT.lookup_table)
        LT.make_lookup_table(T,pv)
        print('given value of T: ', LT.lookup_table[T[0]])        # only possible to lookup for values of T that are already calculated --> need interpolation
        print('interpolated value of T: ', LT.lookup(T, 298.0))
        return




cdef class ClausiusClapeyron_c:
    def __init__(self):
        print('CC c')
        return

    cpdef rhs(self, z, T_):
        sys.path.append('../Thermo/')
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

        LT = LookupTable_h()
        LT.initialize(T,pv)
        print(LT.lookup(298.0))
        print(LT.lookup(296.0))
        print('pv: ', pv)

        LT_c = LookupTable_c()
        return