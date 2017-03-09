# sys.path.append("..")
# from parameters import *
# sys.path.append("../Thermo/")

# cdef struct eos_struct:
#     double T
#     double ql
#     double alpha

#----------------------------------------------------------------------
# cpdef eos_struct sat_adj_fromentropy(double p, double s, double qt)
#
# cpdef eos_struct sat_adj(double p, double s, double qt)

cdef extern from "lookup.h":
    struct LookupStruct:
        int n_table
        double y_max
        double y_min
        double x_max
        double x_min
        double dx
        double dxi
        double* table_values

    void init_table(LookupStruct *LT, long n, double *x, double *y) nogil
    void free_table(LookupStruct *LT) nogil
    inline double lookup(LookupStruct * LT, double x) nogil


cdef class Lookup:
    cdef:
        LookupStruct LookupStructC
    cpdef initialize(self, double[:] x, double[:] y)
    cpdef finalize(self)
    cpdef table_bounds(self)
    cpdef lookup(self, double x)
    cdef:
        inline double fast_lookup(self, double x) nogil


cdef class ClausiusClapeyron:
    cdef Lookup LT
    #def initialize(self,namelist,LatentHeat LH, ParallelMPI.ParallelMPI Par)
    cpdef finalize(self)


cdef class LatentHeat:
    cdef:
        #In the functions pointed to by the function pointer L* must not require gil
        double (*L_fp)(double T, double Lambda) nogil
        double (*Lambda_fp)(double T) nogil
    cpdef L(self, double T, double Lambda)
    cpdef Lambda(self, double T)


cdef inline double lambda_constant(double T) nogil:
    return 1.0

cdef inline double latent_heat_constant(double T, double Lambda) nogil:
    return 2.501e6

cdef inline double latent_heat_variable(double T, double Lambda) nogil:
    cdef:
        double TC = T - 273.15
    return (2500.8 - 2.36 * TC + 0.0016 * TC *
            TC - 0.00006 * TC * TC * TC) * 1000.0


cpdef sat_adj_fromentropy(double p0, double s, double qt, ClausiusClapeyron CC, LatentHeat LH)
cpdef sat_adj_fromthetali_c(double p, double thl, double qt, ClausiusClapeyron CC, LatentHeat LH)



