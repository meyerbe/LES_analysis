sys.path.append("..")
from parameters import *
sys.path.append("../Thermo/")

# cdef struct eos_struct:
#     double T
#     double ql
#     double alpha

#----------------------------------------------------------------------
# cpdef eos_struct sat_adj_fromentropy(double p, double s, double qt)
#
# cpdef eos_struct sat_adj(double p, double s, double qt)

# cpdef ClausiusClapeyron

# cdef aux_c()

cdef class ClausiusClapeyron_c:
    cdef:
        double Tmin
        double Tmax
        long n_lookup
        double [:] pv

    cpdef rhs(self, z, T_)
    cpdef initialize(self)


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

cdef class LookupTable_c:
    cdef:
        dict lookup_table
        double x_min, x_max, y_min, y_max
        double dx, dxi
    cdef make_lookup_table(self, double[:] x, double[:] y)
    cdef lookup(self, double[:] x, double x_in)

cdef class LookupTable_h:
    cdef:
        LookupStruct LookupStructC
    cpdef initialize(self, double[:] x, double[:] y)
    cpdef finalize(self)
    cpdef table_bounds(self)
    cpdef lookup(self, double x)
    cdef:
        inline double fast_lookup(self, double x) nogil
