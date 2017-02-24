from CC_thermodynamics_c cimport Lookup

# cdef extern from "thermodynamics_sa.h":
# cdef extern from "../Thermo/thermodynamics_sa.h":
    # inline double alpha_c(double p0, double T, double qt, double qv) nogil
    # void eos_c(Lookup.LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double p0, double s, double qt, double *T, double *qv, double *ql, double *qi) nogil
    # void eos_update(Grid.DimStruct *dims, Lookup.LookupStruct *LT, double(*lam_fp)(double), double(*L_fp)(double, double), double *p0, double *s, double *qt, double *T,
    #                 double * qv, double * ql, double * qi, double * alpha)
    # void buoyancy_update_sa(Grid.DimStruct *dims, double *alpha0, double *alpha, double *buoyancy, double *wt)
    # void bvf_sa(Grid.DimStruct * dims, Lookup.LookupStruct * LT, double(*lam_fp)(double), double(*L_fp)(double, double), double *p0, double *T, double *qt, double *qv, double *theta_rho, double *bvf)
    # void thetali_update(Grid.DimStruct *dims, double (*lam_fp)(double), double (*L_fp)(double, double), double *p0, double *T, double *qt, double *ql, double *qi, double *thetali)
    # void clip_qt(Grid.DimStruct *dims, double  *qt, double clip_value)

cdef extern from "thermodynamic_functions.h":
# cdef extern from "../Thermo/thermodynamic_functions.h":
    # Dry air partial pressure
    inline double pd_c(double p0, double qt, double qv) nogil
    # Water vapor partial pressure
    inline double pv_c(double p0, double qt, double qv) nogil