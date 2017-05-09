from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron

cdef class CloudClosure:
    cdef:
        str path_ref
        str path_out
        double [:] p_ref
        double [:] z_ref
        double [:] zrange

    cpdef initialize(self, krange, path, case_name)
    cpdef predict_pdf(self, files, path, int n_sample, ncomp_range, dk_range, int [:] krange_, nml)
    cpdef sample_pdf(self, data, clf, double ql_mean_ref, double cf_ref, double pref,
                     ClausiusClapeyron CC, LatentHeat LH, ncomp_range, n_sample, nml)







cpdef Gaussian_bivariate(ncomp_, data, var_name1, var_name2, time, z, path)