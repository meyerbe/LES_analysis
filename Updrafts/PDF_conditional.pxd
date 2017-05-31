from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron

cdef class PDF_conditional:
    cdef:
        str path_ref
        str path_out
        double [:] p_ref
        double [:] z_ref
        double [:] zrange
        LatentHeat LH
        ClausiusClapeyron CC
        dict nml

    cpdef initialize(self, krange, path, case_name)
    cpdef predict_pdf(self, files, path, int n_sample, ncomp_range, Lx_, Ly_, dk_, int [:] krange_, nml)
