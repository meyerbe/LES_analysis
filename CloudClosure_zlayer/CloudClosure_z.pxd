

cdef class CloudClosure:
    cdef str path_ref
    cdef double [:] p_ref
    cdef double [:] z_ref
    cdef double [:] zrange

    cpdef initialize(self, path, case_name)
    cpdef predict_pdf(self, files, path, ncomp_range, dz_range_, krange_, nml)






cpdef Gaussian_bivariate(ncomp_, data, var_name1, var_name2, time, z, path)