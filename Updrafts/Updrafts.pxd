
cdef class CloudClosure:
    cpdef initialize(self)


cdef class Updrafts:
    cdef:
        str path_ref
        str path_out
        double [:] p_ref
        double [:] z_ref
        int [:] zrange

    cpdef initialize(self, krange, path, case_name)
    cpdef predict_pdf(self, files, path, ncomp_range, dz_range_, krange_, nml)

    # cpdef sort_PDF(self, clf):
    cpdef sort_PDF(self, means_, covariance_, weights_, labels_)
    cpdef sort_PDF_allk(self, means_, covariance_, weights_, labels_)

    cpdef create_updrafts_file(self, path, file_name, time, ncomp, nvar, nz_, nml)
    # cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data)
    cpdef write_updrafts_field(self, path, file_name, data, nml)
    cpdef create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_)