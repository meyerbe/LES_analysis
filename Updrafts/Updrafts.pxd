
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

    cpdef create_updrafts_file(self, path, file_name, time, ncomp, nvar, nz_, nml)
    # cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data)
    cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data, nml)
    cpdef create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_)