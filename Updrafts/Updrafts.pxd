
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron

# cdef class CloudClosure:
#     cpdef initialize(self)
import numpy as np
cimport numpy as np

cdef class Updrafts:
    cdef:
        str path_ref
        str path_out
        double [:] p_ref
        double [:] z_ref
        LatentHeat LH
        ClausiusClapeyron CC
        int [:] krange
        int [:] zrange
        int [:] range
        int [:] times_tracers

    cpdef initialize(self, krange, path, case_name)
    # cpdef update(self, path, path_tr)
    cpdef update(self, files, ncomp_range, dz_range, krange, nml, path, path_tr)
    # ----------------------------------------------------------------------
    #               PDF Model
    # ----------------------------------------------------------------------
    cpdef predict_PDF(self, s_, qt_, T_, ql_, path, int ncomp, dz_range_, krange_, time, nml)

    # cpdef sort_PDF(self, clf):
    cpdef sort_PDF(self, means_, covariance_, weights_, labels_)
    cpdef sort_PDF_allk(self, means_, covariance_, weights_, labels_)

    cpdef create_updrafts_file(self, path, file_name, time, ncomp, nvar, nz_, nml)
    # cpdef write_updrafts_field(self, path, file_name, double[:,:,:] data)
    cpdef write_updrafts_field(self, path, file_name, data, nml)
    cpdef create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_)


    # ----------------------------------------------------------------------
    #                Tracer Model
    # ----------------------------------------------------------------------
    cpdef read_in_updrafts_colleen(self, str type, t, path_)

