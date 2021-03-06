import numpy as np
cimport numpy as np
from cpython cimport array
import array
import cython
from cython cimport floating

cdef class CloudClosure:
    cdef:
        double [:] p_ref
        double [:] zrange
        # array.array p_ref
        # np.ndarray[floating, ndim=3] p_ref
    # cdef double [:] p_ref_a
    # # cdef np.ndarray [:] p_ref_d       # does NOT work; need to give type of elements
    # cdef double [:] p_ref_d
    # cdef double [:] p_ref_e


    cpdef initialize(self, path, path_ref, case_name)
    cpdef verification_CC(self, files, path, path_ref)
    cpdef predict_pdf(self, files, path_in, path_out, path_ref, int n_sample, ncomp_, krange_, nml)
    # cpdef predict_pdf_multilevel(self, files, path_in, path_out, path_ref, ncomp_, krange_, nml)
    # def create_statistics_file(self, path, file_name, time, ncomp, nvar, nz_)


# cpdef do_everything(path)