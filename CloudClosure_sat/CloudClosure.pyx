import os
include 'parameters.pxi'
import numpy as np
cimport numpy as np
import netCDF4 as nc
import json as simplejson
import pylab as plt
import time

from cpython cimport array
import array
import cython
from cython cimport floating
from sklearn import mixture
from sklearn.preprocessing import StandardScaler


import CC_thermodynamics_c
from CC_thermodynamics_c cimport LatentHeat, ClausiusClapeyron
from CC_thermodynamics_c import sat_adj_fromentropy, sat_adj_fromthetali
# from plotting_functions import plot_PDF_samples_qt, plot_PDF_samples, plot_sat_adj, plot_PDF
from plotting_functions import plot_PDF, plot_samples, plot_hist, plot_error_vs_ncomp


cdef class CloudClosure:
    def __init__(self):
        print('Cloud Closure')
        self.p_ref = None
        return
    # cpdef do_everything(path, path_ref):
    #     path_fields = os.path.join(path, 'fields', )
    #
    #     # (1) Reference State
    #     try:
    #         p_ref = read_in_netcdf('p0_half', 'reference', path_ref)
    #     except:
    #         print('no p0_half profile')
    #         p_ref = read_in_netcdf('p0', 'reference', path_ref)
    #
    #     # (2) Fields
    #     cdef:
    #         double [:,:,:] s_
    #         double [:,:,:] qt_
    #         double [:,:,:] T_
    #         double [:,:,:] ql_
    #     files = os.listdir(os.path.join(path, 'fields'))
    #     N = len(files)
    #     print('Found the following directories', files, N)
    #     nc_file_name = str(files[0])
    #     path_fields = os.path.join(path_fields, nc_file_name)
    #     s_ = read_in_netcdf('s', 'fields', path_fields)
    #     qt_ = read_in_netcdf('qt', 'fields', path_fields)
    #     T_ = read_in_netcdf('temperature', 'fields', path_fields)
    #     ql_ = read_in_netcdf('ql', 'fields', path_fields)
    #     print('')
    #
    #     cdef:
    #         int i, j, k
    #         double s, qt, T, ql
    #         double pv, pd
    #         double T_unsat
    #     i = 10
    #     j = 10
    #     for k in [10]:
    #         s = s_[i,j,k]
    #         qt = qt_[i,j,k]
    #         ql = ql_[i,j,k]
    #         T = T_[i,j,k]
    #         pref = p_ref[k]
    #         qv = np.double(qt, copy = True)
    #
    #         pv = pref * eps_vi * qv / (1.0 - qt + eps_vi * qv)
    #         pd = pref - pv
    #
    #         T_unsat = T_tilde * np.exp((s -
    #                                     (1.0 - qt)*(sd_tilde - Rd*np.log(pd/p_tilde))
    #                                     - qt*(sv_tilde - Rv*np.log(pv/p_tilde)) )
    #                                     / ( (1.0-qt)*cpd + qt*cpv))
    #
    #         print('difference: ', T - T_unsat)
    #         print('ql=', ql)
    #         print('')
    #         # print('')
    #
    #     return

    cpdef initialize(self, path, path_ref, case_name):
        print('')
        print('--- Initialize ---')
        print('nml: ', os.path.join(path, case_name+'.in'))
        nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
        cdef int nx, ny, nz, gw
        nx = nml['grid']['nx']
        ny = nml['grid']['ny']
        nz = nml['grid']['nz']
        gw = nml['grid']['gw']

        # (1) Reference State

        # # p_ref_a: [75, 0, 0, 0, 0, 0, 0, 0]
        # self.p_ref_a = np.zeros(shape = (nz),dtype=np.double,order='c')                 # <type 'CloudClosure._memoryviewslice'>
        # # p_ref_b: [75, 0, 0, 0, 0, 0, 0, 0]
        # cdef double [:] p_ref_b = np.zeros(shape = (nz),dtype=np.double,order='c')      # <type 'CloudClosure._memoryviewslice'>
        # # p_ref_c: (75,)
        # cdef np.ndarray p_ref_c = np.zeros([nz], dtype=np.double,order='c')             # <type 'numpy.ndarray'>
        # # p_ref_d: (75,)
        # self.p_ref_d = np.zeros(shape = np.shape(p_ref_c),dtype=np.double,order='c')    # <type 'CloudClosure._memoryviewslice'>
        # # p_ref_e: (75,)
        # self.p_ref_e = np.zeros([nz],dtype=np.double,order='c')                         # <type 'CloudClosure._memoryviewslice'>
        # # cdef np.ndarray self.p_ref_d = np.zeros([nz], dtype=np.double,order='c')  # does NOT work
        # # cdef array.array p_ref
        # # cdef np.ndarray[floating, ndim=3] p_ref = np.empty((nz))
        #
        # print('ref', type(self.p_ref_a), type(p_ref_b), type(p_ref_c), type(self.p_ref_d), type(self.p_ref_e))
        #
        # print(self.p_ref_a.shape, np.shape(self.p_ref_a))
        # print(p_ref_b.shape, np.shape(p_ref_b))
        # print(np.shape(p_ref_c))
        # print(np.shape(self.p_ref_d))
        # print(np.shape(self.p_ref_e))
        #
        # self.p_ref_a = read_in_netcdf('p0', 'reference', path_ref)
        # p_ref_b = read_in_netcdf('p0', 'reference', path_ref)
        # p_ref_c = read_in_netcdf('p0', 'reference', path_ref)
        # self.p_ref_d = read_in_netcdf('p0', 'reference', path_ref)
        # self.p_ref_e = read_in_netcdf('p0', 'reference', path_ref)
        #
        # print(self.p_ref_a.shape)
        # print(p_ref_b.shape)
        # print(np.shape(p_ref_c))
        # print(np.shape(self.p_ref_d))
        # print(np.shape(self.p_ref_e))
        #
        # # try:
        # #     print('diff ab', self.p_ref_a-p_ref_b)
        # # except:
        # #     pass
        # # try:
        # #     print('diff ac', self.p_ref_a-p_ref_c)
        # # except:
        # #     pass
        #
        # print('diff bc', p_ref_b-p_ref_c)
        # # print('diff bc', self.p_ref_d-self.p_ref_e)
        # print('diff d e', self.p_ref_d[10] - self.p_ref_e[10])
        #
        # for k in range(p_ref_c.shape[0]):
        #     self.p_ref_a[k] = p_ref_c[k]
        #
        # print(self.p_ref_a.shape)
        # print(p_ref_b.shape)
        # print(np.shape(p_ref_c))
        # print(np.shape(self.p_ref_d))
        # print(np.shape(self.p_ref_e))

        self.p_ref = np.zeros([nz],dtype=np.double,order='c')                         # <type 'CloudClosure._memoryviewslice'>
        try:
            self.p_ref = read_in_netcdf('p0_half', 'reference', path_ref)
        except:
            print('no p0_half profile')
            self.p_ref = read_in_netcdf('p0', 'reference', path_ref)[:]#
        # print('p_ref', np.shape(self.p_ref), self.p_ref.shape, type(self.p_ref[0]), self.p_ref[0])

        return


    # ----------------------------------------------------------------------------------------------------
    cpdef verification_CC(self, files, path, path_ref):
        print('Verification')
        cdef extern from "thermodynamics_sa.h":
            inline double temperature_no_ql(double pd, double pv, double s, double qt)
        cdef extern from "thermodynamic_functions.h":
            inline double pd_c(const double p0,const double qt, const double qv)
            inline double pv_c(const double p0, const double qt, const double qv)
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)

        path_fields = os.path.join(path, 'fields', )
        cdef:
            int i, j, k

        # (0) Namelist File
        nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
        # cdef int nx, ny, nz, gw
        cdef:
            nx = nml['grid']['nx']
            ny = nml['grid']['ny']
            nz = nml['grid']['nz']
            gw = nml['grid']['gw']

        # (1) Reference State
        cdef double [:] p_ref = np.zeros([nz],dtype=np.double,order='c')                         # <type 'CloudClosure._memoryviewslice'>
        for k in range(nz):
            p_ref[k] = self.p_ref[k]
        # p_ref = self.p_ref
        print('p_ref: ', np.shape(self.p_ref), self.p_ref.shape, p_ref.shape, nz, gw)


        # (2) Fields
        cdef:
            double [:,:,:] s_
            double [:,:,:] qt_
            double [:,:,:] T_
            double [:,:,:] ql_
        # files = os.listdir(os.path.join(path, 'fields'))
        # N = len(files)
        # print('Found the following directories', files, N)
        nc_file_name = str(files[0])
        path_fields = os.path.join(path_fields, nc_file_name)
        var_list = ['s', 'qt', 'temperature', 'ql']
        s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)
        # s_ = read_in_netcdf('s', 'fields', path_fields)
        # qt_ = read_in_netcdf('qt', 'fields', path_fields)
        # T_ = read_in_netcdf('temperature', 'fields', path_fields)
        # ql_ = read_in_netcdf('ql', 'fields', path_fields)
        print('')


        # (3) temperature no ql
        print('--- temperature no ql ---')
        print('')
        cdef:
            double s, qt, T, ql, qv
            double pv, pd, pref
            double T_unsat
            double min_T_unsat = 9999.9
            double max_T_unsat = -9999.9
        for i in range(nx):
            for j in range(ny):
                # for k in range(nz):
                for k in range(10,20):
                    s = s_[i,j,k]
                    qt = qt_[i,j,k]
                    ql = ql_[i,j,k]
                    T = T_[i,j,k]
                    pref = p_ref[k]
                    qv = qt

                    pv = pv_c(pref, qt, qv)
                    pd = pd_c(pref, qt, qv)
                    # pd = pref - pv

                    T_unsat = temperature_no_ql(pd, pv, s, qt)
                    if ql == 0.0:
                        if (T-T_unsat) > max_T_unsat:
                            max_T_unsat = T-T_unsat
                        elif (T-T_unsat) < min_T_unsat:
                            min_T_unsat = T-T_unsat

        print('difference T unsat min: ', min_T_unsat)
        print('difference T unsat max: ', max_T_unsat)
        print('')
        print('')


        ###############################################
        # Saturation Adjustment
        # (3) saturation adjustment
        cdef:
            LatentHeat LH = CC_thermodynamics_c.LatentHeat(nml)
            ClausiusClapeyron CC = CC_thermodynamics_c.ClausiusClapeyron()
        print('done Lookup table')
        CC.initialize(nml, LH)
        alpha_ = np.zeros(shape=(nx,ny,nz))
        T_comp = np.zeros(shape=(nx,ny,nz))
        ql_comp = np.zeros(shape=(nx,ny,nz))
        alpha = np.zeros(shape=(nx,ny,nz))
        theta_l = np.zeros(shape=(nx,ny,nz))
        # double [:,:,:] theta_l = np.zeros([nx,ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
        T_comp_thl = np.zeros(shape=(nx,ny,nz))
        ql_comp_thl = np.zeros(shape=(nx,ny,nz))
        alpha_thl = np.zeros(shape=(nx,ny,nz))
        qi_ = 0
        n_alpha = 0
        # print('alpha shape', alpha_.shape)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if ql_[i,j,k] != 0.0:
                        alpha_[i,j,k] = 1.0
                        n_alpha += 1

        cdef:
            double max_T_sat = -9999.9
            double max_ql = -9999.9
            double min_ql = 9999.9
            double max_T_sat_thl = -9999.9
            double max_T_unsat_thl = -9999.9
            double max_ql_thl = -9999.9
            double min_ql_thl = 9999.9
            # double s, qt
            double al
            double one, two, three


        print('')
        print('--- Saturation Adjustment ---')
        time1 = time.clock()
        for i in range(nx):
            # print('i', i, nx)
            for j in range(ny):
                for k in range(15,16):
                # for k in range(nz):
                    s = s_[i,j,k]
                    qt = qt_[i,j,k]
                    ql = ql_[i,j,k]
                    T = T_[i,j,k]
                    Lv = LH.L(T_[i,j,k],LH.Lambda_fp(T_[i,j,k]))
                    theta_l[i,j,k] = thetali_c(p_ref[k], T_[i,j,k], qt_[i,j,k], ql_[i,j,k], qi_, Lv)

                    T_comp[i,j,k], ql_comp[i,j,k], alpha[i,j,k] = sat_adj_fromentropy(p_ref[k], s, qt, CC, LH)
                    T_comp_thl[i,j,k], ql_comp_thl[i,j,k], alpha_thl[i,j,k] = sat_adj_fromthetali(p_ref[k], theta_l[i,j,k], qt_[i,j,k], CC, LH)

                    if (ql_comp[i,j,k] - ql) > max_ql:
                        max_ql = (ql_comp[i,j,k] - ql)
                    elif (ql_comp[i,j,k] - ql) < min_ql:
                        min_ql = (ql_comp[i,j,k] - ql)
                    if ql > 0.0:# and alpha_[i,j,k] > 0.0:
                        if np.abs(T_comp[i,j,k] - T) > max_T_sat:
                            max_T_sat = np.abs(T_comp[i,j,k] - T)
                        if np.abs(T_comp_thl[i,j,k] - T) > max_T_sat_thl:
                            max_T_sat_thl = np.abs(T_comp_thl[i,j,k] - T)
                    # elif alpha_[i,j,k] == 0.0:
                    else:
                        # print('not saturated')
                        if np.abs(T_comp[i,j,k] - T) > max_T_unsat:
                            max_T_unsat = np.abs(T_comp[i,j,k] - T)
                        if np.abs(T_comp_thl[i,j,k]-T) > max_T_unsat_thl:
                            max_T_unsat_thl = np.abs(T_comp_thl[i,j,k]-T)
                    if (ql_comp_thl[i,j,k] - ql) > max_ql_thl:
                        max_ql_thl = ql_comp_thl[i,j,k] - ql
                    elif (ql_comp_thl[i,j,k] - ql) < min_ql_thl:
                            min_ql_thl = ql_comp_thl[i,j,k] - ql
        time2 = time.clock()

        print('')
        print('From Entropy (CC):')
        print('max T sat: ', max_T_sat)             # max_T_sat = 0.096
        print('max T unsat: ', max_T_unsat)         # max_T_unsat = 0.05
        print('max ql:', max_ql)                    # max_ql = 4.4e-5
        print('min ql:', min_ql)                    # min_ql = -6.7e-5
        print('n_alpha: ', n_alpha)
        print('')
        print('From Thetali:')
        print('max T sat: ', max_T_sat_thl)         # max_T_sat = 0.12
        print('max T unsat: ', max_T_unsat_thl)     # max_T_unsat = 0.013
        print('max ql:', max_ql_thl)                # max_ql = 3.31e-5
        print('min ql:', min_ql_thl)                # min_ql = 0.0
        print('')
        print('execution time: ', time2-time1)
        print('')
        return



    # ----------------------------------------------------------------------------------------------------
    cpdef predict_pdf(self, files, path_in, path_out, path_ref, ncomp_range, krange_, nml):
        print('--- PDF Prediction ---')
        print('')
        # cdef extern from "thermodynamics_sa.h":
        #     inline double temperature_no_ql(double pd, double pv, double s, double qt)
        cdef extern from "thermodynamic_functions.h":
            # inline double pd_c(const double p0,const double qt, const double qv)
            # inline double pv_c(const double p0, const double qt, const double qv)
            inline double thetali_c(const double p0, const double T, const double qt, const double ql, const double qi, const double L)
        time1 = time.clock()


        # ________________________________________________________________________________________
        '''Set up files etc.'''
        cdef:
            int ncomp
            int nvar = 2
            Py_ssize_t [:] krange = krange_
            int nk = len(krange)
            Py_ssize_t k, iz
            str d
            int nx = nml['grid']['nx']
            int ny = nml['grid']['ny']
            int nz = nml['grid']['nz']
            int dz = nml['grid']['dz']
            int dt_stat = nml['stats_io']['frequency']
            int nt
            double [:] p_ref = np.zeros([nz],dtype=np.double,order='c')       # <type 'CloudClosure._memoryviewslice'>

        time_profile = read_in_netcdf('t', 'timeseries', path_ref)
        nt = time_profile.shape[0]
        print('time in stats file: ', time_profile.shape, nt, dt_stat)

        # files = os.listdir(os.path.join(path_in,'fields'))
        N = len(files)
        print('Found the following fields: ', files, N)

        for k in range(nz):
            p_ref[k] = self.p_ref[k]
        ql_profile = read_in_netcdf('ql_mean', 'profiles', path_ref)
        # print('ql profile:', ql_profile[krange[0:nk]])#ql_profile.shape, ql_profile[krange[0]])
        # print(ql_profile[krange[1]])

        '''Initialize Latent Heat and ClausiusClapeyron'''
        print('')
        print('Initializing Clausius Clapeyron')
        cdef:
            LatentHeat LH = CC_thermodynamics_c.LatentHeat(nml)
            ClausiusClapeyron CC = CC_thermodynamics_c.ClausiusClapeyron()
        CC.initialize(nml, LH)
        print('')

        # ________________________________________________________________________________________
        '''(A) Compute PDF f(s,qt) from LES data'''
        #       - read in fields
        #       - compute theta_l
        #       - compute PDFs f(s,qt), g(th_l, qt)
        '''(1) Statistics File'''
        nc_file_name_out = 'CC_alltime.nc'
        # create_statistics_file(os.path.join(fullpath_out, 'CloudClosure'), nc_file_name_out, ncomp, nvar, len(zrange))

        '''(2) Read in Fields'''
        # for PDF construction
        cdef:
            double [:,:,:] s_
            double [:,:,:] T_
            double [:,:,:] qt_
            double [:,:] qt = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:,:,:] ql_
            double [:,:] ql = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double qi_ = 0.0
            double [:,:] theta_l = np.zeros([nx*ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            # double [:,:,:] T_comp = np.zeros([nx,ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            # double [:,:,:] ql_comp = np.zeros([nx,ny,nk],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            int ishift = ny
            int i, j, ij
            double Lv
        # theta_l_norm = np.zeros([nx*ny*N],dtype=np.double,order='c')               # <type 'CloudClosure._memoryviewslice'>
        # qt_norm = np.zeros([nx*ny*N],dtype=np.double,order='c')                    # <type 'CloudClosure._memoryviewslice'>

        # for PDF sampling
        cdef:
            int n_sample = np.int(1e6)
            double [:] T_comp = np.zeros([n_sample],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:] ql_comp = np.zeros([n_sample],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            # int [:] alpha_comp = np.zeros([n_sample],dtype=np.int_,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:] T_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            double [:] ql_comp_thl = np.zeros([n_sample],dtype=np.double,order='c')         # <type 'CloudClosure._memoryviewslice'>
            # int [:] alpha_comp_thl = np.zeros([n_sample],dtype=np.int_,order='c')         # <type 'CloudClosure._memoryviewslice'
            double [:,:] Th_l= np.zeros([n_sample, nvar], dtype=np.double, order='c')
            double ql_mean = 0.0
            double [:,:] error_array = np.zeros(shape=(nk,len(ncomp_range)))
            double [:,:] rel_error_array = np.zeros(shape=(nk,len(ncomp_range)))
        alpha_comp = np.zeros(n_sample)
        alpha_comp_thl = np.zeros(n_sample)

        # for <ql> intercomparison
        cdef:
            double ql_mean_profile = 0.0
            double ql_mean_field
            int n, n1, n2
            double error
            double rel_error = 0.0
            int count_ncomp = 0

        var_list = ['s', 'qt', 'temperature', 'ql']
        data = np.ndarray(shape=((nx * ny), nvar))
        for ncomp in ncomp_range:
            means_ = np.ndarray(shape=(nk, ncomp, nvar))
            covariance_ = np.zeros(shape=(nk, ncomp, nvar, nvar))
            for k in range(nk):
                iz = krange[k]
                print('- z = '+str(iz*dz)+ ', ncomp = '+str(ncomp)+' -')
                data_all = np.ndarray(shape=(0, nvar))
                ql_mean_field = 0.0
                time_a = time.clock()
                for d in files:
                    nc_file_name = d
                    path_fields = os.path.join(path_in, 'fields', nc_file_name)
                    print('path_fields', path_fields)
                    s_, qt_, T_, ql_ = read_in_fields('fields', var_list, path_fields)
                    '''(3) Compute liquid potential temperature from temperature and moisture'''
                    for i in range(nx):
                        for j in range(ny):
                            ij = i*ishift + j
                            Lv = LH.L(T_[i,j,iz],LH.Lambda_fp(T_[i,j,iz]))
                            theta_l[ij,k] = thetali_c(p_ref[iz], T_[i,j,iz], qt_[i,j,iz], ql_[i,j,iz], qi_, Lv)
                            qt[ij,k] = qt_[i,j,iz]
                            ql[ij,k] = ql_[i,j,iz]
                            ql_mean_field += ql_[i,j,iz]
                    del s_, ql_, qt_, T_
                    ql_mean_field /= (nx*ny)
                    data[:, 0] = theta_l[:, k]
                    data[:, 1] = qt[:, k]
                    # data[:, 0] = theta_l_norm[:, k]
                    # data[:, 1] = qt_norm[:, k]
                    data_all = np.append(data_all, data, axis=0)
                time_b = time.clock()
                print('time to read in data_all for z='+str(iz*dz)+': '+str(time_b-time_a))
                time_a = time.clock()

                '''(4) Normalise Data'''
                # X = StandardScaler().fit_transform(X)
                # theta_l_norm = StandardScaler().fit_transform(data_all[:,0])#[:,0]
                # qt_norm = StandardScaler().fit_transform(data[:,1])#[:,0]
                scaler = StandardScaler()
                data_all_norm = scaler.fit_transform(data_all)
                # print('normalise')
                # print(theta_l.shape, theta_l_norm.shape, qt.shape, qt_norm.shape)
                # print(theta_l.shape, qt.shape, data_all.shape, data_all_norm.shape)
                print('')
                '''(5) Compute bivariate PDF'''
                #   (a) for (s,qt)
                #           ...
                #   (b) for (th_l,qt)
                # clf_thl = Gaussian_bivariate(ncomp, data_all, 'T', 'qt', np.int(d[0:-3]), iz * dz, path_in)
                clf_thl_norm = Gaussian_bivariate(ncomp, data_all_norm, 'T', 'qt', np.int(d[0:-3]), iz * dz, path_in)
                # # means_[k, :, :] = clf_thl.means_[:, :]
                # # covariance_[k,:,:,:] = clf_thl.covariances_[:,:,:]
                means_[k, :, :] = clf_thl_norm.means_[:, :]
                covariance_[k,:,:,:] = clf_thl_norm.covariances_[:,:,:]
                time_b = time.clock()
                # print('time for GMM for z='+str(iz*dz)+': '+str(time_b-time_a))

                # '''(6) Compute Kernel-Estimate PDF '''
                # # kde, kde_aux = Kernel_density_estimate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz, path_in)
                #
                # # relative_entropy(data, clf, kde)
                # print('')



                '''(B) Compute mean liquid water <ql> from PDF f(s,qt)'''
                #       1. sample (th_l, qt) from PDF (Monte Carlo ???
                #       2. compute ql for samples
                #       3. consider ensemble average <ql> = domain mean representation???
                time_a = time.clock()

                '''(1) Draw samples'''
                # S, y = clf.sample(n_samples=n_sample)
                # S, y = clf_s.sample(n_samples=nn)
                # S_norm, y = clf_s_norm.sample(n_samples=nn)
                # print('clf samples: ', S.shape, y.shape)
                # Th_l, y = clf_thl.sample(n_samples=n_sample)
                # print('clf thl samples: ', Th_l.shape, y.shape)
                Th_l_norm, y_norm = clf_thl_norm.sample(n_samples=n_sample)
                # print('clf thl samples: ', Th_l_norm.shape, y_norm.shape)

                time_b = time.clock()
                # print('time clf.sample: ', time_b-time_a)

                '''(2) Rescale theta_l and qt'''
                # print(scaler.get_params())       # {'copy': True, 'with_mean': True, 'with_std': True}
                Th_l = scaler.inverse_transform(Th_l_norm)      # Inverse Normalisation

                '''(3) Compute ql (saturation adjustment)'''
                print('!!!!', n_sample, np.shape(T_comp_thl), np.shape(ql_comp_thl), np.shape(alpha_comp_thl), np.shape(Th_l), np.shape(Th_l_norm))
                for i in range(n_sample-2):
                    # T_comp[i], ql_comp[i], alpha_comp[i] = sat_adj_fromentropy(p_ref[iz-1], S[i,0],S[i,1])
                    T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)
                    ql_mean = ql_mean + ql_comp_thl[i]
                ql_mean = ql_mean / n_sample

                '''(C) compare to mean profile'''
                print('')
                # find correct times in profile
                time1_ = np.int(files[0][0:-3])
                time2_ = np.int(files[-1][0:-3])
                n = 0
                while(time_profile[n]<nt):
                    while(time_profile[n]<time1_):
                        n+=1
                    n1 = n
                    while(time_profile[n]<time2_):
                        n+=1
                    n2 = n

                if time1_ == time2_:
                    ql_mean_profile = ql_profile[n1,iz]
                else:
                    for n in range(n1,n2):
                        ql_mean_profile += ql_profile[n,iz]
                    ql_mean_profile /= (n2-n1)

                error_array[k,count_ncomp] = ql_mean - ql_mean_field
                if ql_mean_field > 0.0:
                    rel_error_array[k,count_ncomp] = (ql_mean - ql_mean_field) / ql_mean_field
                else:
                    rel_error_array[k,count_ncomp] = 0.0

                print('<ql> from CloudClosure Scheme: ', ql_mean)
                print('<ql> from ql fields: ', ql_mean_field)
                print('<ql> from ql_profile[k]: ', ql_mean_profile)
                print('error (<ql>_CC - <ql>_field): '+str(error_array[k,count_ncomp]))
                print('rel err: '+ str(rel_error_array[k,count_ncomp]))
                # print('times files', time1_, time2_)
                # print('times profile', n1, n2, time_profile[n1], time_profile[n2])
                print('')

                '''(D) plotting'''
                # # samples for ql-shading
                # xmin = np.min([np.amin(data_all[:, 0]),np.amin(Th_l[:, 0])])
                # xmax = np.max([np.amax(data_all[:, 0]),np.amax(Th_l[:, 0])])
                # ymin = np.min([np.amin(data_all[:, 1]),np.amin(Th_l[:, 1])])
                # ymax = np.max([np.amax(data_all[:, 1]),np.amax(Th_l[:, 1])])
                # qt_arr = np.linspace(xmin, xmax, 1e2)
                # thl_arr = np.linspace(ymin, ymax, 1e2)
                # T_comp_thl[i], ql_comp_thl[i], alpha_comp_thl[i] = sat_adj_fromthetali(p_ref[iz], Th_l[i, 0], Th_l[i, 1], CC, LH)

                # plotting
                plot_PDF(data_all, data_all_norm, 'thl', 'qt', nvar, clf_thl_norm, scaler, ncomp, error_array[k,count_ncomp], iz*dz, path_out)
                plot_samples('norm', data_all_norm, ql[:,k], Th_l_norm, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                plot_samples('original', data_all, ql[:,k], Th_l, ql_comp_thl, 'thl', 'qt', scaler, ncomp, iz*dz, path_out)
                plot_hist(ql, path_out)
                print('')
                print('')

                del data_all, data_all_norm
            count_ncomp += 1
        plot_error_vs_ncomp(error_array, rel_error_array, ncomp_range, krange, dz, path_out)
        time2 = time.clock()



        print('total time: ', time2-time1)


        '''(6) Save Gaussian Mixture PDFs '''
        # dump_variable(os.path.join(path_out, 'CloudClosure', nc_file_name_out), 'means', means_, 'qtT', ncomp, nvar, len(zrange))
        # dump_variable(os.path.join(path_out, 'CloudClosure', nc_file_name_out), 'covariances', covariance_, 'qtT', ncomp, nvar, len(zrange))


        return




# #----------------------------------------------------------------------
def relative_entropy(data, clf, kde):
    print('Relative Entropy')
    # rel_ent_clf = D(p_clf || p_kde)
    # rel_ent_kde = D(p_kde || p_clf )

    rel_ent_clf = 0
    rel_ent_kdf = 0
    n_sample = 50
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    X, Y = np.meshgrid(x_, y_)
    XX = np.array([X.ravel(), Y.ravel()]).T
    Z_clf = np.exp(clf.score_samples(XX)).reshape(X.shape)
    Z_kde = np.exp(kde.score_samples(XX)).reshape(X.shape)
    for i in range(n_sample):
        for j in range(n_sample):
            rel_ent_clf += Z_clf[i,j] * np.log(Z_clf[i,j] / Z_kde[i,j])
            rel_ent_kdf += Z_kde[i, j] * np.log(Z_kde[i, j] / Z_clf[i, j])

    print('rel entr D(clf || kdf): ', rel_ent_clf)
    print('rel entr D(kdf || clf): ', rel_ent_kdf)
    # rel_ent_clf_np = np.sum(Z_clf * np.log(Z_clf / Z_kde))
    # rel_ent_kdf_np = np.sum(Z_kde * np.log(Z_kde / Z_clf))

    return
#----------------------------------------------------------------------
cpdef Gaussian_bivariate(ncomp_, data, var_name1, var_name2, time, z, path):
    cdef int ncomp = ncomp_
    clf = mixture.GaussianMixture(n_components=ncomp,covariance_type='full')
    # clf = sklearn.mixture.GaussianMixture(n_components=2, covariance_type='full')
    clf.fit(data)
    # print('')

    # nvar = 2
    # plot_PDF(data, 'thl', 'qt', nvar, clf, path)
    # if var_name1 == 'qt' or var_name2 == 'qt':
    #     plot_PDF_samples_qt(data, var_name1, var_name2, clf, time, z, path)
    #     print('!!!! qt: factor 100')
    # else:
    #     plot_PDF_samples(data, var_name1, var_name2, clf, time, z, path)
    # plot_PDF_samples(data, var_name1, var_name2, clf, time, z)

    return clf

#----------------------------------------------------------------------
cdef Kernel_density_estimate(data, var_name1, var_name2, time_, z, fullpath_out):
    from sklearn.neighbors.kde import KernelDensity
    from plotting_functions import labeling
    ''' Kerne Density Estimation:
    from sklearn.neighbors import KernelDensity

    Parameters:
    - bandwidth: The bandwidth here acts as a smoothing parameter, controlling the tradeoff between bias and variance
    in the result. A large bandwidth leads to a very smooth (i.e. high-bias) density distribution.
    A small bandwidth leads to an unsmooth (i.e. high-variance) density distribution.
    'metric': 'euclidean' (distance metric to use. Note that not all metrics are valid with all algorithms.)
    'atol': 0 (The desired absolute tolerance of the result.)
    'leaf_size': 40
    'kernel': 'gaussian'
    'rtol': 0 (The desired relative tolerance of the result. )
    'breadth_first': True
    'metric_params': None
    'algorithm': 'auto'
    '''

    print('kde: ', data.shape)
    cdef:
        double [:,:] data_aux = data
        int amp = 100
        int ij
        int ij_max = data.shape[0]
        double bw
    # data_aux = np.ndarray(shape=(data.shape))
    # data_aux[:, 0] = data[:, 0]
    print('start loop')
    for ij in range(ij_max):
        data_aux[ij, 1] = data[ij, 1] * amp

    print('end loop')
    # construct a kernel density estimate of the distribution
    # kde = KernelDensity(bandwidth=0.04, metric='haversine',
    #                     kernel='gaussian', algorithm='ball_tree')
    # kde.fit(Xtrain[ytrain == i])
    bw = 5e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    print('after kde')
    # kde.score_samples(data)
    # print('after kde score')

    # Plotting
    n_sample = 100
    x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    X, Y = np.meshgrid(x_, y_)
    XX = np.array([X.ravel(), Y.ravel()]).T
    x_aux = np.linspace(np.amin(data_aux[:, 0]), np.amax(data_aux[:, 0]), n_sample)
    y_aux = np.linspace(np.amin(data_aux[:, 1]), np.amax(data_aux[:, 1]), n_sample)
    X_aux, Y_aux = np.meshgrid(x_aux, y_aux)
    XX_aux = np.array([X_aux.ravel(), Y_aux.ravel()]).T

    print('start figure')
    fig = plt.figure(figsize=(12, 16))
    plt.subplot(3, 2, 1)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    time3 = time.clock()
    print('before score')
    Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    time4 = time.clock()
    # print('after score, time = ', time4-time3)

    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    # ax1 = plt.contour(X_aux, Y_aux, Z_log)
    # plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = '+str(bw))

    plt.subplot(3, 2, 2)
    bw = 3e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
    # Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
    # ax1 = plt.contour(X_aux, Y_aux, Z_log)
    # plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 3)
    bw = 1e-2
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
#     Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
#     ax1 = plt.contour(X_aux, Y_aux, Z_log)
#     plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 4)
    bw = 8e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
#     Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
#     ax1 = plt.contour(X_aux, Y_aux, Z_log)
#     plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 5)
    bw = 5e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
#     Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
#     ax1 = plt.contour(X_aux, Y_aux, Z_log)
#     plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    plt.subplot(3, 2, 6)
    bw = 2e-3
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(data_aux)
    # kde.score_samples(data)
    # Z = np.exp(kde.score_samples(XX_aux)).reshape(X.shape)
#     Z_log = kde.score_samples(XX_aux).reshape(X.shape)
    plt.scatter(data_aux[:, 0], data_aux[:, 1], s=5, alpha=0.2)
#     ax1 = plt.contour(X_aux, Y_aux, Z_log)
#     plt.colorbar(ax1, shrink=0.8)
    labeling(var_name1, var_name2, amp)
    plt.title('bw = ' + str(bw))

    fig.suptitle('Cloud Closure: Kernel Density Estimate (gaussian)', fontsize=20)
    plt.savefig(os.path.join(fullpath_out,'CloudClosure_alltimes_figures','CC_' + var_name1 + '_' + var_name2 + '_z' + str(np.int(z)) + 'm_KDE_alltime.png'))
    plt.close()
#     plt.show()

    print('KDE shapes: ', kde.score_samples(XX).shape, X.shape)
    print(kde.get_params())

    return kde, kde


#----------------------------------------------------------------------
def read_in_fields(group_name, var_list, path):
    rootgrp = nc.Dataset(path, 'r')
    grp = rootgrp.groups[group_name]
    if group_name == 'reference':
        var1 = grp.variables[var_list[0]][:]
        var2 = grp.variables[var_list[1]][:]
        var3 = grp.variables[var_list[2]][:]
        var4 = grp.variables[var_list[3]][:]
        rootgrp.close()
    elif group_name == 'fields':
        var1 = grp.variables[var_list[0]][:,:,:]
        var2 = grp.variables[var_list[1]][:,:,:]
        var3 = grp.variables[var_list[2]][:,:,:]
        var4 = grp.variables[var_list[3]][:,:,:]
        rootgrp.close()
    return var1, var2, var3, var4


def read_in_netcdf(var_name, group_name, path):
    print('read in '+ var_name + ': ' +path)
    rootgrp = nc.Dataset(path, 'r')
    grp = rootgrp.groups[group_name]
    if group_name == 'reference' or group_name == 'timeseries':
        var = grp.variables[var_name][:]
        rootgrp.close()
        return var[:]
    elif group_name == 'profiles':
        var = grp.variables[var_name][:,:]
        rootgrp.close()
        return var[:,:]
    elif group_name == 'fields':
        var = grp.variables[var_name][:,:,:]
        rootgrp.close()
        return var[:,:,:]





#----------------------------------------------------------------------
def create_statistics_file(path,file_name, ncomp, nvar, nz_):
    # ncomp: number of Gaussian components in EM
    # nvar: number of variables of multi-variate Gaussian components
    global time, zrange
    print('create file:', path, file_name)
    rootgrp = nc.Dataset(os.path.join(path,file_name), 'w', format='NETCDF4')
    dimgrp = rootgrp.createGroup('dims')
    means_grp = rootgrp.createGroup('means')
    means_grp.createDimension('nz', nz_)
    means_grp.createDimension('ncomp', ncomp)
    means_grp.createDimension('nvar', nvar)
    cov_grp = rootgrp.createGroup('covariances')
    cov_grp.createDimension('nz', nz_)
    cov_grp.createDimension('ncomp', ncomp)
    cov_grp.createDimension('nvar', nvar)
    weights_grp = rootgrp.createGroup('weights')
    weights_grp.createDimension('nz', nz_)
    weights_grp.createDimension('EM2', 2)
    ts_grp = rootgrp.createGroup('time')
    ts_grp.createDimension('nt',len(time)-1)
    var = ts_grp.createVariable('t','f8',('nt'))
#     for i in range(len(time)-1):
#         var[i] = time[i+1]
#     z_grp = rootgrp.createGroup('z-profile')
#     z_grp.createDimension('nz', len(zrange))
#     var = z_grp.createVariable('height', 'f8', ('nz'))
#     for i in range(len(zrange)):
#         var[i] = zrange[i]
#     rootgrp.close()
#     # print('create file end')
    return

# def dump_variable(path, group_name, data_, var_name, ncomp, nvar, nz_):
#     print('-------- dump variable --------', var_name, group_name, path)
#     # print('dump variable', path, group_name, var_name, data_.shape, ncomp, nvar)
#     if group_name == 'means':
#         add_means(path, var_name, ncomp, nvar)
#         data = np.empty((nz_,ncomp,nvar), dtype=np.double, order='c')
#         for i in range(nz_):
#             for j in range(ncomp):
#                 for k in range(nvar):
#                     data[i,j,k] = data_[i,j,k]
#         write_mean(path, group_name, data, var_name)
#
#     elif group_name == 'covariances':
#         add_covariance(path, var_name, ncomp, nvar)
#         data = np.empty((nz_, ncomp, nvar, nvar), dtype=np.double, order='c')
#         for i in range(nz_):
#             for j in range(ncomp):
#                 for k1 in range(nvar):
#                     for k2 in range(nvar):
#                         data[i, j, k1, k2] = data_[i, j, k1, k2]
#         write_covar(path, group_name, data, var_name)
#
#     elif group_name == 'weights':
#         add_weights(path, var_name, ncomp, nvar)
#         data = np.empty((nz_, ncomp), dtype=np.double, order='c')
#         for i in range(nz_):
#             for j in range(ncomp):
#                 data[i, j] = data_[i, j]
#         write_weights(path, group_name, data, var_name)
#
#     # write_field(path, group_name, data, var_name)
#     # print('--------')
#     return
#
# def add_means(path, var_name, ncomp, nvar):
#     print('add means: ', var_name, path)
#     # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     rootgrp = nc.Dataset(path, 'r+')
#     group = rootgrp.groups['means']
#     var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar'))
#     rootgrp.close()
#     return
#
# def add_covariance(path, var_name, ncomp, nvar):
#     # print('add covariance: ', var_name, path)
#     # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     rootgrp = nc.Dataset(path, 'r+')
#     group = rootgrp.groups['covariances']
#     var = group.createVariable(var_name, 'f8', ('nz', 'ncomp', 'nvar', 'nvar'))
#     rootgrp.close()
#     return
#
# def add_weights(path, var_name, ncomp, nvar):
#     # print('add weights: ', var_name, path)
#     # rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     rootgrp = nc.Dataset(path, 'r+')
#     group = rootgrp.groups['weights']
#     var = group.createVariable(var_name, 'f8', ('nz', 'EM2'))
#     rootgrp.close()
#     return
#
#
# def write_mean(path, group_name, data, var_name):
#     print('write mean:', path, var_name, data.shape)
#     rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     fieldgrp = rootgrp.groups[group_name]
#     var = fieldgrp.variables[var_name]
#     var[:, :, :] = data[:,:,:]
#     rootgrp.close()
#     return
#
# def write_covar(path, group_name, data, var_name):
#     # print('write covar:', path, var_name, data.shape)
#     rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     fieldgrp = rootgrp.groups[group_name]
#     var = fieldgrp.variables[var_name]
#     var[:, :, :,:] = data[:, :, :,:]
#     rootgrp.close()
#     return
#
# def write_weights(path, group_name, data, var_name):
#     # print('write weights:', path, var_name, data.shape)
#     rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     fieldgrp = rootgrp.groups[group_name]
#     var = fieldgrp.variables[var_name]
#     print(var.shape, data.shape)
#     var[:, :] = data[:,:]
#     rootgrp.close()
#     return
#
# def write_field(path, group_name, data, var_name):
#     # print('')
#     # print('write field:', path, var_name, data.shape, group_name)
#     rootgrp = nc.Dataset(path, 'r+', format='NETCDF4')
#     fieldgrp = rootgrp.groups[group_name]
#     # print('fieldgrp', fieldgrp)
#     var = fieldgrp.variables[var_name]
#     # var = data
#        # var[:] = np.array(data)
#     # var[:, :, :] = data
#     rootgrp.close()
#     return
