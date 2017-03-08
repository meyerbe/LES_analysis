import os
include 'parameters.pxi'
import numpy as np
cimport numpy as np
import netCDF4 as nc
import json as simplejson
import time

from cpython cimport array
import array
import cython
from cython cimport floating


import CC_thermodynamics_c
from CC_thermodynamics_c import sat_adj_fromentropy_c

cdef class CloudClosure:
    def __init__(self):
        print('Cloud Closure')
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

    cpdef initialize(self, path, path_ref):
        nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
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
        print('p_ref', np.shape(self.p_ref), self.p_ref.shape, type(self.p_ref[0]), self.p_ref[0])
        print('')

        return



    cpdef verification_CC(self, path, path_ref):
        print('Verification')
        cdef extern from "thermodynamics_sa.h":
            inline double temperature_no_ql(double pd, double pv, double s, double qt)
        cdef extern from "thermodynamic_functions.h":
            inline double pd_c(const double p0,const double qt, const double qv)
            inline double pv_c(const double p0, const double qt, const double qv)

        path_fields = os.path.join(path, 'fields', )
        cdef:
            int i, j, k

        # (0) Namelist File
        nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
        # cdef int nx, ny, nz, gw
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
        files = os.listdir(os.path.join(path, 'fields'))
        N = len(files)
        print('Found the following directories', files, N)
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
            double pv, pd
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
        LH = CC_thermodynamics_c.LatentHeat(nml)
        print('done Lookup table')
        CC_c = CC_thermodynamics_c.ClausiusClapeyron_c()
        print('done CC __init__')
        CC_c.initialize(nml, LH)
        print('done CC')
        # CC = CC_thermodynamics_c.ClausiusClapeyron()
        # CC.initialize(nml, LH)
        alpha_ = np.zeros(shape=(nx,ny,nz))
        T_comp = np.zeros(shape=(nx,ny,nz))
        ql_comp = np.zeros(shape=(nx,ny,nz))
        alpha = np.zeros(shape=(nx,ny,nz))
        theta_l = np.zeros(shape=(nx,ny,nz))
        T_comp_thl = np.zeros(shape=(nx,ny,nz))
        ql_comp_thl = np.zeros(shape=(nx,ny,nz))
        alpha_thl = np.zeros(shape=(nx,ny,nz))
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
                    one, two, three = sat_adj_fromentropy_c(p_ref[k], s, qt, CC_c, LH)
                    T_comp[i,j,k], ql_comp[i,j,k], alpha[i,j,k] = sat_adj_fromentropy_c(p_ref[k], s, qt, CC_c, LH)
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
        print('From Entropy (CC_c):')
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




    cpdef predict_pdf(self, path, path_ref, ncomp_):
        print('Predict PDF')
        # global ncomp
        # global nvar
        cdef:
            int ncomp = ncomp_
            int nvar = 2
        # create_statistics_file(os.path.join(fullpath_out, 'CloudClosure'), nc_file_name_out, ncomp, nvar, len(zrange))

        # for i in range(len(zrange)):
        #     iz = zrange[i]
        #     for d in files:
        #         '''(1) compute liquid potential temperature from temperature and moisture'''
        #         p0 = 1e5
        #         # T, ql, qi = sat_adj(p0, 6500, 1e-3)
        #
        #         nc_file_name = str(d)
        #         fullpath_in = os.path.join(in_path, 'fields', nc_file_name)
        #         print('fullpath_in', fullpath_in)
        #         T = read_in_netcdf('temperature', 'fields', fullpath_in)
        #         qt = read_in_netcdf('qt', 'fields', fullpath_in)
        #         ql = read_in_netcdf('ql', 'fields', fullpath_in)
        #         qi = np.zeros(shape=T.shape)
        #         theta_l = thetali(p0,T,qt,ql,qi)
        #
        #         data = np.ndarray(shape=((nx * ny), nvar))
        #         means_ = np.ndarray(shape=(len(zrange), ncomp, nvar))
        #         covariance_ = np.zeros(shape=(len(zrange), ncomp, nvar, nvar))
        #
        #         data1_ = theta_l.reshape((nx * ny), nz)
        #         data2_ = qt.reshape((nx * ny), nz)
        #         data[:, 0] = data1_[:, iz]
        #         data[:, 1] = data2_[:, iz]
        #         data_all = np.append(data_all, data, axis=0)
        #
        #     '''(2) Compute bivariate Gaussian PDF (theta_l, qt) '''
        #     # means, covariance, weights = Gaussian_mixture_bivariate(data, var1, var2, np.int(d[0:-3]), iz*dz)
        #     clf = Gaussian_bivariate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz)
        #     means_[i, :, :] = clf.means_[:, :]
        #     covariance_[i,:,:,:] = clf.covariances_[:,:,:]
        #
        #     '''(3) Compute Kernel-Estimate PDF '''
        #     kde, kde_aux = Kernel_density_estimate(data, 'T', 'qt', np.int(d[0:-3]), iz * dz)
        #
        #     relative_entropy(data, clf, kde)
        #
        #     '''(4) Save Gaussian Mixture PDFs '''
        # dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name_out), 'means', means_, 'qtT', ncomp, nvar, len(zrange))
        # dump_variable(os.path.join(fullpath_out, 'CloudClosure', nc_file_name_out), 'covariances', covariance_, 'qtT', ncomp, nvar, len(zrange))
        return





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
    if group_name == 'reference':
        var = grp.variables[var_name][:]

        rootgrp.close()
        return var[:]
    elif group_name == 'fields':
        var = grp.variables[var_name][:,:,:]
        rootgrp.close()
        return var[:,:,:]

