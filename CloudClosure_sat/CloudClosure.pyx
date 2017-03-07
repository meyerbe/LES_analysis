import os
include 'parameters.pxi'
import numpy as np
cimport numpy as np
import netCDF4 as nc
import json as simplejson
import time


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

    cpdef verification_CC(self, path, path_ref):
        cdef extern from "thermodynamics_sa.h":
            inline double temperature_no_ql(double pd, double pv, double s, double qt)
        cdef extern from "thermodynamic_functions.h":
            inline double pd_c(const double p0,const double qt, const double qv)
            inline double pv_c(const double p0, const double qt, const double qv)

        path_fields = os.path.join(path, 'fields', )

        # (0) Namelist File
        nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
        cdef int nx, ny, nz, gw
        nx = nml['grid']['nx']
        ny = nml['grid']['ny']
        nz = nml['grid']['nz']
        gw = nml['grid']['gw']

        # (1) Reference State
        try:
            p_ref = read_in_netcdf('p0_half', 'reference', path_ref)
        except:
            print('no p0_half profile')
            p_ref = read_in_netcdf('p0', 'reference', path_ref)
        # p_ref = read_in_netcdf('p0', 'reference', path_ref)
        print('p_ref: ', p_ref.shape, nz, gw)

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
        cdef:
            int i, j, k
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

    cpdef do_everything_with_pycles(self, path, path_ref):

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
    elif group_name == 'fields':
        var = grp.variables[var_name][:,:,:]
        rootgrp.close()
    return var

