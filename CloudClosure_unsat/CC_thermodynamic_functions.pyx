import os
include 'parameters.pxi'
import numpy as np
cimport numpy as np
import netCDF4 as nc
import json as simplejson

cpdef do_everything(path, path_ref):
    path_fields = os.path.join(path, 'fields', )

    # (1) Reference State
    try:
        p_ref = read_in_netcdf('p0_half', 'reference', path_ref)
    except:
        print('no p0_half profile')
        p_ref = read_in_netcdf('p0', 'reference', path_ref)

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
    s_ = read_in_netcdf('s', 'fields', path_fields)
    qt_ = read_in_netcdf('qt', 'fields', path_fields)
    T_ = read_in_netcdf('temperature', 'fields', path_fields)
    ql_ = read_in_netcdf('ql', 'fields', path_fields)
    print('')

    cdef:
        int i, j, k
        double s, qt, T, ql
        double pv, pd
        double T_unsat
    i = 10
    j = 10
    for k in [10]:
        s = s_[i,j,k]
        qt = qt_[i,j,k]
        ql = ql_[i,j,k]
        T = T_[i,j,k]
        pref = p_ref[k]
        qv = np.double(qt, copy = True)

        pv = pref * eps_vi * qv / (1.0 - qt + eps_vi * qv)
        pd = pref - pv

        T_unsat = T_tilde * np.exp((s -
                                    (1.0 - qt)*(sd_tilde - Rd*np.log(pd/p_tilde))
                                    - qt*(sv_tilde - Rv*np.log(pv/p_tilde)) )
                                    / ( (1.0-qt)*cpd + qt*cpv))

        print('difference: ', T - T_unsat)
        print('ql=', ql)
        print('')
        print('')

    return



cpdef do_everything_with_pycles(path, path_ref):
    cdef extern from "thermodynamics_sa.h":
        inline double temperature_no_ql(double pd, double pv, double s, double qt)
    cdef extern from "thermodynamic_functions.h":
        inline double pd_c(const double p0,const double qt, const double qv)
        inline double pv_c(const double p0, const double qt, const double qv)

    path_fields = os.path.join(path, 'fields', )

    # (0) Namelist File
    nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    gw = nml['grid']['gw']


    # (1) Reference State
    # try:
    #     p_ref = read_in_netcdf('p0_half', 'reference', path_ref)
    # except:
    #     print('no p0_half profile')
    #     p_ref = read_in_netcdf('p0', 'reference', path_ref)
    p_ref = read_in_netcdf('p0', 'reference', path_ref)
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
    s_ = read_in_netcdf('s', 'fields', path_fields)
    qt_ = read_in_netcdf('qt', 'fields', path_fields)
    T_ = read_in_netcdf('temperature', 'fields', path_fields)
    ql_ = read_in_netcdf('ql', 'fields', path_fields)
    print('')

    cdef:
        int i, j, k
        double s, qt, T, ql, qv
        double pv, pd
        double T_unsat
        double T_max = -9999.9
    i = 10
    j = 10

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
    # for k in [10,20,30]:
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
                    if np.abs(T_unsat - T) > T_max:
                        T_max = np.abs(T_unsat - T)

    print('difference: ', T_max)
    # print('difference: ', T - T_unsat)
    # print('ql=', ql)
    # print('T', T, s)
    print('')
    print('')



    # p0 from Stats-file (k = 0,...,gw-1,gw)
    pref_ = np.ones(10)
    pref_[0] = 101045.425934746
    pref_[1] = 100141.496848432
    pref_[2] = 99244.4861572627
    pref_[3] = 98354.3517151152
    pref_[4] = 97471.0459866453
    pref_[5] = 96594.5278115847
    pref_[6] = 95724.7506264438
    pref_[7] = 94861.6724312234

    ########### output at ijk=0 #################
    s = 6960.2884872
    qt = 0.0160628857142
    T = 295.153291065
    qv = qt
    # pref = 95724.7506264    # p0_half[0]
    # pref = 96158.7990115    # p0[0]
    for k in range(gw+1):
        pref = pref_[k]
        print('pref=p0['+str(k)+']  '+str(pref))

        pv = pv_c(pref, qt, qv)
        pd = pd_c(pref, qt, qv)
        T_unsat = temperature_no_ql(pd, pv, s, qt)
        print('directly from LES: ijk=0')
        print('difference: ', T - T_unsat, T_unsat)
        print('ql=', ql)
        print('')
    print('')

    ########### output at ijk=imin*istride+jmin*istride+kmin #################
    s = 6964.59308601
    qt = 0.01696061616
    T =  299.593565755
    qv = qt
    # pref = 95724.7506264    # p0_half[0]
    # pref = 96158.7990115    # p0[0]
    for k in range(gw+1):
        pref = pref_[k]
        print('pref=p0['+str(k)+']  '+str(pref))

        pv = pv_c(pref, qt, qv)
        # pd = pd_c(pref, qt, qv)
        pd = pref - pv
        T_unsat = temperature_no_ql(pd, pv, s, qt)
        print('directly from LES: ')
        print('difference: ', T - T_unsat, T_unsat)
        print('ql=', ql)
        print('')
    print('')

    # pref = 96158.7990115    # p0[0]
    # # pref = 95724.7506264    # p0_half[0]
    # pv = pv_c(pref, qt, qv)
    # # pd = pd_c(pref, qt, qv)
    # pd = pref - pv
    # T_unsat = temperature_no_ql(pd, pv, s, qt)
    # print('directly from LES: ')
    # print('difference: ', T - T_unsat, T_unsat)
    # print('ql=', ql)
    # print('')

    ########### output at ijk=imin*istride+jmin*istride+(kmin+1) #################
    s = 6963.49595557
    qt = 0.0167766796152
    qv = qt
    T = 298.79515229
    for k in range(gw+1):
        pref = pref_[k]
        print('pref=p0['+str(k)+']  '+str(pref))

        pv = pv_c(pref, qt, qv)
        # pd = pd_c(pref, qt, qv)
        pd = pref - pv
        T_unsat = temperature_no_ql(pd, pv, s, qt)
        print('directly from LES: ')
        print('difference: ', T - T_unsat, T_unsat)
        print('ql=', ql)
        print('')
    print('')

    return




def read_in_netcdf(var_name, group_name, path):
    print('read in '+ var_name + ': ' +path)
    rootgrp = nc.Dataset(path, 'r')
    grp = rootgrp.groups[group_name]
    # print(grp)
    #var = grp.variables[var_name]
    #shape = var.shape
    #print(var_name)
    # print('shape:',var.shape)
    # data = np.ndarray(shape=var.shape)
    # data = var[:]
    # print('shape', data.shape)

    #print(var)
    if group_name == 'reference':
        var = grp.variables[var_name][:]
        rootgrp.close()
        return var
    elif group_name == 'fields':
        var = grp.variables[var_name][:,:,:]
        rootgrp.close()
        return var

