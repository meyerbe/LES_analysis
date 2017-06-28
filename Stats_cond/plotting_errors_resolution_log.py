# plotting erros for different resolutions
import pylab as plt
import netCDF4 as nc
import sys, os
import json as simplejson
import argparse
import numpy as np

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 21
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 10
plt.rcParams['lines.markeredgewidth'] = 0.0     # the line width around the marker symbol

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("time")
    # parser.add_argument("--Lx")
    parser.add_argument('--Lx', nargs='+', type=int)
    # # This is the correct way to handle accepting multiple arguments (or list)
    # # '+' == 1 or more.
    # # '*' == 0 or more.
    # # '?' == 0 or 1.
    # parser.add_argument("--dz")
    parser.add_argument('--dz', nargs='+', type=int)
    parser.add_argument("--ncomp_max")
    parser.add_argument("--xa_ql")
    parser.add_argument("--xb_ql")
    parser.add_argument("--xa_cf")
    parser.add_argument("--xb_cf")
    parser.add_argument("tracer_type")

    args = parser.parse_args()
    path = args.path


    time_field = np.int(args.time)
    type_ = args.tracer_type
    if args.Lx:
        Lx_range = np.round(np.double(args.Lx),1)
    else:
        Lx_range = [1000.0, 5000.0, 10000.0, 15000.0]
        # Lx_range = [5000.0]
    if args.dz:
        dz_range = args.dz
    else:
        # dz_range = [20, 60, 100]
        dz_range = [20]
    if args.ncomp_max:
        ncomp_max = np.int(args.ncomp_max)
    else:
        ncomp_max = 5
    print('Lx range: ', Lx_range)
    print('dz range: ', dz_range)
    print('ncomp max: ', ncomp_max)


    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    dt_stats = nml['stats_io']['frequency']

    xa = -1.2e-5
    xb = 1e-6
    if args.xa_ql:
        xa = np.double(args.xa_ql)
    if args.xb_ql:
        xb = np.double(args.xb_ql)
    xlimits_ql = [xa, xb]
    xa = -5e-2
    xb = 2e-2
    if args.xa_cf:
        xa = np.double(args.xa_cf)
    if args.xb_cf:
        xb = np.double(args.xb_cf)
    xlimits_cf = [xa, xb]


    ''' Read in zrange from Cloud Closure file '''
    print('')
    print('-- '+type_ + ' --')
    file_name = 'PDF_cond_log_' + type_ + '_allcomp_Lx' + str(Lx_range[0]) + 'Ly' + str(Lx_range[0]) \
                + '_dz' + str(dz_range[0]) + '_time' + str(time_field) + '.nc'
    # file_name = 'PDF_cond_' + type_ + '_Lx10000.0Ly10000.0_dz20' + '_time' + str(time_field) + '.nc'
    path_in = os.path.join(path, 'PDF_cond_log')
    fullpath_in = os.path.join(path_in, file_name)
    print(fullpath_in)
    print('')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    krange_ = rootgrp.groups['profiles'].variables['krange'][:]
    rootgrp.close()

    ''' read in 3D field '''
    root = nc.Dataset(os.path.join(path, 'fields', str(time_field)+'.nc'), 'r')
    ql_3d_ = root.groups['fields'].variables['ql'][:,:,:]
    root.close()

    ''' read in zrange from Stats-file '''
    if time_field > 1e6:
        time_stats = np.int((np.double(time_field)-1e6) / dt_stats)
    else:
        time_stats = np.int(np.double(time_field)/ dt_stats) - 1
    print('time field:        '+ str(time_field))
    print('time Stats file:   ' + str(time_stats) + ' (dt_stats=' + str(dt_stats)+')')
    root = nc.Dataset(os.path.join(path, 'Stats.'+case_name+'.nc'),'r')
    ql_mean_stats_ = root.groups['profiles'].variables['ql_mean'][:]
    cf_stats_ = root.groups['profiles'].variables['cloud_fraction'][:]
    ql_mean_stats = np.zeros(shape=zrange.shape[0])
    cf_stats = np.zeros(shape=zrange.shape[0])
    zrange_stats = root.groups['profiles'].variables['z_half'][:]
    if case_name == 'TRMM_LBA':
        print('extra for TRMM')
        krange = np.asarray([10, 20, 30, 40, 50, 75, 85, 95, 105, 127], dtype=np.int32)
        zrange_stats = root.groups['profiles'].variables['z_half'][:]
    root.close()

    zrange_ = np.zeros(shape=krange_.shape)
    krange = np.zeros(shape=krange_.shape, dtype=np.int)
    ql_3d = np.zeros(shape=(nx, ny, zrange.shape[0]))
    ql_mean_les = np.zeros(shape=krange_.shape, dtype=np.double)
    cf_les = np.zeros(shape=krange_.shape, dtype=np.double)
    for k in range(zrange.shape[0]):
        iz = np.int(krange_[k])
        krange[k] = iz
        zrange_[k] = zrange_stats[iz]
        ql_mean_stats[k] = ql_mean_stats_[time_stats, iz]
        cf_stats[k] = cf_stats_[time_stats, iz]
        ql_3d[:,:,k] = ql_3d_[:,:,iz]
        ql_mean_les[k] = np.average(np.average(ql_3d_[:, :, iz], axis=0), axis=0)
        for i in range(nx):
            for j in range(ny):
                if ql_3d[i,j,k] > 0.0:
                    cf_les[k] += 1.
    cf_les /= (nx*ny)
    if zrange_.any() != zrange.any():
        print('!! problem in zrange !! ')
        sys.exit()
    else:
        del zrange_
    # print('zrange: ', zrange)
    del krange_, ql_3d_
    print('')

    ''' Read in error profiles (from CC_alltime_res_error) and mean profiles (from CC_alltime) '''
    nz = zrange.shape[0]
    error_ql_env = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_ql_rel_env = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_cf_env = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_cf_rel_env = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_ql_domain = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_ql_rel_domain = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_cf_domain = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    error_cf_rel_domain = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    ql_mean_domain = np.zeros((len(Lx_range), len(dz_range), nz))
    ql_mean_env = np.zeros((len(Lx_range), len(dz_range), nz))
    ql_mean_up = np.zeros((len(Lx_range), len(dz_range), nz))
    ql_mean_pdf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))
    cf_domain = np.zeros((len(Lx_range), len(dz_range), nz))
    cf_env = np.zeros((len(Lx_range), len(dz_range), nz))
    cf_up = np.zeros((len(Lx_range), len(dz_range), nz))
    cf_pdf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max, 3))

    # print('shapes: ', error_ql_domain.shape, '#Lx: ', len(Lx_range), '#dz: ', len(dz_range), 'nz: ', nz, 'ncomp_max: ', ncomp_max)
    # n_log_type = 0
    type_list = ['', 'log', 'loglog']
    for n_log_type, log_type in enumerate(type_list):
        print('')
        print('--- PDF cond ' + log_type + ', ' + str(n_log_type)+' ---')
        # print(type(log_type))
        for i in range(len(Lx_range)):
            Lx = Lx_range[i]
            for j in range(len(dz_range)):
                print('')
                delta_z = dz_range[j]
                # # file_name = 'PDF_cond_' + log_type + type_ + '_allcomp_Lx' + str(Lx) + 'Ly' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time_field) + '.nc'
                # file_name = 'PDF_cond_log_' + type_ + '_allcomp_Lx' + str(Lx) + 'Ly' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time_field) + '.nc'
                # # print(file_name)
                # fullpath_in = os.path.join(path_in, file_name)
                # print(fullpath_in)
                # rootgrp = nc.Dataset(fullpath_in, 'r')
                # if log_type != '':
                #     error_ql_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_ql_env' + '_' + log_type][:,
                #                                0:ncomp_max]
                #     error_ql_rel_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables[
                #                                        'rel_error_ql_env' + '_' + log_type][:, 0:ncomp_max]
                #     error_ql_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_ql_domain' + '_' + log_type][
                #                                   :, 0:ncomp_max]
                #     error_ql_rel_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables[
                #                                           'rel_error_ql_domain' + '_' + log_type][:, 0:ncomp_max]
                #     error_cf_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_cf_env' + '_' + log_type][:,
                #                                0:ncomp_max]
                #     error_cf_rel_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables[
                #                                        'rel_error_cf_env' + '_' + log_type][:, 0:ncomp_max]
                #     error_cf_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_cf_domain' + '_' + log_type][
                #                                   :, 0:ncomp_max]
                #     error_cf_rel_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables[
                #                                           'rel_error_cf_domain' + '_' + log_type][:, 0:ncomp_max]
                # else:
                #     error_ql_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_ql_env'][:, 0:ncomp_max]
                #     error_ql_rel_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['rel_error_ql_env'][:, 0:ncomp_max]
                #     error_ql_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_ql_domain'][:, 0:ncomp_max]
                #     error_ql_rel_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['rel_error_ql_domain'][:, 0:ncomp_max]
                #     error_cf_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_cf_env'][:, 0:ncomp_max]
                #     error_cf_rel_env[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_env'][:, 0:ncomp_max]
                #     error_cf_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['error_cf_domain'][:, 0:ncomp_max]
                #     error_cf_rel_domain[i, j, :, :, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_domain'][:, 0:ncomp_max]
                # ql_mean_domain[i,j, :] = rootgrp.groups['profiles'].variables['ql_mean_domain'][:]
                # ql_mean_env[i, j, :] = rootgrp.groups['profiles'].variables['ql_mean_environment'][:]
                # ql_mean_up[i, j, :] = rootgrp.groups['profiles'].variables['ql_mean_updraft'][:]
                # cf_domain[i, j, :] = rootgrp.groups['profiles'].variables['cf_domain'][:]
                # cf_env[i, j, :] = rootgrp.groups['profiles'].variables['cf_environment'][:]
                # cf_up[i, j, :] = rootgrp.groups['profiles'].variables['cf_updraft'][:]
                # zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
                # rootgrp.close()

                for ncomp in range(0,ncomp_max):
                    # file_name = 'CC_alltime_ncomp' +str(ncomp+1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(np.int(Lx)) + '_dz' + str(delta_z) + '_time' + time_field + '.nc'
                    # file_name = 'CC_alltime_anomaly_ncomp' + str(ncomp+1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(np.int(Lx)) + '_dz' + str(
                    #     delta_z) + '_time' + time_field + '.nc'
                    # file_name = 'PDF_cond_' + type_ + '_ncomp' + str(ncomp+1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(np.int(Lx)) + \
                    #         '_dz' + str(dz_range[0]) + '_time' + str(time_field) + '.nc'
                    file_name = 'PDF_cond_log_' + type_ + '_ncomp' + str(ncomp + 1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(
                        np.int(Lx)) + '_dz' + str(dz_range[j]) + '_time' + str(time_field) + '.nc'
                    # print(file_name)
                    fullpath_in = os.path.join(path_in, file_name)
                    print('')
                    print(fullpath_in)
                    print('')
                    rootgrp = nc.Dataset(fullpath_in, 'r')
                    try:
                        ql_mean_pdf[i, j, :, ncomp, n_log_type] = rootgrp.groups['profiles'].variables['ql_mean_pdf'][:]
                        cf_pdf[i, j, :, ncomp, n_log_type] = rootgrp.groups['profiles'].variables['cf_pdf'][:]
                    except:
                        ql_mean_pdf[i, j, :, ncomp, n_log_type] = rootgrp.groups['profiles'].variables['ql_mean_comp'][:]
                        cf_pdf[i, j, :, ncomp, n_log_type] = rootgrp.groups['profiles'].variables['cf_comp'][:]
                    print (ql_mean_pdf[i, j, :, ncomp, n_log_type])
                    if not np.any(ql_mean_pdf[i,j,:,ncomp,n_log_type]):
                        print('is zero ----------', i, j, ncomp, log_type,  ql_mean_pdf[i,j,:,ncomp,n_log_type])

                    if log_type != '':
                        error_ql_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_ql_env' + '_' + log_type][:]
                        error_ql_rel_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_ql_env' + '_' + log_type][:]
                        error_ql_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_ql_domain' + '_' + log_type][:]
                        error_ql_rel_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_ql_domain' + '_' + log_type][:]
                        error_cf_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_cf_env' + '_' + log_type][:]
                        error_cf_rel_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_env' + '_' + log_type][:]
                        error_cf_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_cf_domain' + '_' + log_type][:]
                        error_cf_rel_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_domain' + '_' + log_type][:]
                    else:
                        error_ql_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_ql_env'][:]
                        error_ql_rel_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables[
                                                                       'rel_error_ql_env'][:]
                        error_ql_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_ql_domain'][:]
                        error_ql_rel_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_ql_domain'][:]
                        error_cf_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_cf_env'][:]
                        error_cf_rel_env[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_env'][:]
                        error_cf_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['error_cf_domain'][:]
                        error_cf_rel_domain[i, j, :, ncomp, n_log_type] = rootgrp.groups['error'].variables['rel_error_cf_domain'][:]
                    ql_mean_domain[i, j, :] = rootgrp.groups['profiles'].variables['ql_mean_domain'][:]
                    ql_mean_env[i, j, :] = rootgrp.groups['profiles'].variables['ql_mean_environment'][:]
                    ql_mean_up[i, j, :] = rootgrp.groups['profiles'].variables['ql_mean_updraft'][:]
                    cf_domain[i, j, :] = rootgrp.groups['profiles'].variables['cf_domain'][:]
                    cf_env[i, j, :] = rootgrp.groups['profiles'].variables['cf_environment'][:]
                    cf_up[i, j, :] = rootgrp.groups['profiles'].variables['cf_updraft'][:]
                    zrange_ = rootgrp.groups['profiles'].variables['z'][:]
                    rootgrp.close()

                if zrange_.any() != zrange.any():
                    print('differences in zrange')
                    sys.exit()


        markers = ['o', 'v', '*', 'd', 'p']
        if log_type == '':
            path_out = os.path.join(path, 'PDF_cond_log_figures', 'error_profiles', 'normal')
        else:
            path_out = os.path.join(path, 'PDF_cond_log_figures', 'error_profiles', log_type)
        print(path_out)
        save_name = 'error_ql_env'
        # plot_ql_error_allres_ncompmax(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env,
        #                               error_ql_env[:,:,:,:,n_log_type], error_ql_rel_env[:,:,:,:,n_log_type],
        #                               zrange, ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
        #                               markers, time_field, log_type, path_out, save_name)
        # plot_ql_rel_error_allres_ncompmax(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env,
        #                               error_ql_env[:,:,:,:,n_log_type], error_ql_rel_env[:,:,:,:,n_log_type],
        #                               zrange, ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
        #                                   markers, time_field, log_type, path_out)
        # save_name = 'error_cf_env'
        # plot_cf_error_allres_ncompmax(cf_les, cf_stats, cf_domain, cf_env,
        #                               error_cf_env[:,:,:,:,n_log_type], error_cf_rel_env[:,:,:,:,n_log_type],
        #                               zrange, ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
        #                               markers, time_field, log_type, path_out, save_name)
        # plot_cf_rel_error_allres_ncompmax(cf_les, cf_stats, cf_domain, cf_env,
        #                                   error_cf_env[:,:,:,:,n_log_type], error_cf_rel_env[:,:,:,:,n_log_type],
        #                                   zrange, ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
        #                                   markers, time_field, log_type, path_out)



        n_log_type += 1

    path_out = os.path.join(path, 'PDF_cond_log_figures', 'error_profiles')
    plot_ql_fromsampling_allres(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
                                    ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type,
                                    path_out)
    plot_ql_fromsampling_allres_env(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
                                ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type,
                                path_out)
    # plot_ql_fromsampling_allres_domainmod(0.1, ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
    #                             ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type, path_out)
    plot_cf_fromsampling_allres(cf_les, cf_stats, cf_domain, cf_env, cf_pdf, zrange,
                                ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, log_type, path_out)
    plot_cf_fromsampling_allres_env(cf_les, cf_stats, cf_domain, cf_env, cf_pdf, zrange,
                                ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, log_type, path_out)

    # plotting ql_mean_domain, ql_mean_stats, ql_mean_les, ql_mean_env --> the same for all types
    path_out = os.path.join(path, 'PDF_cond_log_figures', 'error_profiles')
    plot_mean_ql_allres(ql_mean_les, ql_mean_domain, ql_mean_stats, ql_mean_env,
                            zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, path_out)
    plot_cf_allres(cf_les, cf_domain, cf_stats, cf_env,
                            zrange, ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)

    # save_name = 'ql_all_types_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    # plot_ql_all_types(ql_mean_pdf,
    #                   zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range,
    #                   markers, time_field, type_list, save_name, path_out)
    # save_name = 'cf_all_types_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    # plot_cf_all_types(cf_pdf,
    #                   zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range,
    #                   markers, time_field, type_list, save_name, path_out)

    return



def plot_mean_ql_allres(ql_mean_les, ql_mean_domain, ql_mean_stats, ql_mean_env,
                        zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time, path_out):
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')
    # path_out = os.path.join(path, 'error_profiles')
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_dz, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot((10 ** 6) * ql_mean_stats[:], zrange[:], 'k--', label='<ql> from stats-file')
        plt.plot((10 ** 6) * ql_mean_les[:], zrange[:], 'k', label='<ql> from 3D field')
        for nl in range(n_lx):
            col = cm1(np.double(nl) / n_lx)
            plt.plot((10**6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=col, linewidth=3,
                     label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=col, linewidth=3,
                     label='<ql>_env (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        plt.xlim(-6, 100)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'<ql> $\cdot 10^-{6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql>' + '(t=' + str(time) + ')')
    save_name = 'ql_fromLES_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot((10 ** 6) * ql_mean_stats[:], zrange[:], 'k--', label='<ql> from stats-file')
        plt.plot((10 ** 6) * ql_mean_les[:], zrange[:], 'k', label='<ql> from 3D field')
        for ndz in range(n_dz):
            col = cm1(np.double(ndz) / n_dz)
            plt.plot((10**6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=col, linewidth=3, label='<ql> (dz=' + str(dz_range[ndz]) + ')')
            plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=col, linewidth=3,
                     label='<ql>_env (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        plt.xlim(-6, 35)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'<ql> $\cdot 10^-{6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql>' + '(t=' + str(time) + ')')
    save_name = 'ql_fromLES_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()
    return



def plot_cf_allres(cf_les, cf_domain, cf_stats, cf_env,
                        zrange, ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time, path_out):
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')
    # path_out = os.path.join(path, 'error_profiles')
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_dz, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot(cf_stats[:], zrange[:], 'k--', linewidth=3, label='CF from stats-file')
        plt.plot(cf_les[:], zrange[:], 'k', linewidth=3, label='CF from 3D field')
        for nl in range(n_lx):
            col = cm1(np.double(nl) / n_lx)
            plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=col, linewidth=3,
                     label='CF (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=col, linewidth=3,
                     label='CF_env (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        # plt.xlim([-0.25, 0.01])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'CF')
        plt.ylabel('height z (m)')
    plt.suptitle('CF' + '(t=' + str(time) + ')')
    save_name = 'cf_fromLES_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot(cf_stats[:], zrange[:], 'k--', linewidth=3, label='CF from stats-file')
        plt.plot(cf_les[:], zrange[:], 'k', linewidth=3, label='CF from 3D field')
        for ndz in range(n_dz):
            col = cm1(np.double(ndz) / n_dz)
            plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=col, linewidth=3, label='CF (dz=' + str(dz_range[ndz]) + ')')
            plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=col, linewidth=3,
                 label='CF_env (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        # plt.xlim([-0.25, 0.01])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$CF$')
        plt.ylabel('height z (m)')
    plt.suptitle('CF' + '(t=' + str(time) + ')')
    save_name = 'cf_fromLES_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()
    return

def plot_ql_fromsampling_allres_env(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
                                    ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type,
                                    path_out):

    print('')
    print('plot error ql')
    print(ql_mean_pdf.shape)

    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    type_list = ['normal', 'log', 'loglog']
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for nl in range(n_lx):
            print('max ' + type_ + 'Lx=' + str(Lx_range[nl]) + ': ', np.amax(ql_mean_pdf[:, :, :, :, n_type]))
            plt.subplot(3, n_lx, n_type * 3 + nl + 1)
            # # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            # plt.plot((10 ** 6) * ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for ndz in range(n_dz):
                n_col = np.double(ndz) / n_dz
                # plt.plot((10 ** 6) * ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                #          label='<ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    print(np.amax(ql_mean_pdf[nl, ndz, :, nc, n_type]))
                    plt.plot((10 ** 6) * ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz],
                             color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (Lx=' + str(Lx_range[nl]) + 'm')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'ql_sampling_env_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * 3 + ndz + 1)
            # # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            # plt.plot((10 ** 6) * ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for nl in range(n_lx):
                n_col = np.double(nl) / n_lx
                # plt.plot((10 ** 6) * ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                #          label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot((10 ** 6) * ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl],
                             color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'ql_sampling_env_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * 3 + ndz + 1)
            # # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            # plt.plot((10 ** 6) * ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for nl in range(n_lx):
                if Lx_range[nl] > 1000:
                    n_col = np.double(nl) / n_lx
                    # plt.plot((10 ** 6) * ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                    #          label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
                    plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                             label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    for nc in range(ncomp_max):
                        n_col = np.double(nc) / ncomp_max
                        plt.plot((10 ** 6) * ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl],
                                 color=cm1(n_col),
                                 label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'ql_sampling_env_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return

def plot_ql_fromsampling_allres(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
                                  ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type, path_out):
    print('')
    print('plot error ql')
    print(ql_mean_pdf.shape)

    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    type_list = ['normal', 'log', 'loglog']
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for nl in range(n_lx):
            print('max '+type_ + 'Lx='+str(Lx_range[nl]) + ': ', np.amax(ql_mean_pdf[:, :, :, :, n_type]))
            plt.subplot(3, n_lx, n_type*3 + nl + 1)
            # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            plt.plot((10**6)*ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for ndz in range(n_dz):
                n_col = np.double(ndz) / n_dz
                plt.plot((10**6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='<ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot((10**6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    print(np.amax(ql_mean_pdf[nl, ndz, :, nc, n_type]))
                    plt.plot((10**6)*ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (Lx=' + str(Lx_range[nl]) + 'm')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'ql_sampling_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * 3 + ndz + 1)
            # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            plt.plot((10**6)*ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for nl in range(n_lx):
                n_col = np.double(nl) / n_lx
                plt.plot((10**6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot((10**6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot((10**6)*ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'ql_sampling_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * 3 + ndz+ 1)
            # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
            plt.plot((10**6)*ql_mean_les[:], zrange[:], 'b', linewidth=1, label='<ql> (LES)')
            for nl in range(n_lx):
                if Lx_range[nl] > 1000:
                    n_col = np.double(nl) / n_lx
                    plt.plot((10**6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                             label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
                    plt.plot((10**6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                             label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    for nc in range(ncomp_max):
                        n_col = np.double(nc) / ncomp_max
                        plt.plot((10**6)*ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                                 label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc=1)
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel(r'<ql> $\cdot 10^-{6}$')
            plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'ql_sampling_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return

# def plot_ql_fromsampling_allres_domainmod(factor, ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env, ql_mean_pdf, zrange,
#                                     ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, log_type,path_out):
#     print('')
#     # data = (Lx x dz x nz x ncomp)
#
#     import netCDF4 as nc
#
#     cm1 = plt.cm.get_cmap('viridis')
#     cm2 = plt.cm.get_cmap('bone')
#     cm3 = plt.cm.get_cmap('winter')
#     # path_out = os.path.join(path, 'error_profiles')
#     n_lx = len(Lx_range)
#     n_dz = len(dz_range)
#     # print('n_lx', n_lx, 'n_dz', n_dz)
#
#     type_list = ['normal', 'log', 'loglog']
#     plt.figure(figsize=(9 * n_lx, 3*9))
#     for n_type, type_ in enumerate(type_list):
#         for nl in range(n_lx):
#             plt.subplot(3, n_lx, n_type * 3 + nl + 1)
#             # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
#             # plt.plot(factor * (10 ** 6) * ql_mean_les[:], zrange[:], 'b', linewidth=1,
#             #          label=str(factor) + r'$\cdot$<ql> (LES)')
#             for ndz in range(n_dz):
#                 n_col = np.double(ndz) / n_dz
#                 # plt.plot(factor * (10 ** 6) * ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
#                 #          label=str(factor) + r'$\cdot$<ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
#                 # plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
#                 #          label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
#                 for nc in range(1):
#                     n_col = np.double(nc) / ncomp_max
#                     plt.plot((10 ** 6) * ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
#                              label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
#         plt.legend(loc=1)
#         # plt.xlim(0,1.5)
#         plt.title('Lx=' + str(Lx_range[nl]) + 'm')
#         plt.xlabel(r'<ql> $\cdot 10^-{6}$')
#         plt.ylabel('height z (m)')
#     plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
#     save_name = 'ql_sampling_scaled__dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
#     # plt.show()
#     plt.savefig(os.path.join(path_out, save_name))
#     plt.close()
#
#     # (1) one plot per Lx (for all dz)
#     plt.figure(figsize=(9 * n_lx, 3 * 9))
#     for n_type, type_ in enumerate(type_list):
#         for nl in range(n_lx):
#             plt.subplot(3, n_lx, n_type * 3 + nl + 1)
#             # plt.plot((10**6)*ql_mean_stats[:], zrange[:], 'k', label='<ql> (stats)')
#             plt.plot(factor*(10 ** 6) * ql_mean_les[:], zrange[:], 'b', linewidth=1, label=str(factor)+r'$\cdot$<ql> (LES)')
#             for ndz in range(n_dz):
#                 n_col = np.double(ndz) / n_dz
#                 plt.plot(factor*(10 ** 6) * ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
#                          label=str(factor)+r'$\cdot$<ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
#                 plt.plot((10 ** 6) * ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
#                          label='<ql>_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
#                 for nc in range(ncomp_max):
#                     n_col = np.double(nc) / ncomp_max
#                     plt.plot((10 ** 6) * ql_mean_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
#                              label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
#             plt.legend(loc=1)
#             # plt.xlim(0,1.5)
#             plt.title('Lx=' + str(Lx_range[nl]) + 'm')
#             plt.xlabel(r'<ql> $\cdot 10^-{6}$')
#             plt.ylabel('height z (m)')
#     plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
#     save_name = 'ql_sampling_scaled_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
#     # plt.show()
#     plt.savefig(os.path.join(path_out, save_name))
#     plt.close()
#
#     return


def plot_cf_fromsampling_allres_env(cf_les, cf_stats, cf_domain, cf_env, cf_pdf, zrange,
                                ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, log_type, path_out):
    print('')
    print('plot cf from sampling')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 3*10))
    type_list = ['normal', 'log', 'loglog']
    for n_type, type_ in enumerate(type_list):
        print('000', type_, n_type, cf_pdf.shape)
        for nl in range(n_lx):
            plt.subplot(3, n_lx, n_type*n_lx + nl + 1)
            # plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            # plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for ndz in range(n_dz):
                n_col = np.double(ndz) / n_dz
                # plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                #          label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.legend(loc='best')
            # plt.xlim(xlimits_ql)
            plt.title(type_+ ' (Lx=' + str(Lx_range[nl]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_env_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_dz, 3*9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type*n_lx + ndz + 1)
            # # plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            # plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for nl in range(n_lx):
                n_col = np.double(nl) / n_lx
                # plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                #          label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                    plt.legend(loc='best')
            # plt.xlim(xlimits_ql)
            plt.title(type_+' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_env_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type*n_lx + ndz + 1)
            # # plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            # plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for nl in range(n_lx):
                if Lx_range[nl] > 1000:
                    n_col = np.double(nl) / n_lx
                    # plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                    #          label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                             label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    for nc in range(ncomp_max):
                        n_col = np.double(nc) / ncomp_max
                        plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                                 label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc='best')
    #         # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_env_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return

def plot_cf_fromsampling_allres(cf_les, cf_stats, cf_domain, cf_env, cf_pdf, zrange,
                                    ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, log_type,
                                    path_out):


    print('')
    print('plot cf from sampling')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 3 * 10))
    type_list = ['normal', 'log', 'loglog']
    for n_type, type_ in enumerate(type_list):
        for nl in range(n_lx):
            plt.subplot(3, n_lx, n_type * n_lx + nl + 1)
            plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for ndz in range(n_dz):
                n_col = np.double(ndz) / n_dz
                plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.legend(loc='best')
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (Lx=' + str(Lx_range[nl]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_dz, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * n_lx + ndz + 1)
            # plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for nl in range(n_lx):
                n_col = np.double(nl) / n_lx
                plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                    plt.legend(loc='best')
            # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 3 * 9))
    for n_type, type_ in enumerate(type_list):
        for ndz in range(n_dz):
            plt.subplot(3, n_dz, n_type * n_lx + ndz + 1)
            # plt.plot(cf_stats[:], zrange[:], 'k', label='CF (stats)')
            plt.plot(cf_les[:], zrange[:], 'b', label='CF (LES)')
            for nl in range(n_lx):
                if Lx_range[nl] > 1000:
                    n_col = np.double(nl) / n_lx
                    plt.plot(cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                             label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    plt.plot(cf_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                             label='CF_env (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
                    for nc in range(ncomp_max):
                        n_col = np.double(nc) / ncomp_max
                        plt.plot(cf_pdf[nl, ndz, :, nc, n_type], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                                 label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.legend(loc='best')
            #         # plt.xlim(xlimits_ql)
            plt.title(type_ + ' (dz=' + str(dz_range[ndz]) + 'm)')
            plt.xlabel('CF')
            plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + '), Model: ' + log_type)
    save_name = 'cf_sampling_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


def plot_ql_error_allres_ncompmax(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env,
                                  error_ql_env, error_ql_rel_env,
                                  zrange,ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
                                  markers, time, log_type, path_out, save_name_):
    print('')
    print('plot error ql')
    # data = (Lx x dz x nz x ncomp)

    # import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # print('....')
    # print(error_ql_env.shape)
    # print(ql_mean_env.shape)

    # # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9*n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl+1)
        # # plt.plot(-0.4*(1e6)*ql_mean_stats[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (stats)')
        # plt.plot(-0.4*(1e6)*ql_mean_les[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (3D fields)')
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot((1e6) *error_ql_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            n_col = np.double(ndz) / n_dz
            # plt.plot(-(1e6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- <ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
            # plt.plot(-(1e6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1, label='- <ql> env.')
        plt.legend(loc=2)
        # plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)\cdot 10^{-6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + '), Model: ' + log_type)
    # save_name = 'error_ql_allres_dz_nc' + str(ncomp_max)+ '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()



    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_dz, 10))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        # plt.plot(-0.4*(1e6)*ql_mean_stats[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (stats)')
        plt.plot(-0.4 * (1e6) * ql_mean_les[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (3D fields)')
        for nl in range(n_lx):
            n_col = np.double(nl) / n_lx
            plt.plot(-(1e6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='- <ql> (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot(-(1e6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1, label='- <ql> env.')
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot((1e6)*error_ql_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
        plt.legend(loc=2)
    #     plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)\cdot 10^{-6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + '), Model: ' + log_type)
    # save_name = 'error_ql_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 10))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        # plt.plot(-0.4*(1e6)*ql_mean_stats[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (stats)')
        plt.plot(-0.4 * (1e6) * ql_mean_les[:], zrange[:], 'k', linewidth=1, label=r'- 0.4$\cdot$ <ql> (3D fields)')
        for nl in range(n_lx):
            if Lx_range[nl]>1000:
                n_col = np.double(nl) / n_lx
                plt.plot(-(1e6)*ql_mean_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='- <ql> (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot(-(1e6)*ql_mean_env[nl, ndz, :], zrange[:], '-', color=cm2(n_col), linewidth=1,
                         label='- <ql> env.')
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot((1e6)*error_ql_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
        plt.legend(loc=2)
        # plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)\cdot 10^{-6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + '), Model: ' + log_type)
    # save_name = 'error_ql_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return




def plot_ql_rel_error_allres_ncompmax(ql_mean_les, ql_mean_stats, ql_mean_domain, ql_mean_env,
                                  error_ql_env, error_ql_rel_env,
                                  zrange,ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
                                markers, time, log_type, path_out):
    print('')
    print('plot relative error ql')
    # data = (Lx x dz x nz x ncomp)
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')
    # path_out = os.path.join(path, 'error_profiles')
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='$\epsilon(ql)/<ql>$, ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm3(np.double(nc) / ncomp_max)
                plt.plot(error_ql_env[nl, ndz, :, nc] / ql_mean_domain[nl, ndz, :], zrange[:], '-', linewidth=1, marker=markers[ndz], color=col,
                         label=r'$\epsilon(ql)/<ql>$_domain, ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
                # plt.plot(ql_mean_domain[nl, ndz, :], zrange[:], '-')
            # plt.plot(-ql_mean_domain[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.xlim([-1.1, 1.1])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>) / <ql>$_env')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_rel_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                # plt.plot(error_ql_rel_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                #          label='$\epsilon(ql)$, ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
                # col = cm3(np.double(nc) / ncomp_max)
                plt.plot(error_ql_env[nl, ndz, :, nc] / ql_mean_domain[nl, ndz, :], zrange[:], '-',
                         marker=markers[ndz], color=col,
                         label=r'$\epsilon(ql)/<ql>$_domain, ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
                # plt.plot(ql_mean_domain[nl, ndz, :], zrange[:], '-')
                # plt.plot(-ql_mean_domain[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        # lim = np.amax(ql_mean_env[nl, :, 2:] / ql_mean_domain[2:])
        plt.xlim([-0.2,0.2])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>) / <ql>$_domain')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + '), Model: ' + log_type)
    save_name = 'error_ql_rel_domain_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()




    # # (2) one plot per dz (for all Lx)
    # plt.figure(figsize=(9 * n_dz, 12))
    # for ndz in range(n_dz):
    #     plt.subplot(1, n_dz, ndz + 1)
    #     plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
    #     for nl in range(n_lx):
    #         for nc in range(ncomp_max):
    #             col = cm1(np.double(nc) / ncomp_max)
    #             plt.plot(error_ql_rel_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
    #                      label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
    #         # plt.plot(-ql_mean_domain[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
    #     plt.legend(loc=1)
    #     #plt.xlim(xlimits_ql)
    #     plt.title('dz=' + str(dz_range[ndz]) + 'm')
    #     plt.xlabel(r'$\epsilon(<ql>)$')
    #     plt.ylabel('height z (m)')
    # plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    # save_name = 'error_ql_rel_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.savefig(os.path.join(path_out, save_name))
    # # plt.show()
    # plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_ql_rel_env[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-ql_mean_domain[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.xlim([-1.1,1.1])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>) / <ql>$_env')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + '), Model: ' + log_type)
    save_name = 'error_ql_rel_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return




def plot_cf_error_allres_ncompmax(cf_les, cf_stats, cf_domain, cf_env, error_cf, error_cf_rel, zrange,
                                  ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
                                  markers, time, log_type, path_out, save_name_):
    print('')
    print('plot error CF')
    # print('CF: ', cf_domain)
    # print('e(CF): ', error_cf)
    # print('ql: ', ql_mean_field)
    # print('e(ql): ', error_ql)
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)
    # markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot(-0.4*cf_stats, zrange, 'k', linewidth=1, label='-0.4*CF stats')
        plt.plot(-0.4*cf_les, zrange, 'b', linewidth=1, label='-0.4*CF (3D fields)')
        for ndz in range(n_dz):

            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            n_col = np.double(ndz) / n_dz
            plt.plot(-0.4*cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='-0.4*CF (dz=' + str(dz_range[ndz]) + ')')
            plt.plot(-cf_env[nl, ndz, :], zrange[:], '-', linewidth=1, color=cm2(n_col),
                     label='- CF_env (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=2)
        # plt.xlim(xlimits_cf)
        plt.xlim(-0.03,0.01)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<CF>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + ')')
    # save_name = 'error_cf_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_dz, 10))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot(-0.4*cf_stats, zrange, 'k', linewidth=1, label=r'- 0.4$\cdot$CF stats')
        plt.plot(-0.4*cf_les, zrange, 'b', linewidth=1, label='- 0.4*CF (3D fields)')
        for nl in range(n_lx):
            n_col = np.double(nl) / n_lx
            plt.plot(-0.4*cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- CF (Lx=' + str(Lx_range[nl]) + ')')
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
        plt.legend(loc=2)
        # plt.xlim(-1.1,0.3)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<CF>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + '), Model: ' + log_type)
    # save_name = 'error_cf_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 10))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot(-cf_stats, zrange, 'k', linewidth=1, label='- CF stats')
        plt.plot(-cf_les, zrange, 'b', linewidth=1, label='- CF (3D fields)')
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                n_col = np.double(nl) / n_lx
                plt.plot(-cf_domain[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- CF (Lx=' + str(Lx_range[nl]) + ')')
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
        plt.legend(loc=2)
        plt.xlim(xlimits_cf)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(CF)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + '), Model: ' + log_type)
    # save_name = 'error_cf_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    save_name = save_name_ + '_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


def plot_cf_rel_error_allres_ncompmax(cf_les, cf_stats, cf_domain, cf_env, error_cf_env, error_cf_rel, zrange,
                                      ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range,
                                      markers, time, log_type, path_out):
    print('')
    print('plot relative error CF')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm3(np.double(nc) / ncomp_max)
                plt.plot(error_cf_env[nl, ndz, :, nc] / cf_domain[nl, ndz, :], zrange[:], '-', linewidth=1,
                             marker=markers[ndz], color=col,
                             label=r'$\epsilon(CF)/CF$_domain, ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            # plt.plot(-cf_domain[nl, ndz, :], zrange[:], 'k--', label='- CF (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        plt.xlim([-1.1, 1.1])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(CF) / CF$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error CF' + '(t=' + str(time) + '), Model: ' + log_type)
    save_name = 'error_cf_rel_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # # (2) one plot per dz (for all Lx)
    # plt.figure(figsize=(9 * n_dz, 12))
    # for ndz in range(n_dz):
    #     plt.subplot(1, n_dz, ndz + 1)
    #     for nl in range(n_lx):
    #         for nc in range(ncomp_max):
    #             col = cm1(np.double(nc) / ncomp_max)
    # #             plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
    # #                      label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
    # #         plt.plot(-cf_domain[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
    # #     plt.legend(loc=1)
    #     #plt.xlim(xlimits_cf)
    #     plt.title('dz=' + str(dz_range[ndz]) + 'm')
    #     plt.xlabel(r'$\epsilon(CF)$')
    #     plt.ylabel('height z (m)')
    # plt.suptitle('error CF' + '(t=' + str(time) + ')')
    # save_name = 'error_cf_rel_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.savefig(os.path.join(path_out, save_name))
    # # plt.show()
    # plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_dz, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-cf_domain[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        plt.xlim([-1.1, 1.1])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(CF) / CF$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error CF' + '(t=' + str(time) + '), Model: ' + log_type)
    save_name = 'error_cf_rel_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '_' + log_type + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return



def plot_ql_all_types(ql_mean_pdf,
                      zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range,
                      markers, time_field, type_list, save_name, path_out):
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_lx, 10))
    ncomp = 1
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            n_col = np.double(ndz) / n_dz
            for type in range(len(type_list)):
                plt.plot((1e6)*ql_mean_pdf[nl, ndz, :, ncomp, type], zrange,
                         label='type: '+type_list[type] + ' (dz='+str(dz_range[ndz])+')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'<ql> $\cdot 10^-{6}$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    # save_name = 'ql_all_types_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    return

def plot_cf_all_types(cf_pdf,
                       zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range,
                       markers, time_field, type_list, save_name, path_out):


    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_lx, 10))
    ncomp = 1
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            n_col = np.double(ndz) / n_dz
            for type in range(len(type_list)):
                plt.plot((1e6) * cf_pdf[nl, ndz, :, ncomp, type], zrange,
                         label='type: ' + type_list[type] + ' (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'CF')
        plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + ')')
    # save_name = 'cf_all_types_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    return


if __name__ == '__main__':
    main()