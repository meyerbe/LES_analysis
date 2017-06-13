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
plt.rcParams['legend.fontsize'] = 10
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

    args = parser.parse_args()
    path = args.path


    time_field = args.time
    if args.Lx:
        Lx_range = np.round(np.double(args.Lx),1)
    else:
        Lx_range = [1000.0, 5000.0, 10000.0, 15000.0]
        # Lx_range = [5000.0]
    if args.dz:
        dz_range = args.dz
    else:
        # dz_range = [20, 60, 100]
        dz_range = [40]
    if args.ncomp_max:
        ncomp_max = np.int(args.ncomp_max)
    else:
        ncomp_max = 3
    print('Lx range: ', Lx_range)
    print('dz range: ', dz_range)
    print('ncomp max: ', ncomp_max)


    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
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
    file_name = 'CC_alltime_anomaly_res_error' + '_Lx' + str(Lx_range[0]) + 'Ly' + str(Lx_range[0]) \
                + '_dz' + str(dz_range[0]) + '_time' + time_field + '.nc'
    # file_name = 'CC_alltime_res_error' + '_Lx10000.0Ly10000.0_dz20' + '_time' + time_field + '.nc'
    path_in = os.path.join(path, 'CloudClosure_res_anomaly_5h')
    # path_in = os.path.join(path, 'CloudClosure_res')
    fullpath_in = os.path.join(path_in, file_name)
    print(path_in)
    print('')
    rootgrp = nc.Dataset(fullpath_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    rootgrp.close()

    ''' read in zrange from Stats-file '''
    if time_field > 1e6:
        time_stats = np.int((np.double(time_field)-1e6) / dt_stats)
    else:
        time_stats = np.int(np.double(time_field)/ dt_stats)
    print('time field: '+ str(time_field))
    print('time Stats file: ' + str(time_stats) + ' (dt_stats=' + str(dt_stats)+')')
    root = nc.Dataset(os.path.join(path, 'Stats.'+case_name+'.nc'),'r')
    ql_mean_stats_ = root.groups['profiles'].variables['ql_mean'][:]
    cf_stats_ = root.groups['profiles'].variables['cloud_fraction'][:]
    ql_mean_stats = np.zeros(shape=zrange.shape[0])
    cf_stats = np.zeros(shape=zrange.shape[0])
    if case_name == 'TRMM_LBA':
        krange = np.asarray([10, 20, 30, 40, 50, 75, 85, 95, 105, 127], dtype=np.int32)
        zrange_stats = root.groups['profiles'].variables['z_half'][:]
        for k in range(zrange.shape[0]):
            zrange[k] = zrange_stats[krange[k]]
            ql_mean_stats[k] = ql_mean_stats_[time_stats, krange[k]]
            cf_stats[k] = cf_stats_[time_stats, krange[k]]
    root.close()
    print('zrange: ', zrange)
    print('')

    ''' Read in error profiles (from CC_alltime_res_error) and mean profiles (from CC_alltime) '''
    nz = zrange.shape[0]
    error_ql = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_ql_rel = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_cf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_cf_rel = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    ql_mean_field = np.zeros((len(Lx_range), len(dz_range), nz))
    ql_mean_pdf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    cf_field = np.zeros((len(Lx_range), len(dz_range), nz))
    cf_pdf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))

    print('shapes: ', error_ql.shape, '#Lx: ', len(Lx_range), '#dz: ', len(dz_range), 'nz: ', nz, 'ncomp_max: ', ncomp_max)
    for i in range(len(Lx_range)):
        Lx = Lx_range[i]
        for j in range(len(dz_range)):
            print('')
            delta_z = dz_range[j]
            # file_name = 'CC_alltime_res_error' + '_Lx' + str(Lx) + 'Ly' + str(Lx) + '_dz' + str(delta_z) + '_time' + time_field + '.nc'
            file_name = 'CC_alltime_anomaly_res_error' + '_Lx' + str(Lx) + 'Ly' + str(Lx) + '_dz' + str(
                delta_z) + '_time' + time_field + '.nc'
            print(file_name)
            fullpath_in = os.path.join(path_in, file_name)
            print(fullpath_in)
            rootgrp = nc.Dataset(fullpath_in, 'r')
            error_ql[i, j, :, :] = rootgrp.groups['error'].variables['error_ql'][:, 0:ncomp_max]
            error_ql_rel[i, j, :, :] = rootgrp.groups['error'].variables['rel_error_ql'][:, 0:ncomp_max]
            ql_mean_field[i,j, :] = rootgrp.groups['profiles'].variables['ql_mean_field'][:]
            cf_field[i, j, :] = rootgrp.groups['profiles'].variables['cf_field'][:]
            error_cf[i, j, :, :] = rootgrp.groups['error'].variables['error_cf'][:, 0:ncomp_max]
            error_cf_rel[i, j, :, :] = rootgrp.groups['error'].variables['rel_error_cf'][:, 0:ncomp_max]
            zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
            rootgrp.close()

            for ncomp in range(0,ncomp_max):
                # file_name = 'CC_alltime_ncomp' +str(ncomp+1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(np.int(Lx)) + '_dz' + str(delta_z) + '_time' + time_field + '.nc'
                file_name = 'CC_alltime_anomaly_ncomp' + str(ncomp+1) + '_Lx' + str(np.int(Lx)) + 'Ly' + str(np.int(Lx)) + '_dz' + str(
                    delta_z) + '_time' + time_field + '.nc'
                # print(file_name)
                fullpath_in = os.path.join(path_in, file_name)
                print('')
                print(fullpath_in)
                print('')
                rootgrp = nc.Dataset(fullpath_in, 'r')
                # ql_mean_pdf[i, j, :, ncomp] = rootgrp.groups['profiles'].variables['ql_mean_pdf'][:]
                ql_mean_pdf[i, j, :, ncomp] = rootgrp.groups['profiles'].variables['ql_mean_comp'][:]
                # cf_pdf[i, j, :, ncomp] = rootgrp.groups['profiles'].variables['cf_pdf'][:]
                cf_pdf[i, j, :, ncomp] = rootgrp.groups['profiles'].variables['cf_comp'][:]
                rootgrp.close()
            if zrange_.any() != zrange.any():
                print('differences in zrange')
                sys.exit()

    markers = ['o', 'v', '*', 'd', 'p']
    # path_out = os.path.join(path, 'CloudClosure_res_figures')
    path_out = os.path.join(path, 'CloudClosure_res_anomaly_figures_5h')
    # plot_ql_error_allres_ncompmax(ql_mean_stats, ql_mean_field, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
    #                     ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)
    # plot_ql_rel_error_allres_ncompmax(ql_mean_stats, ql_mean_field, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
    #                            ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)
    # plot_cf_error_allres_ncompmax(cf_stats, cf_field, error_cf, error_cf_rel, zrange,
    #                               ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)
    # plot_cf_rel_error_allres_ncompmax(cf_stats, cf_field, error_cf, error_cf_rel, zrange,
    #                               ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)

    plot_ql_fromsampling_allres(ql_mean_stats, ql_mean_field, ql_mean_pdf, zrange,
                                  ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, path_out)
    plot_cf_fromsampling_allres(cf_stats, cf_field, cf_pdf, zrange,
                                ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)

    plot_mean_ql_allres(ql_mean_field, ql_mean_stats, zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, path_out)
    plot_mean_cf_allres(cf_field, cf_stats, zrange, ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path_out)

    return



def plot_mean_ql_allres(ql_mean_field, ql_mean_stats, zrange, ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time, path):
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')
    path_out = os.path.join(path, 'error_profiles')
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_dz, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            col = cm1(np.double(nl) / n_lx)
            plt.plot((10**6)*ql_mean_field[nl, ndz, :], zrange[:], '-', color=col, linewidth=3, label='- <ql> (Lx=' + str(Lx_range[nl]) + ')')
        plt.plot((10 ** 6) * ql_mean_stats[:], zrange[:], 'k--', label='- <ql> from stats-file')
        plt.legend(loc=1)
        plt.xlim(-6, 35)
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
        for ndz in range(n_dz):
            col = cm1(np.double(ndz) / n_dz)
            plt.plot((10**6)*ql_mean_field[nl, ndz, :], zrange[:], '-', color=col, linewidth=3, label='- <ql> (dz=' + str(dz_range[ndz]) + ')')
        plt.plot((10 ** 6) * ql_mean_stats[:], zrange[:], 'k--', label='- <ql> from stats-file')
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



def plot_mean_cf_allres(cf_field, cf_stats, zrange, ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time, path):
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')
    path_out = os.path.join(path, 'error_profiles')
    n_lx = len(Lx_range)
    n_dz = len(dz_range)

    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            col = cm1(np.double(nl) / n_lx)
            plt.plot(cf_field[nl, ndz, :], zrange[:], '-', color=col, linewidth=3, label='- CF (Lx=' + str(Lx_range[nl]) + ')')
        plt.plot(cf_stats[:], zrange[:], 'k--', linewidth=3, label='- CF from stats-file')
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
        for ndz in range(n_dz):
            col = cm1(np.double(ndz) / n_dz)
            plt.plot(cf_field[nl, ndz, :], zrange[:], '-', color=col, linewidth=3, label='- CF (dz=' + str(dz_range[ndz]) + ')')
        plt.plot(cf_stats[:], zrange[:], 'k--', linewidth=3, label='- CF from stats-file')
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





def plot_ql_fromsampling_allres(ql_mean_stats, ql_mean_field, ql_mean_pdf, zrange,
                                  ncomp_max, xlimits_ql, dz, Lx_range, dz_range, markers, time_field, path):
    print('')
    print('plot error ql')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(ql_mean_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='<ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
            plt.plot(ql_mean_stats[:], zrange[:], 'k--', label='<ql> (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$<ql>$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'ql_sampling_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(ql_mean_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot(ql_mean_stats[:], zrange[:], 'k--', label='<ql> (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$<ql>$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'ql_sampling_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(ql_mean_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='<ql> (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot(ql_mean_stats[:], zrange[:], 'k--', label='<ql> (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$<ql>$')
        plt.ylabel('height z (m)')
    plt.suptitle('<ql> from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'ql_sampling_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


def plot_cf_fromsampling_allres(cf_stats, cf_field, cf_pdf, zrange,
                                ncomp_max, xlimits_cf, dz, Lx_range, dz_range, markers, time_field, path):
    print('')
    print('plot cf from sampling')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(cf_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(cf_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='CF (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
            plt.plot(cf_stats[:], zrange[:], 'k--', label='CF (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel('CF')
        plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'cf_sampling_dz_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(cf_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(cf_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                     label='CF (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot(cf_stats[:], zrange[:], 'k--', label='CF (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel('CF')
        plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'cf_sampling_Lx__nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(cf_pdf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(cf_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col),
                         label='CF (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot(cf_stats[:], zrange[:], 'k--', label='CF (stats)')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel('CF')
        plt.ylabel('height z (m)')
    plt.suptitle('CF from PDF sampling ' + '(t=' + str(time_field) + ')')
    save_name = 'cf_sampling_Lx_nc' + str(ncomp_max) + '_time' + str(time_field) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return




def plot_ql_error_allres_ncompmax(ql_mean_stats, ql_mean_field, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                        ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time, path):
    print('')
    print('plot error ql')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9*n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl+1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- <ql> (Lx=' + str(Lx_range[nl]) + ' dz=' + str(dz_range[ndz]) + ')')
            plt.plot(-ql_mean_stats[:], zrange[:], 'k--', label='- <ql> (stats)')
        plt.legend(loc=2)
        plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_allres_dz_nc' + str(ncomp_max)+ '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()



    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            for nc in range(ncomp_max):
                n_col = np.double(nc) / ncomp_max
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- <ql> (Lx=' + str(Lx_range[nl]) + ')')
            plt.plot(-ql_mean_stats[:], zrange[:], 'k--', label='- <ql> (stats)')
        plt.legend(loc=2)
        plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl]>1000:
                for nc in range(ncomp_max):
                    n_col = np.double(nc) / ncomp_max
                    plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=cm1(n_col),
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], '--', color=cm2(n_col), label='- <ql> (Lx=' + str(Lx_range[nl]) + ')')
                plt.plot(-ql_mean_stats[:], zrange[:], 'k--', label='- <ql> (stats)')
        plt.legend(loc=2)
        plt.xlim(xlimits_ql)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


def plot_ql_rel_error_allres_ncompmax(ql_mean_stats, ql_mean_field, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                                   ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time, path):
    print('')
    print('plot relative error ql')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.xlim([-1.1, 1.1])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql> / <ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_rel_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # # (2) one plot per dz (for all Lx)
    # plt.figure(figsize=(9 * n_lx, 12))
    # for ndz in range(n_dz):
    #     plt.subplot(1, n_dz, ndz + 1)
    #     for nl in range(n_lx):
    #         for nc in range(ncomp_max):
    #             col = cm1(np.double(nc) / ncomp_max)
    #             plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
    #                      label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
    #         plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
    #     plt.legend(loc=1)
    #     #plt.xlim(xlimits_ql)
    #     plt.title('dz=' + str(dz_range[ndz]) + 'm')
    #     plt.xlabel(r'$\epsilon(<ql>)$')
    #     plt.ylabel('height z (m)')
    # plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    # save_name = 'error_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.savefig(os.path.join(path_out, save_name))
    # # plt.show()
    # plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-ql_mean_field[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.xlim([-1.1,1.1])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql> / <ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_rel_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return




def plot_cf_error_allres_ncompmax(cf_stats, cf_field, error_cf, error_cf_rel, zrange,
                                  ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time, path):
    print('')
    print('plot error CF')
    # print('CF: ', cf_field)
    # print('e(CF): ', error_cf)
    # print('ql: ', ql_mean_field)
    # print('e(ql): ', error_ql)
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)
    # markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=2)
        plt.xlim(xlimits_cf)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<CF>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + ')')
    save_name = 'error_cf_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # (2) one plot per dz (for all Lx)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=2)
        # plt.xlim(-1.1,0.3)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<CF>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + ')')
    save_name = 'error_cf_allres_Lx__nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    # (2) one plot per dz (for all Lx>1000m)
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=2)
        plt.xlim(xlimits_cf)
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(CF)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error CF ' + '(t=' + str(time) + ')')
    save_name = 'error_cf_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


def plot_cf_rel_error_allres_ncompmax(cf_stats, cf_field, error_cf, error_cf_rel, zrange,
                                      ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, markers, time, path):
    print('')
    print('plot relative error CF')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    # print('n_lx', n_lx, 'n_dz', n_dz)

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        plt.xlim([-1.1, 1.1])
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(CF) / CF$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error CF' + '(t=' + str(time) + ')')
    save_name = 'error_cf_rel_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    # plt.show()
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    # # (2) one plot per dz (for all Lx)
    # plt.figure(figsize=(9 * n_lx, 12))
    # for ndz in range(n_dz):
    #     plt.subplot(1, n_dz, ndz + 1)
    #     for nl in range(n_lx):
    #         for nc in range(ncomp_max):
    #             col = cm1(np.double(nc) / ncomp_max)
    # #             plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
    # #                      label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
    # #         plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
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
    plt.figure(figsize=(9 * n_lx, 12))
    for ndz in range(n_dz):
        plt.subplot(1, n_dz, ndz + 1)
        for nl in range(n_lx):
            if Lx_range[nl] > 1000:
                for nc in range(ncomp_max):
                    col = cm1(np.double(nc) / ncomp_max)
                    plt.plot(error_cf_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                             label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
                plt.plot(-cf_field[nl, ndz, :], zrange[:], 'k--', label='- CF (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_cf)
        plt.xlim([-1.1, 1.1])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(CF) / CF$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error CF' + '(t=' + str(time) + ')')
    save_name = 'error_cf_rel_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


if __name__ == '__main__':
    main()