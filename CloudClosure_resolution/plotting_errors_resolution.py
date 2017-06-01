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
plt.rcParams['lines.markersize'] = 9
plt.rcParams['lines.markeredgewidth'] = 0.0     # the line width around the marker symbol

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("time")
    parser.add_argument("--Lx")
    parser.add_argument("--dz")
    parser.add_argument("--ncomp_max")
    parser.add_argument("--xa_ql")
    parser.add_argument("--xb_ql")

    args = parser.parse_args()
    path = args.path


    time_field = args.time
    if args.Lx:
        Lx_range = np.round(np.double(args.Lx),1)
    else:
        Lx_range = [1000.0, 5000.0, 10000.0, 15000.0]
    if args.dz:
        dz_range = np.int(args.dz)
    else:
        dz_range = [20, 60, 100]
    if args.ncomp_max:
        ncomp_max = np.int(args.ncomp_max)
    else:
        ncomp_max = 3


    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']



    xa = -1.2e-5
    xb = 1e-6
    if args.xa_ql:
        xa = np.double(args.xa_ql)
    if args.xb_ql:
        xb = np.double(args.xb_ql)
    xlimits_ql = [xa, xb]
    xlimits_cf = [-5e-2, 2e-2]


    # Read in
    file_name = 'CC_alltime_res_error' + '_Lx10000.0Ly10000.0_dz20' + '_time' + time_field + '.nc'
    path_in = os.path.join(path, 'CloudClosure_res', file_name)
    print(path_in)

    rootgrp = nc.Dataset(path_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    nz = zrange.shape[0]
    error_ql = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_ql_rel = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_cf = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    error_cf_rel = np.zeros((len(Lx_range), len(dz_range), nz, ncomp_max))
    ql_mean_fields = np.zeros((len(Lx_range), len(dz_range), nz))
    print('..', error_ql.shape)
    for i in range(len(Lx_range)):
        Lx = Lx_range[i]
        for j in range(len(dz_range)):
            delta_z = dz_range[j]
            file_name = 'CC_alltime_res_error' + '_Lx' + str(Lx) + 'Ly' + str(Lx) + '_dz' + str(delta_z) + '_time' + time_field + '.nc'
            print(file_name)
            path_in = os.path.join(path, 'CloudClosure_res', file_name)
            rootgrp = nc.Dataset(path_in, 'r')
            error_ql[i, j, :, :] = rootgrp.groups['error'].variables['error_ql'][:, 0:ncomp_max]
            error_ql_rel[i, j, :, :] = rootgrp.groups['error'].variables['rel_error_ql'][:, 0:ncomp_max]
            ql_mean_fields[i,j, :] = rootgrp.groups['profiles'].variables['ql_mean_fields'][:]
            zrange_ = rootgrp.groups['profiles'].variables['zrange'][:]
            rootgrp.close()
            if zrange_.any() != zrange.any():
                print('differences in zrange')
                sys.exit()

    path_out = os.path.join(path, 'CloudClosure_res_figures')
    plot_ql_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                        ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time_field, path_out)
    plot_ql_rel_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                               ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time_field, path_out)

    return






def plot_ql_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                        ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time, path):
    print('')
    print('plot error ')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    print('n_lx', n_lx, 'n_dz', n_dz)
    markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9*n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl+1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
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
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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
        for nl in range(1,n_lx):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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




def plot_ql_rel_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                                   ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time, path):
    print('')
    print('plot relative error ')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    print('n_lx', n_lx, 'n_dz', n_dz)
    markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
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
    #         plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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
        for nl in range(1, n_lx):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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




def plot_cf_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                                      ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time, path):
    print('')
    print('plot error ')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    print('n_lx', n_lx, 'n_dz', n_dz)
    markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_cf[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
        plt.legend(loc=2)
        plt.xlim(xlimits_ql)
        plt.title('Lx=' + str(Lx_range[nl]) + 'm')
        plt.xlabel(r'$\epsilon(<ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_allres_dz_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
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
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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
        for nl in range(1, n_lx):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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


def plot_cf_rel_error_allres_ncompmax(ql_mean_fields, error_ql, error_ql_rel, error_cf, error_cf_rel, zrange,
                                      ncomp_max, xlimits_ql, xlimits_cf, dz, Lx_range, dz_range, time, path):
    print('')
    print('plot relative error ')
    # data = (Lx x dz x nz x ncomp)

    import netCDF4 as nc

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_out = os.path.join(path, 'error_profiles')

    n_lx = len(Lx_range)
    n_dz = len(dz_range)
    print('n_lx', n_lx, 'n_dz', n_dz)
    markers = ['o', 'v', '*', 'd']

    # (1) one plot per Lx (for all dz)
    plt.figure(figsize=(9 * n_lx, 10))
    for nl in range(n_lx):
        plt.subplot(1, n_lx, nl + 1)
        for ndz in range(n_dz):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[ndz], color=col,
                         label='ncomp=' + str(nc + 1) + ', dz=' + str(dz_range[ndz]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (dz=' + str(dz_range[ndz]) + ')')
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
    #         plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
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
        for nl in range(1, n_lx):
            for nc in range(ncomp_max):
                col = cm1(np.double(nc) / ncomp_max)
                plt.plot(error_ql_rel[nl, ndz, :, nc], zrange[:], '-', marker=markers[nl], color=col,
                         label='ncomp=' + str(nc + 1) + ', Lx=' + str(Lx_range[nl]))
            plt.plot(-ql_mean_fields[nl, ndz, :], zrange[:], 'k--', label='- mean ql (Lx=' + str(Lx_range[nl]) + ')')
        plt.legend(loc=1)
        # plt.xlim(xlimits_ql)
        plt.xlim([-1.1, 1.1])
        plt.title('dz=' + str(dz_range[ndz]) + 'm')
        plt.xlabel(r'$\epsilon(<ql> / <ql>)$')
        plt.ylabel('height z (m)')
    plt.suptitle('relative error <ql>' + '(t=' + str(time) + ')')
    save_name = 'error_ql_rel_allres_Lx_nc' + str(ncomp_max) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    # plt.show()
    plt.close()

    return


if __name__ == '__main__':
    main()