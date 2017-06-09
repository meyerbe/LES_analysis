import numpy as np
import pylab as plt
import netCDF4 as nc
import sys, os
import json as simplejson
import argparse

sys.path.append("..")
from io_read_in_files import read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 2


def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    # parser.add_argument("--file_name")
    parser.add_argument("time")
    parser.add_argument("Lx", nargs='+', type=int)
    parser.add_argument("--dz", nargs='+', type=int)
    # parser.add_argument("time")
    parser.add_argument("--xa_ql")
    parser.add_argument("--xb_ql")
    parser.add_argument("--xa_ql_rel")
    parser.add_argument("--xb_ql_rel")
    parser.add_argument("--xa_cf")
    parser.add_argument("--xb_cf")
    parser.add_argument("--xa_cf_rel")
    parser.add_argument("--xb_cf_rel")

    args = parser.parse_args()
    path = args.path
    case_name = args.casename

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nz, dz
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    dt_stats = nml['stats_io']['frequency']
    t_ref = 0

    time_field = np.int(args.time)
    Lx_range = args.Lx
    if args.dz:
        delta_z_range = args.dz
    else:
        delta_z_range = [dz]
    time_name = '5to6h'
    print('Lx range: ', Lx_range)
    print('dz range: ', delta_z_range)


    path_profile = os.path.join(path, 'Stats.'+case_name+'.nc')

    type_ = 'Couvreux'
    for Lx_ in Lx_range:
        Lx = np.double(Lx_)
        for delta_z in delta_z_range:
            # file_name = args.file_name
            # file_name = 'PDF_cond_alltime_res_error_Lx' + str(np.round(Lx,1)) + 'Ly' + str(np.round(Lx,1)) \
            #             + '_dz' + str(np.int(delta_z)) + '_time' + str(time) + '.nc'

            file_name = 'PDF_cond_alltime_' + type_ + '_allcomp_Lx' + str(np.round(Lx, 1)) + 'Ly' + str(np.round(Lx, 1)) \
                    + '_dz' + str(np.int(delta_z)) + '_time' + str(time_field) + '.nc'
            path_in = os.path.join(path, 'PDF_cond', file_name)
            path_out = os.path.join(path, 'PDF_cond_figures', 'error_profiles')
            # path_in = os.path.join(path, 'PDF_cond_anomaly', file_name)
            # path_out = os.path.join(path, 'PDF_cond_anomaly_figures', 'error_profiles')
            print('')
            print('path: ', path)
            print('path_in: ', path_in)
            print('path_out: ', path_out)
            print('file name: ', file_name)
            print('')


            rootgrp = nc.Dataset(path_in, 'r')
            zrange = rootgrp.groups['profiles'].variables['zrange'][:]
            error_ql_env = rootgrp.groups['error'].variables['error_ql_env'][:, :]
            error_ql_rel_env = rootgrp.groups['error'].variables['rel_error_ql_env'][:, :]
            ql_mean_env = rootgrp.groups['profiles'].variables['ql_mean_env'][:]
            ql_mean_domain = rootgrp.groups['profiles'].variables['ql_mean_domain'][:]
            error_cf_env = rootgrp.groups['error'].variables['error_cf_env'][:, :]
            error_cf_rel_env = rootgrp.groups['error'].variables['rel_error_cf_env'][:, :]
            cf_env = rootgrp.groups['profiles'].variables['cf_env'][:]
            cf_domain = rootgrp.groups['profiles'].variables['cf_domain'][:]
            rootgrp.close()

            xa = -2.5e-6
            xb = 2.5e-6
            if args.xa_ql:
                xa = np.double(args.xa_ql)
            if args.xb_ql:
                xb = np.double(args.xb_ql)
            xlimits_ql = [xa, xb]

            xa = -5e-2
            xb = 1e-1
            if args.xa_cf:
                xa = np.double(args.xa_cf)
            if args.xb_cf:
                xb = np.double(args.xb_cf)
            xlimits_cf = [xa, xb]

            xa = -1.1
            xb = 1.1
            if args.xa_ql_rel:
                xa = np.double(args.xa_ql_rel)
            if args.xb_ql_rel:
                xb = np.doube(args.xb_ql_rel)
            xlimits_ql_rel = [xa, xb]

            xa = -1.1
            xb = 1.1
            if args.xa_cf_rel:
                xa = np.double(args.xa_cf_rel)
            if args.xb_cf_rel:
                xb = np.doube(args.xb_cf_rel)
            xlimits_cf_rel = [xa, xb]

            # # for ncomp_max in [3,5,8]:
            for ncomp_max in [5]:
            #     plot_error_ql_ncompmax(error_ql_env, error_ql_rel_env, ql_mean_env, zrange,
            #                            ncomp_max, xlimits_ql, xlimits_ql_rel, dz, Lx, delta_z, time_name, t_ref, type_, path_out)
            #     plot_error_cf_ncompmax(error_cf_env, error_cf_rel_env, cf_env, zrange,
            #                            ncomp_max, xlimits_cf, xlimits_cf_rel, dz, Lx, delta_z, time_name, t_ref, type_, path_out)

                plot_error_ql_domain(error_ql_env, error_ql_rel_env, ql_mean_env, ql_mean_domain, zrange,
                                       ncomp_max, xlimits_ql, xlimits_ql_rel, dz, Lx, delta_z, time_name, t_ref, type_, path_out)

                plot_error_cf_domain(error_cf_env, error_cf_rel_env, cf_env, cf_domain, zrange,
                                 ncomp_max, xlimits_cf, xlimits_cf_rel, dz, Lx, delta_z, time_name, t_ref, type_,
                                 path_out)


            # print('')
            # print(path_profile)
            # print('')
            # root = nc.Dataset(path_profile, 'r')
            # # cf_stats_ = root.groups['profiles'].variables['cloud_fraction'][:,:]
            # # ql_stats_ = root.groups['profiles'].variables['ql_mean'][:,:]
            # time_stats = root.groups['profiles'].variables['t'][:]
            # zrange_stats = root.groups['profiles'].variables['z'][:]
            # it_ref = np.int(time_field / dt_stats)
            # print('')
            # print('Times: ', time_field, time_stats[it_ref])
            # # print(time_stats.shape, it)
            # # krange = np.zeros(shape=zrange.shape)
            # # ql_stats = np.zeros(shape=zrange.shape)
            # # cf_stats = np.zeros(shape=zrange.shape)
            # # k_count = 0
            # # k = 0
            # # while(k_count < zrange.shape[0] and k < zrange_stats.shape[0]):
            # #     # for k in range(zrange_stats.shape[0]):
            # #     if zrange_stats[k] == zrange[k_count]:
            # #         # print('k', k, k_count)
            # #         # print(zrange_stats[k], zrange[k_count])
            # #         # print('')
            # #         krange[k_count] = k
            # #         ql_stats[k_count] = ql_stats_[it, k]
            # #         cf_stats[k_count] = cf_stats_[it, k]
            # #         k_count += 1
            # #     k += 1
            # # # print('krange: ', krange)
            # # print('')
            #
            #
            # ql_stats, cf_stats = plot_mean_profiles(ql_mean_env, ql_mean_domain, cf_env, cf_domain,
            #                    Lx, zrange, zrange_stats, time_field, it_ref, time_name, type_, path_profile, path_out)

    return




def plot_error_ql_ncompmax(error_ql_env, error_ql_rel_env, ql_mean_env, zrange,
                           ncomp_max, xlimits, xlimits_rel, dz, Lx, delta_z_, time, t_ref, type_, path_out):
    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_in = os.path.join(path, 'PDF_cond', file_name)
    # path_out = os.path.join(path, 'PDF_cond_figures', 'error_profiles')
    # print(path_in)
    # rootgrp = nc.Dataset(path_in, 'r')
    # zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    # error_ql_env = rootgrp.groups['error'].variables['error_ql_env'][:, :]
    # error_ql_rel_env = rootgrp.groups['error'].variables['rel_error_ql_env'][:, :]
    # ql_mean_env = rootgrp.groups['profiles'].variables['ql_mean_env'][:]
    # ql_mean_domain = rootgrp.groups['profiles'].variables['ql_mean_domain'][:]
    # rootgrp.close()

    delta_z = np.int(np.int(delta_z_))
    plt.figure(figsize=(18, 9))
    plt.subplot(1, 2, 1)
    plt.plot([0,0], [zrange[0],zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.plot(-ql_mean_env[:], zrange[:], 'k--', label='- mean ql env.')
    plt.xlim(xlimits)
    plt.legend(loc=2)
    plt.title('error <ql> environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm, t=' + str(time) + ')')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    plt.subplot(1, 2, 2)
    plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql_rel_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.xlim(xlimits_rel)
    plt.legend()
    plt.title('relative error <ql> environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_alltime_error_ql_nc' + str(ncomp_max) + '_dz' + str(delta_z) + '_Lx' + str(Lx) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    return



def plot_error_cf_ncompmax(error_cf_env, error_cf_rel_env, cf_env, zrange,
                           ncomp_max, xlimits, xlimits_rel, dz, Lx, delta_z_, time, t_ref, type_, path_out):
    # import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    # path_in = os.path.join(path, 'PDF_cond', file_name)
    # path_out = os.path.join(path, 'PDF_cond_figures', 'error_profiles')
    # print('')
    # print(path_in)
    # print('')
    # rootgrp = nc.Dataset(path_in, 'r')
    # zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    # error_cf = rootgrp.groups['error'].variables['error_cf_env'][:, :]
    # error_cf_rel = rootgrp.groups['error'].variables['rel_error_cf_env'][:, :]
    # cf_fields = rootgrp.groups['profiles'].variables['cf_domain'][:]
    # cf_env = rootgrp.groups['profiles'].variables['cf_env'][:]
    # rootgrp.close()

    delta_z = np.int(np.int(delta_z_))
    plt.figure(figsize=(18, 9))
    plt.subplot(1, 2, 1)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.plot(-cf_env[:], zrange[:], 'k--', label='- CF env.')
    plt.legend(loc=1)
    plt.xlim(xlimits)
    plt.title('error CF environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm, t=' + str(time) + ')')
    plt.xlabel(r'$\epsilon(CF)$')
    plt.ylabel('height z (m)')
    plt.subplot(1, 2, 2)
    plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf_rel_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.xlim(xlimits_rel)
    plt.legend()
    plt.title('relative error CF environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(CF)/CF$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_alltime_error_cf_nc' + str(ncomp_max) + '_dz' + str(delta_z) + '_Lx' + str(Lx) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    return


def plot_error_ql_domain(error_ql_env, error_ql_rel_env, ql_mean_env, ql_mean_domain, zrange,
                           ncomp_max, xlimits, xlimits_rel, dz, Lx, delta_z_, time, t_ref, type_, path_out):
    import netCDF4 as nc
    global nz

    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    print('..', error_ql_env.shape)
    print(type(ncomp_max))
    print(zrange)
    print(zrange/dz)
    print('')
    delta_z = np.int(np.int(delta_z_))
    krange = zrange / dz
    error_ql_rel_to_domain = np.zeros((len(krange), ncomp_max))
    for k in range(len(krange)):
        for nc in range(ncomp_max):
            if ql_mean_domain[k] > 0:
                error_ql_rel_to_domain[k,nc] = error_ql_env[k,nc] / ql_mean_domain[k]

    plt.figure(figsize=(18, 9))
    plt.subplot(1, 2, 1)
    plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.xlim(xlimits)
    plt.plot(-ql_mean_env[:], zrange[:], 'k--', label='- mean ql env.')
    plt.plot(-ql_mean_domain[:], zrange[:], '--', color='b', label='- mean ql domain')
    plt.legend(loc=2)
    plt.title('error <ql> environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm, t=' + str(time) + ')')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')

    plt.subplot(1, 2, 2)
    plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql_rel_env[:, nc], zrange[:], '-x', linewidth=1, color=col, label='error rel. to env (ncomp=' + str(nc + 1)+')')
        plt.plot(error_ql_rel_to_domain[:, nc], zrange[:], '--', color=col, label='error rel. to domain')
    plt.xlim(xlimits_rel)
    plt.legend()
    plt.title('relative error <ql> environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_domain_error_ql_domain_nc' + str(ncomp_max) + '_dz' + str(delta_z) + '_Lx' + str(Lx) + '_time' + str(
        time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))
    plt.close()

    return


def plot_error_cf_domain(error_cf_env, error_cf_rel_env, cf_env, cf_domain, zrange,
                           ncomp_max, xlimits, xlimits_rel, dz, Lx, delta_z_, time, t_ref, type_, path_out):
    # import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    delta_z = np.int(np.int(delta_z_))
    krange = zrange / dz
    error_cf_rel_to_domain = np.zeros((len(krange), ncomp_max))
    for k in range(len(krange)):
        for nc in range(ncomp_max):
            if cf_domain[k] > 0:
                error_cf_rel_to_domain[k, nc] = error_cf_env[k, nc] / cf_domain[k]

    delta_z = np.int(np.int(delta_z_))
    plt.figure(figsize=(18, 9))
    plt.subplot(1, 2, 1)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf_env[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.plot(-cf_env[:], zrange[:], 'k--', label='- CF env.')
    plt.plot(-cf_domain[:], zrange[:], '--', color='b', label='- CF domain.')
    plt.legend(loc=1)
    plt.xlim(xlimits)
    plt.title('error CF environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm, t=' + str(time) + ')')
    plt.xlabel(r'$\epsilon(CF)$')
    plt.ylabel('height z (m)')
    plt.subplot(1, 2, 2)
    plt.plot([0, 0], [zrange[0], zrange[-1]], 'k', linewidth=0.2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf_rel_env[:, n c], zrange[:], '-x', linewidth=1, color=col, label='error rel. to env (ncomp=' + str(nc + 1)+')')
        plt.plot(error_cf_rel_to_domain[:, nc], zrange[:], '--', color=col, label='error rel. to domain')
    plt.xlim(xlimits_rel)
    plt.legend()
    plt.title('relative error CF environment (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(CF)/CF$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_domain_error_cf_nc' + str(ncomp_max) + '_dz' + str(delta_z) + '_Lx' + str(
        Lx) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    return


def plot_mean_profiles(ql_mean_env, ql_mean_domain, cf_env, cf_domain,
                       Lx, zrange, zrange_stats, time_field, it_ref, time_, type_, path_profile, path_out):
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    root = nc.Dataset(path_profile, 'r')
    cf_stats_ = root.groups['profiles'].variables['cloud_fraction'][:, :]
    ql_stats_ = root.groups['profiles'].variables['ql_mean'][:, :]
    krange = np.zeros(shape=zrange.shape)
    ql_stats = np.zeros(shape=zrange.shape)
    cf_stats = np.zeros(shape=zrange.shape)
    k_count = 0
    k = 0
    while(k_count < zrange.shape[0] and k < zrange_stats.shape[0]):
        # for k in range(zrange_stats.shape[0]):
        if zrange_stats[k] == zrange[k_count]:
            # print('k', k, k_count)
            # print(zrange_stats[k], zrange[k_count])
            # print('')
            krange[k_count] = k
            ql_stats[k_count] = ql_stats_[it_ref, k]
            cf_stats[k_count] = cf_stats_[it_ref, k]
            k_count += 1
        k += 1
    print('')

    plt.figure(figsize=(9,9))
    plt.plot(ql_mean_env[:], zrange[:], label='<ql> env.')
    plt.plot(ql_mean_domain[:], zrange[:], label='<ql> domain.')
    plt.plot(ql_stats[:], zrange[:], label='<ql> stats')
    plt.legend()
    plt.title(' <ql> (Lx=' + str(Lx) + '(t=' + str(time_) + ')')
    plt.xlabel(r'$<ql>$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_ql_Lx' + str(Lx) + '_time' + str(time_) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    # plt.show()
    plt.close()


    plt.figure(figsize=(9, 9))
    plt.plot(cf_env[:], zrange[:], label='CF env.')
    plt.plot(cf_domain[:], zrange[:], label='CF domain.')
    plt.plot(cf_stats[:], zrange[:], label='CF stats')
    plt.legend()
    plt.title('CF (Lx=' + str(Lx) + '(t=' + str(time_) + ')')
    plt.xlabel(r'$CF$')
    plt.ylabel('height z (m)')
    save_name = type_ + '_cf_Lx' + str(Lx) + '_time' + str(time_) + '.pdf'
    plt.savefig(os.path.join(path_out, save_name))

    # plt.show()
    plt.close()

    return ql_stats, cf_stats



if __name__ == '__main__':
    main()
