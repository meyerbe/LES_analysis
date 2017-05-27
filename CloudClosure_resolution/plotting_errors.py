import pylab as plt
import numpy as np
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
    parser.add_argument("Lx")
    parser.add_argument("dz")
    parser.add_argument("time")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    Lx = args.Lx
    delta_z = args.dz
    time_field = args.time

    fullpath = os.path.join(path, 'CloudClosure_res')
    file_name_err = 'CC_res_error_Lx'+ str(Lx)+'.0Ly'+str(Lx)+'.0_dz'+str(delta_z)+'_time'+str(time_field)+'.nc'
    path_ref = os.path.join(path, 'Stats.' + case_name + '.nc')

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    dt_stats = nml['stats_io']['frequency']
    time_ref = read_in_netcdf('t', 'timeseries', path_ref)
    t_ref = np.int(0)
    while ((time_ref[t_ref] - np.int(time_field) ) < dt_stats and t_ref < (time_ref.shape[0]-1) ):
        t_ref = np.int(t_ref + 1)
    print('t_ref: ', t_ref, time_ref.shape[0])

    xlimits_ql = [-9e-6,1e-6]
    xlimits_cf = [-5e-2, 1e-2]
    # plot_error_ql(case_name, path, file_name_err, xlimits, dz, Lx, delta_z, time_field, t_ref)
    for ncomp_max in [3,5]:
        plot_error_ql_ncompmax(case_name, path, file_name_err, ncomp_max, xlimits_ql, dz, Lx, delta_z, time_field, t_ref)
        plot_error_cf_ncompmax(fullpath, file_name_err, ncomp_max, xlimits_cf, dz, Lx, delta_z, time_field, t_ref)


    return




def plot_error_ql_ncompmax(case_name, path, file_name, ncomp_max, xlimits, dz, Lx, delta_z_, time, t_ref):
    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_in = os.path.join(path, 'CloudClosure_res', file_name)
    path_out = os.path.join(path, 'CloudClosure_res_figures')
    path_ref = os.path.join(path, 'Stats.' + case_name + '.nc')
    print('')
    print('plot error: ', path_in)
    print('')
    rootgrp = nc.Dataset(path_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    delta_z = np.int(np.int(delta_z_) + dz)
    error_ql = rootgrp.groups['error'].variables['error_ql'][:,:]
    error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:,:]
    ql_mean_fields = rootgrp.groups['profiles'].variables['ql_mean_fields'][:]
    rootgrp.close()

    print('')
    [nz, ncomp] = error_ql.shape
    print('zrange: '+str(zrange), ' nz '+ str(nz))
    print('ncomp ' + str(ncomp_max))
    plt.figure(figsize=(18, 9))
    plt.subplot(1, 2, 1)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.legend()
    plt.title('error <ql> (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm, t=' + str(time) + ')')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    plt.subplot(1, 2, 2)
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.plot(-ql_mean_fields[:], zrange[:], 'k--', label='- mean ql')
    plt.xlim(xlimits)
    plt.legend()
    plt.title('error <ql> (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    save_name = 'error_ql_nc' + str(ncomp_max) + '_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))
    print(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_ql_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
    plt.legend()
    plt.title('relative error <ql> (Lx='+str(Lx)+', dz='+ str(delta_z) +'m)')
    plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
    plt.ylabel('height z (m)')
    save_name = 'error_ql_rel_nc'+str(ncomp_max)+ '_Lx'+str(Lx) + '_dz' + str(delta_z) + '_time' + str(time) +'.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))

    print('')
    return

def plot_error_ql(case_name, path, file_name, xlimits, dz, Lx, delta_z_, time, t_ref):
    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_in = os.path.join(path, 'CloudClosure_res', file_name)
    path_out = os.path.join(path, 'CloudClosure_res_figures')
    path_ref = os.path.join(path, 'Stats.' + case_name + '.nc')
    print('')
    print('plot error: ', path_in)
    print('')
    rootgrp = nc.Dataset(path_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    delta_z = np.int(np.int(delta_z_) + dz)
    error_ql = rootgrp.groups['error'].variables['error_ql'][:,:]
    error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:,:]
    ql_mean_fields = rootgrp.groups['profiles'].variables['ql_mean_fields'][:]
    # ql_mean_pdf = rootgrp.groups['profiles'].variables['ql_mean_pdf'][:]
    # ql_mean_domain = rootgrp.groups['profiles'].variables['ql_mean_domain'][:]
    rootgrp.close()



    print('')
    [nz, ncomp] = error_ql.shape
    print('zrange: '+str(zrange), ' nz '+ str(nz))
    print('ncomp ' + str(ncomp))

    plt.figure(figsize=(18,9))
    plt.subplot(1,2,1)
    for nc in range(ncomp):
        col = cm1(np.double(nc)/10)
        plt.plot(error_ql[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
    plt.legend()
    plt.title('error <ql> (Lx='+str(Lx)+', dz='+ str(delta_z) +'m, t=' + str(time)+')')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    plt.subplot(1, 2, 2)
    for nc in range(ncomp):
        col = cm1(np.double(nc) / 10)
        plt.plot(error_ql[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc + 1))
    plt.plot(-ql_mean_fields[:], zrange[:], 'k--', label='- mean ql')
    plt.legend()
    plt.title('error <ql> (Lx=' + str(Lx) + ', dz=' + str(delta_z) + 'm)')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    save_name = 'error_ql_Lx'+str(Lx)+'_dz'+ str(delta_z) +'_time'+str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))
    print(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp):
        col = cm1(np.double(nc) / 10)
        plt.plot(error_ql_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
    plt.legend()
    plt.title('relative error <ql> (Lx='+str(Lx)+', dz='+ str(delta_z) +'m)')
    plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_ql_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)  + '.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))

    print('')
    return

def plot_error_cf_ncompmax(fullpath, file_name, ncomp_max, xlimits, dz, Lx, delta_z_, time, t_ref):
    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_in = os.path.join(fullpath, file_name)
    path_out = fullpath + '_figures'
    print('')
    print('plot error: ', path_in)
    rootgrp = nc.Dataset(path_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    delta_z = np.int(np.int(delta_z_) + dz)
    error_ql = rootgrp.groups['error'].variables['error_ql'][:]
    error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:]
    error_cf = rootgrp.groups['error'].variables['error_cf'][:]
    error_cf_rel = rootgrp.groups['error'].variables['rel_error_cf'][:]
    rootgrp.close()

    [nz, ncomp] = error_ql.shape
    print('zrange: '+str(zrange), ' nz '+ str(nz))
    print('ncomp ' + str(ncomp))

    plt.figure(figsize=(10,8))
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
    plt.legend()
    plt.xlim(xlimits)
    plt.title('error CF (Lx='+str(Lx)+', dz='+ str(delta_z) +'m, t=' + str(time)+')')
    plt.xlabel(r'$\epsilon(CF)$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_cf_nc'+str(ncomp_max) + '_Lx'+str(Lx)+'_dz'+ str(delta_z) +'_time'+str(time) + '.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))
    # print(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp_max):
        col = cm1(np.double(nc) / ncomp_max)
        plt.plot(error_cf_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
    plt.legend()
    plt.xlim([-1,1e-2])
    plt.title('relative error CF (Lx='+str(Lx)+', dz='+ str(delta_z) +'m)')
    plt.xlabel(r'$\epsilon(CF)/CF$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_cf_rel_nc' + str(ncomp_max) + '_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)  + '.pdf'
    plt.savefig(os.path.join(path_out, 'error', save_name))

    print('')
    return

# def plot_error(fullpath, file_name, xlimits, dz, Lx, delta_z, time):
#     import netCDF4 as nc
#     cm1 = plt.cm.get_cmap('viridis')
#     cm2 = plt.cm.get_cmap('bone')
#     cm3 = plt.cm.get_cmap('winter')
#
#     path_in = os.path.join(fullpath, file_name)
#     path_out = fullpath + '_figures'
#     print('')
#     print('plot error: ', path_in)
#     print('')
#     rootgrp = nc.Dataset(path_in, 'r')
#     zrange = rootgrp.groups['profiles'].variables['zrange'][:]
#     error_ql = rootgrp.groups['error'].variables['error_ql'][:]
#     error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:]
#     error_cf = rootgrp.groups['error'].variables['error_cf'][:]
#     error_cf_rel = rootgrp.groups['error'].variables['rel_error_cf'][:]
#     rootgrp.close()
#
#     print('')
#     [nz, ncomp] = error_ql.shape
#     print('zrange: '+str(zrange), ' nz '+ str(nz))
#     print('ncomp ' + str(ncomp))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc)/10)
#         plt.plot(error_ql[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
#     plt.legend()
#     plt.title('error <ql> (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(<ql>)$')
#     plt.ylabel('height z (m)')
#     save_name = 'error_ql_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#     print(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_ql_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
#     plt.legend()
#     plt.title('relative error <ql> (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_ql_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_cf[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
#     plt.legend()
#     plt.title('error CF (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(CF)$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_cf_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#     print(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_cf_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
#     plt.legend()
#     plt.title('relative error CF (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(CF)/CF$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_cf_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#
#     print('')
#     return

# def plot_error(fullpath, file_name, xlimits, dz, Lx, delta_z, time):
#     import netCDF4 as nc
#     cm1 = plt.cm.get_cmap('viridis')
#     cm2 = plt.cm.get_cmap('bone')
#     cm3 = plt.cm.get_cmap('winter')
#
#     path_in = os.path.join(fullpath, file_name)
#     path_out = fullpath + '_figures'
#     print('')
#     print('plot error: ', path_in)
#     print('')
#     rootgrp = nc.Dataset(path_in, 'r')
#     zrange = rootgrp.groups['profiles'].variables['zrange'][:]
#     error_ql = rootgrp.groups['error'].variables['error_ql'][:]
#     error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:]
#     error_cf = rootgrp.groups['error'].variables['error_cf'][:]
#     error_cf_rel = rootgrp.groups['error'].variables['rel_error_cf'][:]
#     rootgrp.close()
#
#     print('')
#     [nz, ncomp] = error_ql.shape
#     print('zrange: '+str(zrange), ' nz '+ str(nz))
#     print('ncomp ' + str(ncomp))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc)/10)
#         plt.plot(error_ql[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
#     plt.legend()
#     plt.title('error <ql> (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(<ql>)$')
#     plt.ylabel('height z (m)')
#     save_name = 'error_ql_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#     print(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_ql_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
#     plt.legend()
#     plt.title('relative error <ql> (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_ql_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_cf[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
#     plt.legend()
#     plt.title('error CF (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(CF)$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_cf_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#     print(os.path.join(path_out, save_name))
#
#     plt.figure()
#     for nc in range(ncomp):
#         col = cm1(np.double(nc) / 10)
#         plt.plot(error_cf_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
#     plt.legend()
#     plt.title('relative error CF (Lx='+str(Lx)+', dz='+str(delta_z)+'m)')
#     plt.xlabel(r'$\epsilon(CF)/CF$')
#     plt.ylabel('height z (m)')
#     # plt.show()
#     save_name = 'error_cf_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
#     plt.savefig(os.path.join(path_out, save_name))
#
#     print('')
#     return


if __name__ == '__main__':
    main()
