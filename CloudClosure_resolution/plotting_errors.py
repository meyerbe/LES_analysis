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
plt.rcParams['axes.labelsize'] = 10
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
    time = args.time

    fullpath = os.path.join(path, 'CloudClosure_res')
    file_name_err = 'CC_res_error_Lx'+ str(Lx)+'.0Ly'+str(Lx)+'.0_dz'+str(delta_z)+'_time'+str(time)+'.nc'

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']

    plot_error(fullpath, file_name_err, dz, Lx, delta_z, time)


    return


def plot_error(fullpath, file_name, dz, Lx, delta_z, time):
    import netCDF4 as nc
    cm1 = plt.cm.get_cmap('viridis')
    cm2 = plt.cm.get_cmap('bone')
    cm3 = plt.cm.get_cmap('winter')

    path_in = os.path.join(fullpath, file_name)
    path_out = fullpath + '_figures'
    print('')
    print('plot error: ', path_in)
    print('')
    rootgrp = nc.Dataset(path_in, 'r')
    zrange = rootgrp.groups['profiles'].variables['zrange'][:]
    error_ql = rootgrp.groups['error'].variables['error_ql'][:]
    error_ql_rel = rootgrp.groups['error'].variables['rel_error_ql'][:]
    error_cf = rootgrp.groups['error'].variables['error_cf'][:]
    error_cf_rel = rootgrp.groups['error'].variables['rel_error_cf'][:]
    rootgrp.close()

    print('')
    [nz, ncomp] = error_ql.shape
    print('zrange: '+str(zrange), ' nz '+ str(nz))
    print('ncomp ' + str(ncomp))


    plt.figure()
    for nc in range(ncomp):
        col = cm1(np.double(nc)/10)
        plt.plot(error_ql[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
    plt.legend()
    plt.title('error <ql>')
    plt.xlabel(r'$\epsilon(<ql>)$')
    plt.ylabel('height z (m)')
    save_name = 'error_ql_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
    plt.savefig(os.path.join(path_out, save_name))
    print(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp):
        col = cm1(np.double(nc) / 10)
        plt.plot(error_ql_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
    plt.legend()
    plt.title('relative error <ql>')
    plt.xlabel(r'$\epsilon(<ql>)/<ql>$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_ql_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
    plt.savefig(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp):
        col = cm1(np.double(nc) / 10)
        plt.plot(error_cf[:,nc], zrange[:], '-x', color=col, label='ncomp='+str(nc+1))
    plt.legend()
    plt.title('error CF')
    plt.xlabel(r'$\epsilon(CF)$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_cf_Lx'+str(Lx)+'_dz'+str(delta_z)+'_time'+str(time)
    plt.savefig(os.path.join(path_out, save_name))
    print(os.path.join(path_out, save_name))

    plt.figure()
    for nc in range(ncomp):
        col = cm1(np.double(nc) / 10)
        plt.plot(error_cf_rel[:, nc], zrange[:], '-x', color=col, label='ncomp=' + str(nc+1))
    plt.legend()
    plt.title('relative error CF')
    plt.xlabel(r'$\epsilon(CF)/CF$')
    plt.ylabel('height z (m)')
    # plt.show()
    save_name = 'error_cf_rel_Lx' + str(Lx) + '_dz' + str(delta_z) + '_time' + str(time)
    plt.savefig(os.path.join(path_out, save_name))

    print('')
    return


if __name__ == '__main__':
    main()
