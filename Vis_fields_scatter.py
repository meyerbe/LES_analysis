import os, sys
import argparse
import json as  simplejson
import numpy as np
import pylab as plt
import netCDF4 as nc
from matplotlib.colors import LogNorm



label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['figure.titlesize'] = 21
plt.rcParams['lines.linewidth'] = 2


def main():
    sys.path.append('/Volumes/Data/ClimatePhysics/LES/output/')

    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("time")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    time_field = np.int(args.time)

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']


    print('')
    zrange = np.arange(400, 2500, 400)
    print(zrange)
    krange = zrange / dz
    nk = len(krange)
    print(krange)
    print('')

    path_fields = os.path.join(path, 'fields', str(time_field)+'.nc')
    print(path_fields)

    rootgrp = nc.Dataset(path_fields, 'r')
    variable_name = 'qt'
    qt_ = rootgrp.groups['fields'].variables[variable_name][:, :, :]
    variable_name = 'ql'
    ql_ = rootgrp.groups['fields'].variables[variable_name][:, :, :]
    variable_name = 'thetali'
    thetali_ = rootgrp.groups['fields'].variables[variable_name][:, :, :]

    qt = np.zeros((len(krange), nx*ny), dtype=np.double)
    ql = np.zeros((len(krange), nx * ny), dtype=np.double)
    thetali = np.zeros((len(krange), nx * ny), dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            for k in range(len(krange)):
                iz = np.int(krange[k])
                qt[k, i*ny+j] = qt_[i,j,iz]
                thetali[k, i * ny + j] = thetali_[i, j, iz]
                ql[k,i*ny+j] = ql_[i,j,iz]
    ql_min = 0.0
    ql_max = np.amax(ql)
    print('ql: ', ql_max)

    # plt.figure(figsize=(25,5))
    # for k in range(nk):
    #     plt.subplot(1,nk,k+1)
    #     plt.scatter(thetali[k,:], qt[k,:], s=3, alpha=0.1)
    #     plt.title('z='+str(zrange[k])+'m')
    #     plt.xlabel(r'$\theta_l$')
    #     plt.ylabel(r'$q_t$')
    #     # plt.xlim([298, 305])
    #     # plt.ylim([0.007,0.018])
    # plt.suptitle('LES Data (t='+str(np.round((time_field/3600),1))+'h)')
    # # plt.suptitle('LES Data (t=' + 'h)')
    # plt.savefig(os.path.join(path, 'scatter_plot.png'))
    # plt.close()
    #
    # plt.figure(figsize=(25, 5))
    # for k in range(nk):
    #     # ql_max = np.amax(ql[k,:])
    #     # print('ql: ', ql_max)
    #     plt.subplot(1, nk, k + 1)
    #     plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    #     # plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none')
    #     plt.colorbar(shrink=0.6)
    #     plt.title('z=' + str(zrange[k]) + 'm')
    #     plt.xlabel(r'$\theta_l$')
    #     plt.ylabel(r'$q_t$')
    #     # plt.xlim([298, 305])
    #     # plt.ylim([0.007,0.018])
    # plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    # plt.savefig(os.path.join(path, 'scatter_plot_ql.png'))
    # plt.close()

    plt.figure(figsize=(25, 5))
    for k in range(nk):
        plt.subplot(1, nk, k + 1)
        plt.hist2d(thetali[k,:], qt[k,:], bins=100, normed=True)
        plt.colorbar()
        plt.title('z=' + str(zrange[k]) + 'm')
        plt.xlabel(r'$\theta_l$')
        plt.ylabel(r'$q_t$')
        # plt.xlim([298, 305])
        # plt.ylim([0.007,0.018])
    plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    plt.savefig(os.path.join(path, 'hist2d_plot.png'))
    plt.close()

    plt.figure(figsize=(25, 5))
    for k in range(nk):
        plt.subplot(1, nk, k + 1)
        plt.hist2d(thetali[k, :], qt[k, :], bins=100, norm = LogNorm(), normed = True)
        plt.colorbar()
        plt.title('z=' + str(zrange[k]) + 'm')
        plt.xlabel(r'$\theta_l$')
        plt.ylabel(r'$q_t$')
        # plt.xlim([298, 305])
        # plt.ylim([0.007,0.018])
    plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    plt.savefig(os.path.join(path, 'hist2d_plot_log.png'))
    plt.close()



    return


if __name__=='__main__':
    main()