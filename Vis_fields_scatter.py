import os, sys
import argparse
import json as  simplejson
import numpy as np
import pylab as plt
import netCDF4 as nc
from matplotlib.colors import LogNorm

sys.path.append("./Thermo/")
from thermodynamic_functions import theta_li

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['figure.titlesize'] = 18
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
    nz = nml['grid']['nz']
    print('nz: '+str(nz))



    print('')
    path_fields = os.path.join(path, 'fields', str(time_field) + '.nc')
    print(path_fields)
    print('')

    if case_name != 'TRMM_LBA':
        zrange = np.arange(400, 2500, 400)
        print(zrange)
        krange = zrange / dz
        nk = len(krange)
    else:
        krange = np.arange(0,nz,50)
        # krange = np.asarray([0,20,40,85, 100, 120])
        krange = np.asarray([0, 5, 10, 20, 30, 40, 50, 60, 68, 75, 85, 95, 105, 120, 130, 140], dtype=np.int32)

        nk = len(krange)
        root = nc.Dataset(path_fields, 'r')
        zrange_ = root.groups['fields'].variables['z'][:]
        root.close()
        zrange = np.zeros(shape=krange.shape)
        for k in range(nk):
            zrange[k] = zrange_[krange[k]]


    print('krange: ' + str(krange))
    print('zrange: ' + str(zrange))
    print('')


    ''' read in reference pressure'''
    path_stats = os.path.join(path, 'Stats.'+case_name+'.nc')
    root = nc.Dataset(path_stats, 'r')
    p_ref = root.groups['reference'].variables['p0'][:]
    root.close()

    ''' read in Fields'''
    root = nc.Dataset(path_fields, 'r')
    variable_name = 'qt'
    qt_ = root.groups['fields'].variables[variable_name][:, :, :]
    variable_name = 'ql'
    ql_ = root.groups['fields'].variables[variable_name][:, :, :]
    variable_name = 'thetali'
    try:
        thetali_ = root.groups['fields'].variables[variable_name][:, :, :]
        thetali_flag = True
    except:
        print('!!! thetali not in fields variable !!!')
        thetali_flag = False
        T_ = root.groups['fields'].variables['temperature'][:, :, :]
        print('')
    root.close()

    qt = np.zeros((len(krange), nx*ny), dtype=np.double)
    ql = np.zeros((len(krange), nx * ny), dtype=np.double)
    thetali = np.zeros((len(krange), nx * ny), dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            for k in range(len(krange)):
                iz = np.int(krange[k])
                qt[k, i*ny+j] = qt_[i,j,iz]
                ql[k,i*ny+j] = ql_[i,j,iz]
                if thetali_flag:
                    thetali[k, i * ny + j] = thetali_[i, j, iz]
                else:
                    qi_ = 0
                    thetali[k, i*ny + j] = theta_li(p_ref[iz], T_[i,j,iz], qt_[i,j,iz], ql_[i,j,iz], qi_)
    ql_min = 0.0
    ql_max = np.amax(ql)

    # plt.figure(figsize=(25,5))
    # for k in range(nk):
    #     plt.subplot(1, nk, k + 1)
    #     plt.hist2d(thetali[k,:], qt[k,:], bins=100, normed=True)
    #     plt.colorbar()
    #     plt.title('z=' + str(zrange[k]) + 'm')
    #     plt.xlabel(r'$\theta_l$')
    #     plt.ylabel(r'$q_t$')
    #     # plt.xlim([298, 305])
    #     # plt.ylim([0.007,0.018])
    # plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    # plt.savefig(os.path.join(path, 'hist2d_plot.png'))
    # plt.close()
    #
    # plt.figure(figsize=(25, 5))
    # for k in range(nk):
    #     plt.subplot(1, nk, k + 1)
    #     plt.hist2d(thetali[k, :], qt[k, :], bins=100, norm = LogNorm(), normed = True)
    #     plt.colorbar()
    #     plt.title('z=' + str(zrange[k]) + 'm')
    #     plt.xlabel(r'$\theta_l$')
    #     plt.ylabel(r'$q_t$')
    #     # plt.xlim([298, 305])
    #     # plt.ylim([0.007,0.018])
    # plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    # plt.savefig(os.path.join(path, 'hist2d_plot_log.png'))
    # plt.close()

    if case_name != 'TRMM_LBA':
        scatter_plot(thetali, qt, ql, krange, zrange, time_field, path)
    else:
        # krange = np.asarray([0, 5, 10, 20, 30, 40, 50, 60, 68, 75, 85, 95, 105, 120, 130, 140], dtype=np.int32)
        # krange1 = np.asarray([0, 5, 10, 20, 30, 40, 50, 68, 75, 85, 95, 105, 120, 130], dtype=np.int32)
        krange1 = [0,1,2,3,4,5,6,8,9,10,11,12,13,14]
        scatter_plot_trmm(thetali, qt, ql, zrange, krange, krange1, time_field, path)
        krange = np.asarray([0, 5, 10, 20, 30, 40, 50, 60, 68, 75, 85, 95, 105, 120, 130], dtype=np.int32)
        # krange2 = np.asarray([10, 20, 40, 50, 60, 68, 105, 120], dtype=np.int32)
        krange2 = np.asarray([2, 3, 5, 9, 12, 13], dtype=np.int32)
        scatter_plot_trmm_small(thetali, qt, ql, zrange, krange, krange2, time_field, path)
    return





def scatter_plot(thetali, qt, ql, krange, zrange, time_field, path):
    nk = len(krange)
    if nk <= 6:
        plt.figure(figsize=(5*nk,5))
        for k in range(nk):
            plt.subplot(1,nk,k+1)
            plt.scatter(thetali[k,:], qt[k,:], s=3, alpha=0.1)
            plt.title('z='+str(np.int(zrange[k]))+'m')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            # plt.xlim([298, 305])
            # plt.ylim([0.007,0.018])
        plt.suptitle('LES Data (t='+str(np.round((time_field/3600),1))+'h)')
        # plt.suptitle('LES Data (t=' + 'h)')
        plt.savefig(os.path.join(path, 'scatter_plot.png'))
        plt.close()

        plt.figure(figsize=(25, 5))
        for k in range(nk):
            # ql_max = np.amax(ql[k,:])
            # print('ql: ', ql_max)
            plt.subplot(1, nk, k + 1)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
            # plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none')
            plt.colorbar(shrink=0.6)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            # plt.xlim([298, 305])
            # plt.ylim([0.007,0.018])
        plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
        plt.savefig(os.path.join(path, 'scatter_plot_ql.png'))
        plt.close()

    else:
        plt.figure(figsize=(5 * nk/2, 15))
        for k in range(nk):
            print(nk, '2', nk/2, k+1)
            plt.subplot(3, nk/2, k + 1)
            plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            # plt.xlim([298, 305])
            # plt.ylim([0.007,0.018])
        plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
        # plt.suptitle('LES Data (t=' + 'h)')
        plt.savefig(os.path.join(path, 'scatter_plot.png'))
        plt.close()

        plt.figure(figsize=(5*nk/2, 15))
        for k in range(nk):
            # ql_max = np.amax(ql[k,:])
            # print('ql: ', ql_max)
            plt.subplot(3, nk / 2, k + 1)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,
                        vmax=ql_max)  # , cmap = cm)
            # plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none')
            plt.colorbar(shrink=0.6)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            # plt.xlim([298, 305])
            # plt.ylim([0.007,0.018])
        plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
        plt.savefig(os.path.join(path, 'scatter_plot_ql.png'))
        plt.close()

    return


def scatter_plot_trmm_small(thetali, qt, ql, zrange, krange, krange_plot, time_field, path):
    nk = len(krange_plot)
    ql_min = 0.0
    ql_max = np.amax(ql)

    plt.figure(figsize=(5 * nk, 5))
    for k_ in range(nk):
        k = krange_plot[k_]
        plt.subplot(1, nk, k_ + 1)
        # k = krange[k_]
        plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
        plt.title('z=' + str(np.int(zrange[k])) + 'm')
        plt.xlabel(r'$\theta_l$')
        plt.ylabel(r'$q_t$')
    plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    # plt.suptitle('LES Data (t=' + 'h)')
    plt.savefig(os.path.join(path, 'scatter_plot_small.png'))
    plt.close()

    plt.figure(figsize=(5 * nk, 5))
    for k_ in range(nk):
        plt.subplot(1, nk, k_ + 1)
        k = krange_plot[k_]
        plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
        # plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none')
        plt.colorbar(shrink=0.6)
        plt.title('z=' + str(np.int(zrange[k])) + 'm')
        plt.xlabel(r'$\theta_l$')
        plt.ylabel(r'$q_t$')
    plt.suptitle('LES Data (t=' + str(np.round((time_field / 3600), 1)) + 'h)')
    plt.savefig(os.path.join(path, 'scatter_plot_ql_small.png'))
    plt.close()
    return


def scatter_plot_trmm(thetali, qt, ql, zrange, krange, krange_plot, time_field, path):
    nk = len(krange_plot)
    ql_min = 0.0
    ql_max = np.amax(ql)

    n_subc = 0
    n_shal = 0
    n_deep = 0
    k = 0
    while (zrange[krange_plot[k]] < 600):
        n_subc+=1
        k+=1
    while (zrange[krange_plot[k]] < 2500):
        n_shal += 1
        k += 1
    while (zrange[krange_plot[k]] < 8000):
        n_deep += 1
        k += 1
    nplot = max(n_subc, n_shal, n_deep)
    # print(nplot)

    plt.figure(figsize=(5 * nplot, 22))
    # fig1 = plt.figure(figsize=(5 * nplot, 20))
    n_subc = 1
    n_shal = 1
    n_deep = 1
    n_ft = 1
    print('krange', krange)
    print('krange_plot', krange_plot)
    print('zrange', zrange)

    for k_ in range(nk):
        k = krange_plot[k_]
        print(k_, krange_plot[k_], k, krange[k_], zrange[k_])
        if zrange[k] < 600:
            plt.subplot(4, nplot, k_ + 1)
            # fig1.subplot(4, nplot, k + 1)
            plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_subc += 1
        elif zrange[k] < 2500:
            plt.subplot(4, nplot, nplot + n_shal)
            plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_shal += 1
        elif zrange[k] < 8000:
            plt.subplot(4, nplot, 2*nplot + n_deep)
            plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_deep += 1
        else:
            plt.subplot(4, nplot, 3 * nplot + n_ft)
            plt.scatter(thetali[k, :], qt[k, :], s=3, alpha=0.1)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_ft += 1
    plt.suptitle('LES Data (t=' + str(np.round(((time_field-1e6)/ 3600), 1)) + 'h)')
    plt.savefig(os.path.join(path, 'scatter_plot.png'))
    plt.close()
    # fig1.close()


    plt.figure(figsize=(5 * nplot, 22))
    n_shal = 1
    n_deep = 1
    n_ft = 1
    for k_ in range(nk):
        k = krange_plot[k_]
        if zrange[k] < 600:
            plt.subplot(4, nplot, k_ + 1)
            # fig1.subplot(4, nplot, k + 1)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
        elif zrange[k] < 2500:
            plt.subplot(4, nplot, nplot + n_shal)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_shal += 1
        elif zrange[k] < 8000:
            plt.subplot(4, nplot, 2 * nplot + n_deep)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_deep += 1
        else:
            plt.subplot(4, nplot, 3 * nplot + n_ft)
            plt.scatter(thetali[k, :], qt[k, :], c=ql[k, :], s=5, alpha=0.2, edgecolors='none', vmin=ql_min,vmax=ql_max)  # , cmap = cm)
            plt.title('z=' + str(np.int(zrange[k])) + 'm')
            plt.xlabel(r'$\theta_l$')
            plt.ylabel(r'$q_t$')
            n_ft += 1
    plt.suptitle('LES Data (t=' + str(np.round(((time_field-1e6) / 3600), 1)) + 'h)')
    plt.savefig(os.path.join(path, 'scatter_plot_ql.png'))
    plt.close()
    return

if __name__=='__main__':
    main()