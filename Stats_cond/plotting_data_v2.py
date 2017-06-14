# Plotting Environment Data separately
# (1) Read in Couvreux Updraft Labels & PDF Labels
# (2) Select enviornment / updraft points


import os, sys
import argparse
import json as simplejson
import pylab as plt
import numpy as np
import netCDF4 as nc
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import time


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 10
plt.rcParams['lines.linewidth'] = 0.8




def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("time_field")
    parser.add_argument("--type", nargs='+', type=str)
    args = parser.parse_args()
    path = args.path
    # path_out = os.path.join(path, 'PDF_cond_figures', 'scatter_plots')
    path_out = os.path.join(path, 'PDF_cond_t5h_figures', 'scatter_plots_log')
    case_name = args.casename
    time_field = args.time_field
    print('')
    print('path out: ' + path_out)

    if args.type:
        type_list = args.type
    else:
        type_list = ['Couvreux', 'PDF']


    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    print('')
    print('Files Dimensions: ' + str(nx) + ', ' + str(ny) + ', ' + str(nz))


    krange, files = set_zrange(case_name)
    if case_name != 'TRMM_LBA':
        zrange = krange * dz
    else:
        print('!!! Case Name TRMM_LBA: problem with zrange')
        sys.exit()
    nk = len(krange)
    print('time field: ' + str(time_field))
    # print('Selected Files: ', files)
    print('Selected Levels: nk=' + str(nk))
    print('krange: ', krange)
    print('zrange: ', zrange)
    print('')






    # (1) read in 3d Files: qt, ql, thetal data
    path_files = os.path.join(path, 'fields', time_field+'.nc')
    print(path_files)
    root = nc.Dataset(path_files, 'r')
    qt_ = root.groups['fields'].variables['qt'][:,:,:]
    ql_ = root.groups['fields'].variables['ql'][:, :, :]
    w_ = root.groups['fields'].variables['w'][:, :, :]
    try:
        buoyancy_ = root.groups['fields'].variables['buoyancy'][:,:,:]
        buoyancy_flag = True
    except:
        buoyancy_flag = False
        print('!! buoyancy not in fields')
    try:
        thetal_ = root.groups['fields'].variables['thetali'][:, :, :]
        thetal_flag = True
    except:
        # thetal = np.zeros(shape=qt_.shape)
        # for i in range(nx):
        #     for j in range(ny):
        #         for k in range(nz):
        #             aux = thetali_c(p_ref[iz + k_], T_[i, j, iz + k_], qt_[i, j, iz + k_], ql_[i, j, iz + k_], qi_, Lv)
        #             theta_l = np.append(theta_l, aux)
        thetal_flag = False
        print('!! thetal not in fields')
    root.close()

    # (2) read in Labels
    for type_ in type_list:
        print('')
        print('--- '+type_ + ' ---')
        print('')
        if type_=='PDF':
            path_labels = os.path.join(path, 'Updrafts', 'Labeling_t' + str(time_field) + '.nc')
            print(path_labels)
            root = nc.Dataset(path_labels, 'r')
            labels = root.groups['fields'].variables['labels'][:, :, :]
        else:
            path_labels = os.path.join(path, 'tracer_fields')
            labels = read_in_updrafts_colleen(type_, time_field, path_labels)


        # (3) select the datapoints that are environemnt / updraft
        ini = time.time()
        qt = np.zeros(shape=(nx*ny,nk))
        ql = np.zeros(shape=(nx * ny, nk))
        if thetal_flag:
            thetal = np.zeros(shape=(nx * ny, nk))
        w = np.zeros(shape=(nx * ny, nk))
        if buoyancy_flag:
            buoy = np.zeros(shape=(nx * ny, nk))

        n_up = np.zeros(nk, dtype=np.int32)
        n_env = np.zeros(nk, dtype=np.int32)

        for k in range(nk):
            iz = krange[k]
            print('k', k, iz, nz)
            qt_env = np.zeros(shape=(0), dtype=np.double)
            qt_up = np.zeros(shape=(0), dtype=np.double)
            ql_env = np.zeros(shape=(0), dtype=np.double)
            ql_up = np.zeros(shape=(0), dtype=np.double)
            if thetal_flag:
                th_env = np.zeros(shape=(0), dtype=np.double)
                th_up = np.zeros(shape=(0), dtype=np.double)
            w_env = np.zeros(shape=(0), dtype=np.double)
            w_up = np.zeros(shape=(0), dtype=np.double)
            if buoyancy_flag:
                buoy_env = np.zeros(shape=(0), dtype=np.double)
                buoy_up = np.zeros(shape=(0), dtype=np.double)

            for i in range(nx):
                ishift = i*ny
                for j in range(ny):
                    ij = ishift + j
                    qt[ij,k] = qt_[i,j,iz]
                    ql[ij, k] = ql_[i, j,iz]
                    w[ij,k] = w_[i,j,iz]
                    if thetal_flag:
                        thetal[ij,k] = thetal_[i,j,iz]
                    if buoyancy_flag:
                        buoy[ij,k] = buoyancy_[i,j,iz]
                    if labels[i,j,iz] == 0:
                        n_env[k] += 1
                        qt_env = np.append(qt_env, qt_[i, j, iz])
                        ql_env = np.append(ql_env, ql_[i, j, iz])
                        w_env = np.append(w_env, w_[i,j,iz])
                        if thetal_flag:
                            th_env = np.append(th_env, thetal_[i, j, iz])
                        if buoyancy_flag:
                            buoy_env = np.append(buoy_env, buoyancy_[i, j, iz])
                    else:
                        n_up[k] += 1
                        qt_up = np.append(qt_up, qt_[i, j, iz])
                        ql_up = np.append(ql_up, ql_[i, j, iz])
                        w_up = np.append(w_up, w_[i,j,iz])
                        if thetal_flag:
                            th_up = np.append(th_up, thetal_[i, j, iz])
                        if buoyancy_flag:
                            buoy_up = np.append(buoy_up, buoyancy_[i, j, iz])

            ql_mean_up = np.average(ql_up[:])
            ql_mean_env = np.average(ql_env[:])
            qt_mean_up = np.average(qt_up[:])
            qt_mean_env = np.average(qt_env[:])
            w_mean_up = np.average(w_up[:])
            w_mean_env = np.average(w_env[:])
            b_mean_up = np.average(buoy_up[:])
            b_mean_env = np.average(buoy_env[:])


            # Plot log(qt) vs. qt
            file_name = 'data_log_' + type_  + '_z'+str(np.int(iz*dz)) + '_t' + str(time_field) + '.png'
            plot_data_scatter(qt[:,k], thetal[:,k], ' ', ' ',
                              qt_env, th_env, ' ', ' ',
                              np.log(qt_env), th_env, 'log(qt)', 'th_l',
                              nx * ny, n_env[k], n_up[k], type_, iz * dz, time_field, path_out, file_name)

            # Plot normal data
            # file_name = 'data_env_vs_all_' + type_  + '_z'+str(np.int(iz*dz)) + '_t' + str(time_field) + '.png'
            # plot_data_scatter(qt[:, k], thetal[:, k], ' ', ' ',
            #                   qt_env, th_env, ' ', ' ',
            #                   qt_up, th_up, ' ', ' ',
            #                nx * ny, n_env[k], n_up[k], type_ , iz * dz, time_field, path_out, file_name)
            # # plot data scattered, colored by ql
            # file_name = 'data_env_vs_all_' + type_  + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '_colored_ql.png'
            # plot_data_scatter_colored(qt[:, k], qt_env, qt_up, thetal[:, k], th_env, th_up,
            #                           ql[:, k], ql_env, ql_up, 'ql', ql_mean_up, ql_mean_env,
            #                             nx * ny, n_env[k], n_up[k], type_ , iz * dz, time_field, path_out, file_name)
            # # plot data scattered, colored by w
            # file_name = 'data_env_vs_all_' + type_ + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '_colored_w.png'
            # plot_data_scatter_colored(qt[:, k], qt_env, qt_up, thetal[:, k], th_env, th_up,
            #                       w[:, k], w_env, w_up, 'w', w_mean_up, w_mean_env,
            #                       nx * ny, n_env[k], n_up[k], type_, iz * dz, time_field, path_out, file_name)
            # if buoyancy_flag:
            #     # plot data scattered, colored by buoyancy
            #     file_name = 'data_env_vs_all_' + type_ + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '_colored_b.png'
            #     plot_data_scatter_colored(qt[:, k], qt_env, qt_up, thetal[:, k], th_env, th_up,
            #                           buoy[:, k], buoy_env, buoy_up, 'buoyancy', b_mean_up, b_mean_env,
            #                           nx * ny, n_env[k], n_up[k], type_, iz * dz, time_field, path_out, file_name)


    return


#----------------------------------------------------------------------
def plot_data_scatter(qt_, thetal_, name1a, name2a,
                      qt_env_tr, th_env_tr, name1b, name2b,
                      qt_up_tr, th_up_tr, name1c, name2c,
                      n_tot, n_env_tr, n_up_tr, type_, z, t, path_, file_name):
    x_min = np.amin(thetal_)
    x_max = np.amax(thetal_)
    y_min = np.amin(qt_)
    y_max = np.amax(qt_)

    # (A) Figure without coloring
    plt.figure(figsize=(12,6))
    plt.subplot(1,3,1)
    plt.scatter(thetal_, qt_, s=5, alpha=0.2, edgecolors='none')
    if name1a == ' ' or name2a == ' ':
        labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    else:
        labeling(name1a, name2a, x_min, x_max, y_min, y_max)
    plt.title('all data (n='+str(n_tot)+')')

    plt.subplot(1,3,2)
    plt.scatter(th_env_tr, qt_env_tr, s=5, alpha=0.2, edgecolors='none')
    if name1b == ' ' or name2b == ' ':
        labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    else:
        labeling(name1b, name2b, x_min, x_max, y_min, y_max)
    plt.title('environment (n='+str(n_env_tr)+')')

    plt.subplot(1,3,3)
    plt.scatter(th_up_tr, qt_up_tr, s=5, alpha=0.2, edgecolors='none')
    if name1c == ' ' or name2c == ' ':
        labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    else:
        plt.xlabel(name2c)
        plt.ylabel(name1c)
        # pass
        # labeling(name1c, name2c, x_min, x_max, y_min, y_max)
    plt.title('updrafts (n='+str(n_up_tr)+')')


    plt.suptitle('Updraft selection ' + type_ + ' (z='+str(z)+'m, t=' + str(t) +')')
    # # print('save:', os.path.join(path_, file_name))
    plt.savefig(os.path.join(path_, file_name))
    plt.close()

    return



def plot_data_scatter_colored(qt_, qt_env_tr, qt_up_tr, thetal_, th_env_tr, th_up_tr,
                              ql_, ql_env_tr, ql_up_tr, color_name, color_mean_up,  color_mean_env,
                          n_tot, n_env_tr, n_up_tr, type_, z, t, path_, file_name):
    # (B) Figure with ql-coloring
    ql_min = np.amin(ql_)
    ql_max = np.amax(ql_)
    x_min = np.amin(thetal_)
    x_max = np.amax(thetal_)
    y_min = np.amin(qt_)
    y_max = np.amax(qt_)

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 3, 1)
    plt.scatter(thetal_, qt_, c=ql_[:], s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.colorbar(shrink=0.75)
    plt.title('all data (n=' + str(n_tot) + ')', fontsize=10)

    plt.subplot(1, 3, 2)
    plt.scatter(th_env_tr, qt_env_tr, c=ql_env_tr[:], s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    plt.colorbar(shrink=0.75)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    if color_name == 'ql':
        mean = str(np.round(color_mean_env * 1e5, 2)) + r'$\cdot 10^{-5}$'
    else:
        mean = str(np.round(color_mean_env, 2))
    plt.title('environment, mean=' + mean +' (n=' + str(n_env_tr) + ')', fontsize=10)

    plt.subplot(1, 3, 3)
    plt.scatter(th_up_tr, qt_up_tr, c=ql_up_tr, s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    plt.colorbar(shrink=0.75)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    if color_name == 'ql':
        mean = str(np.round(color_mean_up*1e5,2))+r'$\cdot 10^{-5}$'
    else:
        mean = str(np.round(color_mean_up,2))
    plt.title('updrafts, mean='+ mean + ' (n=' + str(n_up_tr) + ')', fontsize=10)

    plt.suptitle('Updraft selection ' + type_ + ', color: '+  color_name + ' (z=' + str(z) + 'm, t=' + str(t) + ')')
    # # print('save:', os.path.join(path_, file_name))
    plt.savefig(os.path.join(path_, file_name))
    # plt.show()
    # print('saved')
    plt.close()
    return


def plot_hist2d(qt_, qt_env_tr, qt_up_tr, thetal_, th_env_tr, th_up_tr, ql_, ql_env_tr, ql_up_tr,
                                  n_tot, n_env_tr, n_up_tr, type_, z, t, path_, file_name):
    # (B) Figure with ql-coloring
    # ql_min = np.amin(ql_)
    # ql_max = np.amax(ql_)
    # x_min = np.amin(thetal_)
    # x_max = np.amax(thetal_)
    # y_min = np.amin(qt_)
    # y_max = np.amax(qt_)

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 3, 1)
    plt.hist2d(thetal_, qt_)#, c=ql_[:], s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)  # , cmap = cm)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.colorbar(shrink=0.75)
    plt.title('all data (n=' + str(n_tot) + ')')

    # plt.subplot(1, 3, 2)
    # plt.scatter(th_env_tr, qt_env_tr, c=ql_env_tr[:], s=5, alpha=0.5, edgecolors='none', vmin=ql_min,
    #             vmax=ql_max)  # , cmap = cm)
    # plt.colorbar(shrink=0.75)
    # labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    # plt.title('environment (n=' + str(n_env_tr) + ')')
    #
    # plt.subplot(1, 3, 3)
    # plt.scatter(th_up_tr, qt_up_tr, c=ql_up_tr, s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)  # , cmap = cm)
    # plt.colorbar(shrink=0.75)
    # labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    # plt.title('updrafts (n=' + str(n_up_tr) + ')')
    #
    # plt.suptitle('Updraft selection ' + type_ + ' (z=' + str(z) + 'm, t=' + str(t) + ')')
    # print('save:', os.path.join(path_, file_name))
    # plt.savefig(os.path.join(path_, file_name))
    # plt.show()
    # print('saved')
    plt.close()
    return

def labeling(var_name1, var_name2, x_min, x_max, y_min, y_max):
    plt.xlabel(var_name1)
    plt.ylabel(var_name2)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    return

# ----------------------------------
def read_in_updrafts_colleen(type, t, path_):
    print('')
    print('- Updraft Colleen: read in -')
    import pickle


    path = path_
    files = os.listdir(path)
    print(path_)
    print('time: ', t)

    if type == 'Cloud':
        root = 'Bomex_Cloud_'
    elif type == 'Coherent':
        root = 'Bomex_Coherent_'
    elif type == 'Couvreux':
        root = 'Bomex_Couvreux_'
    elif type == 'Core':
        root = 'Bomex_Core_'
    print(root)

    # print('')
    # path = os.path.join(path_, root + 'updraft.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())
    #
    # print('')
    # path = os.path.join(path_, root + 'environment.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())

    path = os.path.join(path_, root + 'time_'+ str(t) + '_Grid.pkl')
    print(path)
    print('')
    # print('')
    f = open(path, 'r')
    labels = pickle.load(f)
    # labels = pickle.load(open(path))
    return labels
    # return


def set_zrange(case_name):
    if case_name[0:8] == 'ZGILS_S6':
        files = ['1382400.nc']
        krange = np.asarray([25, 25, 40, 50, 60, 65], dtype=np.int32)
    elif case_name[0:9] == 'ZGILS_S12':
        files = ['86400.nc']
        krange = np.asarray([35,40,45])
    elif case_name == 'DYCOMS_RF01':
        # DYCOMS RF01 large
        krange = np.asarray([135, 140, 145, 150, 155, 160, 166, 180], dtype=np.int32)
        # files = ['3600.nc']
        # DYCOMS RF01
        # krange = np.asarray([140, 150, 160, 166, 180], dtype=np.int32)
        files = ['10800.nc', '14400.nc']
    elif case_name == 'DYCOMS_RF02':
        krange = np.asarray([120, 140, 150, 160, 170, 200], dtype=np.int32)
        files = ['14400.nc']
    elif case_name == 'Bomex':
        ## Bomex 170314_weno7 (dz=40)
        # krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50], dtype=np.int32)
        # files = ['21600.nc']
        ## Bomex (dz=20)
        # krange = np.asarray([40, 50, 60], dtype=np.int32)
        krange = np.asarray([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 125], dtype=np.int32)
        # files = ['18000.nc', '19800.nc', '21600.nc']
        # files = ['21600.nc']
        files = ['18000.nc']
        # Bomex test
        # files = ['21600.nc']
        # krange = np.asarray([10, 17, 20, 25, 50])
        # krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50])
        # krange = np.asarray([20, 50])
        # krange = np.asarray([18,30,38])
    elif case_name == 'TRMM_LBA':
        # TRMM
        # files = ['1012600.nc', '1014400.nc', '1016200.nc']
        files = ['1014400.nc']
        krange = np.asarray([15, 20, 25, 30, 35, 40, 45], dtype=np.int32)

    return krange, files


if __name__=='__main__':
    main()