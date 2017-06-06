# Plotting Environment Data separately
# (1) Read in Couvreux Updraft Labels & PDF Labels
# (2) Select enviornment / updraft points


import os
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
    args = parser.parse_args()
    path = args.path
    path_out = os.path.join(path, 'PDF_cond_figures', 'scatter_plots')
    case_name = args.casename
    time_field = args.time_field

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    print('Files Dimensions: ' +  str(nx) +', '+ str(ny)+', ' + str(nz))

    krange, files = set_zrange(case_name)
    zrange = krange * dz
    nk = len(krange)
    print('Selected Files: ', files)
    print('Selected Levels: nk='+str(nk))
    print('krange: ', krange)
    print('zrange: ', zrange)
    print('')

    # (1) read in 3d Files: qt, ql, thetal data
    path_files = os.path.join(path, 'fields', time_field+'.nc')
    print(path_files)
    root = nc.Dataset(path_files, 'r')
    qt_ = root.groups['fields'].variables['qt'][:,:,:]
    ql_ = root.groups['fields'].variables['ql'][:, :, :]
    try:
        thetal_ = root.groups['fields'].variables['thetali'][:, :, :]
    except:
        # thetal = np.zeros(shape=qt_.shape)
        # for i in range(nx):
        #     for j in range(ny):
        #         for k in range(nz):
        #             aux = thetali_c(p_ref[iz + k_], T_[i, j, iz + k_], qt_[i, j, iz + k_], ql_[i, j, iz + k_], qi_, Lv)
        #             theta_l = np.append(theta_l, aux)
        print('!! thetal not in fields')
    root.close()

    # (2) read in labels Colleen (Couvreux)
    type = 'Couvreux'
    path_tracers = os.path.join(path, 'tracer_fields')
    labels_tracers = read_in_updrafts_colleen(type, time_field, path_tracers)

    # (3) read in labels PDF
    path_pdf_labels = os.path.join(path, 'Updrafts', 'Labeling_t'+str(time_field)+'.nc')
    print(path_pdf_labels)
    root = nc.Dataset(path_pdf_labels, 'r')
    labels_pdf = root.groups['fields'].variables['labels'][:,:,:]


    # (4) select the datapoints that are environemnt / updraft
    ini = time.time()

    # # aux1 = np.zeros(shape=(0))
    # # aux2 = np.zeros(shape=(0))
    # # aux3 = np.zeros(shape=(0))
    # # aux4 = np.zeros(shape=(0))
    # # aux5 = np.zeros(shape=(0))
    # # aux6 = np.zeros(shape=(0))
    # # qt_env_tr = np.zeros(shape=(0, nz), dtype=np.double)
    # # ql_env_tr = np.zeros(shape=(0, nz), dtype=np.double)
    # # th_env_tr = np.zeros(shape=(0, nz), dtype=np.double)
    # # qt_up_tr = np.zeros(shape=(0, nz), dtype=np.double)
    # # ql_up_tr = np.zeros(shape=(0, nz), dtype=np.double)
    # # th_up_tr = np.zeros(shape=(0, nz), dtype=np.double)
    qt = np.zeros(shape=(nx*ny,nk))
    ql = np.zeros(shape=(nx * ny, nk))
    thetal = np.zeros(shape=(nx * ny, nk))
    # n_up_tr = np.zeros(nk, dtype=np.int32)
    # n_env_tr = np.zeros(nk, dtype=np.int32)
    # n_up_pdf = np.zeros(nk, dtype=np.int32)
    # n_env_pdf = np.zeros(nk, dtype=np.int32)
    #
    for k in range(nk):
    # for k in krange:
    # for k in range(nz):
        iz = krange[k]
        print('k', k, iz, nz)

        n_up_tr = 0
        n_env_tr = 0
        n_up_pdf = 0
        n_env_pdf = 0

        qt_env_tr = np.zeros(shape=(0), dtype=np.double)
        ql_env_tr = np.zeros(shape=(0), dtype=np.double)
        th_env_tr = np.zeros(shape=(0), dtype=np.double)
        qt_up_tr = np.zeros(shape=(0), dtype=np.double)
        ql_up_tr = np.zeros(shape=(0), dtype=np.double)
        th_up_tr = np.zeros(shape=(0), dtype=np.double)
        qt_env_pdf = np.zeros(shape=(0), dtype=np.double)
        ql_env_pdf = np.zeros(shape=(0), dtype=np.double)
        th_env_pdf = np.zeros(shape=(0), dtype=np.double)
        qt_up_pdf = np.zeros(shape=(0), dtype=np.double)
        ql_up_pdf = np.zeros(shape=(0), dtype=np.double)
        th_up_pdf = np.zeros(shape=(0), dtype=np.double)

        for i in range(nx):
            ishift = i*ny
            for j in range(ny):
                ij = ishift + j
                qt[ij,k] = qt_[i,j,iz]
                ql[ij, k] = ql_[i, j,iz]
                thetal[ij,k] = thetal_[i,j,iz]
                if labels_tracers[i,j,iz] == 0:
                    # n_env_tr[k] += 1
                    n_env_tr += 1
                    qt_env_tr = np.append(qt_env_tr, qt_[i, j, iz])
                    ql_env_tr = np.append(ql_env_tr, ql_[i, j, iz])
                    th_env_tr = np.append(th_env_tr, thetal_[i, j, iz])
                else:
                    # n_up_tr[k] += 1
                    n_up_tr += 1
                    qt_up_tr = np.append(qt_up_tr, qt_[i, j, iz])
                    ql_up_tr = np.append(ql_up_tr, ql_[i, j, iz])
                    th_up_tr = np.append(th_up_tr, thetal_[i, j, iz])

                if labels_pdf[i,j,k] == 0:
                    # n_env_pdf[k] += 1
                    n_env_pdf += 1
                    qt_env_pdf = np.append(qt_env_pdf, qt_[i,j,iz])
                    ql_env_pdf = np.append(ql_env_pdf, ql_[i,j,iz])
                    th_env_pdf = np.append(th_env_pdf, thetal_[i,j,iz])
                else:
                    # n_up_pdf[k] += 1
                    n_up_pdf += 1
                    qt_up_pdf = np.append(qt_up_pdf, qt_[i, j, iz])
                    ql_up_pdf = np.append(ql_up_pdf, ql_[i, j, iz])
                    th_up_pdf = np.append(th_up_pdf, thetal_[i, j, iz])



        # type_ = 'Couvreux'
        # file_name = 'data_env_vs_all_' + type_  + '_z'+str(np.int(iz*dz)) + '_t' + str(time_field) + '.pdf'
        # plot_data_scatter(qt[:, k], qt_env_tr, qt_up_tr, thetal[:, k], th_env_tr, th_up_tr, ql[:, k], ql_env_tr, ql_up_tr,
        #               nx * ny, n_env_tr, n_up_tr, type_ , k * dz, time_field, path_out, file_name)
        # file_name = 'data_env_vs_all_' + type_  + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '_colored.pdf'
        # plot_data_scatter_colored(qt[:, k], qt_env_tr, qt_up_tr, thetal[:, k], th_env_tr, th_up_tr,
        #                           ql[:, k], ql_env_tr, ql_up_tr,
        #                             nx * ny, n_env_tr, n_up_tr, type_ , k * dz, time_field, path_out, file_name)
        type_ = 'PDF'
        file_name = 'data_env_vs_all_' + 'pdf' + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '.pdf'
        plot_data_scatter(qt[:, k], qt_env_pdf, qt_up_pdf, thetal[:, k], th_env_pdf, th_up_pdf, ql[:, k], ql_env_pdf, ql_up_pdf,
                      nx * ny, n_env_pdf, n_up_pdf, type_ , k * dz, time_field, path_out, file_name)
        file_name = 'data_env_vs_all_' + 'pdf' + '_z' + str(np.int(iz * dz)) + '_t' + str(time_field) + '_colored.pdf'
        plot_data_scatter_colored(qt[:, k], qt_env_pdf, qt_up_pdf, thetal[:, k], th_env_pdf, th_up_pdf,
                                  ql[:, k], ql_env_pdf, ql_up_pdf,
                      nx * ny, n_env_pdf, n_up_pdf, type_ , k * dz, time_field, path_out, file_name)

        end = time.time()
        print('k=' + str(k) + ', time: ' + str(end - ini))
        # for nk=2: time = end-ini = 12s
        # (a) no plotting
        # with scalar for n_env_tr etc.: for k=1: time = 3.3s, k=2: time=3.5s, k=3: time = 3.5s
        # with array for n_env_tr etc.: for k=1: time = 3.5s, k=2: time = 4s, k=3: time = 4.5s

    # print('shapes', qt_env_tr.shape, ql_env_tr.shape, th_env_tr.shape)



    return


#----------------------------------------------------------------------
def plot_data_scatter(qt_, qt_env_tr, qt_up_tr, thetal_, th_env_tr, th_up_tr, ql_, ql_env_tr, ql_up_tr,
                      n_tot, n_env_tr, n_up_tr, type_, z, t, path_, file_name):
    x_min = np.amin(thetal_)
    x_max = np.amax(thetal_)
    y_min = np.amin(qt_)
    y_max = np.amax(qt_)

    # (A) Figure without coloring
    plt.figure(figsize=(12,6))
    plt.subplot(1,3,1)
    plt.scatter(qt_, thetal_, s=5, alpha=0.2, edgecolors='none')
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.title('all data (n='+str(n_tot)+')')

    plt.subplot(1,3,2)
    plt.scatter(th_env_tr, qt_env_tr, s=5, alpha=0.2, edgecolors='none')
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.title('environment (n='+str(n_env_tr)+')')

    plt.subplot(1,3,3)
    plt.scatter(th_up_tr, qt_up_tr, s=5, alpha=0.2, edgecolors='none')
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.title('updrafts (n='+str(n_up_tr)+')')

    plt.suptitle('Updraft selection ' + type_ + ' (z='+str(z)+'m, t=' + str(t) +')')
    # # print('save:', os.path.join(path_, file_name))
    plt.savefig(os.path.join(path_, file_name))
    plt.close()

    return



def plot_data_scatter_colored(qt_, qt_env_tr, qt_up_tr, thetal_, th_env_tr, th_up_tr, ql_, ql_env_tr, ql_up_tr,
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
    plt.title('all data (n=' + str(n_tot) + ')')

    plt.subplot(1, 3, 2)
    plt.scatter(th_env_tr, qt_env_tr, c=ql_env_tr[:], s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    plt.colorbar(shrink=0.75)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.title('environment (n=' + str(n_env_tr) + ')')

    plt.subplot(1, 3, 3)
    plt.scatter(th_up_tr, qt_up_tr, c=ql_up_tr, s=5, alpha=0.5, edgecolors='none', vmin=ql_min, vmax=ql_max)#, cmap = cm)
    plt.colorbar(shrink=0.75)
    labeling(r'$\theta_l$', r'$q_t$', x_min, x_max, y_min, y_max)
    plt.title('updrafts (n=' + str(n_up_tr) + ')')

    plt.suptitle('Updraft selection ' + type_ + ' (z=' + str(z) + 'm, t=' + str(t) + ')')
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
        files = ['21600.nc']
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