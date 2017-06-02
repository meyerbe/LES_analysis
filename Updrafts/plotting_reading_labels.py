import os
import argparse
import json as simplejson
import numpy as np
import pylab as plt
import netCDF4 as nc

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 'small'
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 0.8
plt.rcParams['image.cmap'] = 'viridis'





def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("--path")
    parser.add_argument("--casename")
    parser.add_argument("--time")
    args = parser.parse_args()
    path = '/Volumes/Data/ClimatePhysics/LES/updrafts_colleen/'
    case_name = 'Bomex'
    time = 21600
    if args.path:
        path = args.path
    if args.casename:
        case_name = args.casename
    if args.time:
        time = np.int(args.time)

    global nx, ny, nz, dx, dy, dz
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']

    # Read in Fields
    path_fields = os.path.join(path, 'fields')
    qt_, ql_, w_, thetali_ = read_in_fields(time, path_fields)
    # Compute horizontal mean profiles
    # qt_mean = np.mean(np.mean(qt_, axis=0), axis=0)
    # ql_mean = np.mean(np.mean(ql_, axis=0), axis=0)
    # w_mean = np.mean(np.mean(w_, axis=0), axis=0)
    # thetali_mean = np.mean(np.mean(thetali_, axis=0), axis=0)
    # qt_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    # ql_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    # w_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    # thetali_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    # for i in range(nx):
    #     for j in range(ny):
    #         for k in range(nz):
    #             qt_anomaly[i,j,k] = qt_[i,j,k] - qt_mean[k]

    # Read in PDF Labels
    path_pdf = os.path.join(path, 'Updrafts')
    labels_pdf, zrange_pdf = read_in_pdf_labels(time, path_pdf)
    krange_pdf = zrange_pdf / dz
    nk = len(krange_pdf)

    # Read in Tracer Labels
    path_tracers = os.path.join(path, 'tracer_fields')
    # type_list = ['Couvreux']
    type_list = ['Couvreux', 'Coherent']
    for type in type_list:
        labels_tracers = read_in_updrafts_colleen(type, time, path_tracers)
        print('labels_tr', labels_tracers.shape)
        print('labels_pdf', labels_pdf.shape)
        print('zrange_pdf', zrange_pdf)
        print('krange_pdf', krange_pdf)


        dk = krange_pdf[1] - krange_pdf[0]
        for i in range(nk-1):
            aux = krange_pdf[i+1] - krange_pdf[i]
            if aux != dk:
                print('!!! zrange is not equidistant !!!', i, aux, dk)
        print('dk: ', dk)
        kmax = np.int(max(krange_pdf))
        kmax = 100
        kmin = 10
        print('kmin:', kmin)
        print('kmax: ', kmax, '(max k: ' + str(max(krange_pdf)) + ')')

        print('')

        # Plotting
        global cm1, cm2, cm3
        cm1 = plt.cm.get_cmap('viridis')
        cm2 = plt.cm.get_cmap('bone')
        cm3 = plt.cm.get_cmap('winter')
        krange_plot = krange_pdf[0:-1:5]
        # for k in range(len(krange_plot)):
        #     print('height: ', krange_plot[k]*dz)
        #     file_name = type + '_hor_t' + str(time) + '_z' + str(np.int(krange_plot[k]*dz)) + 'm.pdf'
        #     plot_labels_comparison_hor(qt_, ql_, w_, thetali_, labels_pdf, labels_tracers, type, time, k, krange_pdf, dz, path, file_name)
        # if dk == 1:
        #     for y0 in [10, 50]:
        #         file_name = type + '_vert_t' + str(time) + '_y' + str(y0) + 'm.pdf'
        #         plot_labels_comparison_vert(qt_, ql_, w_, thetali_, labels_pdf, labels_tracers, type, time, y0, krange_pdf, kmin, kmax, dz, path, file_name)

        plot_anomalies(qt_, ql_, w_, thetali_, labels_pdf, labels_tracers, type, time, krange_pdf, kmin, kmax, dz, path)

    return



def plot_anomalies(qt_, ql_, w_, thetali_, labels_pdf, labels_tracers, type, time, krange_pdf, kmin, kmax, dz, path):
    print('--- plot Anomalies ---')
    # Compute horizontal mean profiles
    qt_mean = np.mean(np.mean(qt_, axis=0), axis=0)
    ql_mean = np.mean(np.mean(ql_, axis=0), axis=0)
    w_mean = np.mean(np.mean(w_, axis=0), axis=0)
    thetali_mean = np.mean(np.mean(thetali_, axis=0), axis=0)
    qt_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    ql_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    w_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    thetali_anomaly = np.ndarray(shape=qt_.shape, dtype=np.double)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                qt_anomaly[i,j,k] = qt_[i,j,k] - qt_mean[k]
                ql_anomaly[i, j, k] = ql_[i, j, k] - ql_mean[k]
                w_anomaly[i, j, k] = w_[i, j, k] - w_mean[k]
                thetali_anomaly[i, j, k] = thetali_[i, j, k] - thetali_mean[k]

    krange_plot = krange_pdf[0:-1:5]
    print('krange plot: ', krange_plot)
    # for k in range(len(krange_plot)):
    #     file_name = type + '_hor_t' + str(time) + '_z' + str(np.int(krange_plot[k]*dz)) + 'm_anomaly.pdf'
    #     plot_labels_comparison_hor(qt_anomaly, ql_anomaly, w_anomaly, thetali_anomaly, labels_pdf, labels_tracers, type, time, k, krange_pdf, dz, path,
    #                            file_name)

    # for y0 in [10, 50]:
    #     file_name = type + '_vert_t' + str(time) + '_y' + str(y0) + 'm_anomaly.pdf'
    #     plot_labels_comparison_vert(qt_anomaly, ql_anomaly, w_anomaly, thetali_anomaly, labels_pdf, labels_tracers, type, time, y0, krange_pdf,
    #                                 kmin, kmax, dz, path, file_name)


    y0 = 10
    k = krange_pdf[0]
    file_name = 'w' + type + '_vert_t' + str(time) + '_y' + str(y0) + 'm_anomaly.pdf'
    plot_variable_vert(w_anomaly, 'w', labels_pdf, labels_tracers, type, time, y0, k, krange_pdf, kmin, kmax, dz, path, file_name)
    file_name = 'thetal' + type + '_vert_t' + str(time) + '_y' + str(y0) + 'm_anomaly.pdf'
    plot_variable_vert(thetali_anomaly, 'thetali', labels_pdf, labels_tracers, type, time, y0, k, krange_pdf, kmin, kmax, dz, path,
                       file_name)

    return






def read_in_pdf_labels(t, path_):
    print('--- read in PDF Labels ---')
    filename = 'Labeling_t'+str(t)+'.nc'
    # filename = 'Labeling_t' + '.nc'
    path = os.path.join(path_, filename)
    print(path)
    root = nc.Dataset(path, 'r')
    labels = root.groups['fields'].variables['labels'][:,:,:]
    zrange = root.groups['profiles'].variables['z'][:]
    root.close()
    print('')
    return labels, zrange


def read_in_updrafts_colleen(type, t, path_):
    print('')
    print('--- Updraft Colleen: read in ---')
    import pickle
    print('')


    path = path_
    files = os.listdir(path)
    print(path_)
    print('time: ', t)
    print('')

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

    print('')
    path = os.path.join(path_, root + 'time_'+ str(t) + '_Grid.pkl')
    print(path)
    print('')
    labels = pickle.load(open(path))
    # print(type(labels))
    # print(path + ': ', labels.shape())

    return labels


def plot_labels_comparison_hor(qt, ql, w, thetali, labels_pdf, labels_tr, type_, time, k, krange, dz, path, file_name):
    print('')
    print('PLOTTING hor: '+ type_)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
    n_pdf = labels_pdf.shape
    n_tr = labels_tr.shape
    print('')

    iz = np.int(krange[k])

    plt.figure(figsize=(24,6))

    plt.subplot(2, 6, 1)
    plt.imshow(w[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, :, k], [0.9], colors='w')
    plt.title(r'w')

    plt.subplot(2, 6, 2)
    plt.imshow(thetali[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, :, k], [0.9], colors='w')
    plt.title(r'$\theta_l$')

    plt.subplot(2, 6, 3)
    plt.imshow(qt[:,:,iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:,:,k], [0.9], colors='w')
    plt.title(r'$q_t$')

    plt.subplot(2, 6, 4)
    plt.imshow(ql[:,:,iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, :, k], [0.9], colors='w')
    plt.title(r'$q_l$')

    plt.subplot(2, 6, 5)
    plt.imshow(labels_pdf[:,:,k])
    plt.colorbar(shrink=0.5)
    plt.title('PDF Labels')

    plt.subplot(2, 6, 7)
    plt.imshow(w[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, :, k], [0.9], colors='w')
    plt.title(r'w')

    plt.subplot(2, 6, 8)
    plt.imshow(thetali[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, :, k], [0.9], colors='w')
    plt.title(r'$\theta_l$')

    plt.subplot(2, 6, 9)
    plt.imshow(qt[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, :, iz], [0.9], colors='w')
    plt.title(r'$q_t$')

    plt.subplot(2, 6, 10)
    plt.imshow(ql[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, :, iz], [0.9], colors='w')
    plt.title(r'$q_l$')

    plt.subplot(2, 6, 11)
    plt.imshow(labels_tr[:, :, iz])
    plt.colorbar(shrink=0.5)
    plt.title('Tracer Labels: '+type_, fontsize=10)

    if nx_pdf == nx_tr and ny_pdf == ny_tr:
        arr = np.zeros(shape=labels_pdf.shape, dtype=np.int16)
        n_env = 0
        n_updr = 0
        n_i = 0
        n_j = 0
        for i in range(nx_pdf):
            for j in range(ny_pdf):
                if labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 0:
                    arr[i,j,k] = 0
                    n_env += 1.
                    # print('label 0')
                elif labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 1:
                    arr[i, j, k] = 1
                    n_updr += 1.
                    # print('label 1')
                elif labels_pdf[i,j,k] == 1 and labels_tr[i,j,iz] == 0:
                    arr[i,j,k] = 2
                    # print('label 2')
                    n_i += 1.
                elif labels_pdf[i, j, k] == 0 and labels_tr[i, j, iz] == 1:
                    arr[i, j, k] = 3
                    # print('label 3')
                    n_j += 1.
        n_env /= (nx_pdf*ny_pdf)
        n_updr /= (nx_pdf * ny_pdf)
        n_i /= (nx_pdf * ny_pdf)
        n_j /= (nx_pdf * ny_pdf)
        # print('total (z='+str(iz*dz)+'): ', n_env, n_updr, n_i, n_j, n_env+n_updr+n_i+n_j)

        plt.subplot(2,6,6)
        plt.imshow(arr[:,:,k], clim=(0,3), interpolation="nearest")
        plt.title('env both: '+str(np.round(n_env*100,1))+'%, updr both: '+str(np.round(n_updr*100,1))+'%', fontsize=10)
        plt.colorbar(shrink=0.5)
        aux = np.zeros((8,8))
        aux[2:4] = 1
        aux[4:6] = 2
        aux[6:8] = 3
        plt.subplot(2, 6, 12)
        plt.imshow(aux, interpolation="nearest")
        plt.title('blue: good, yellow: PDF only, red: tracer only', fontsize=10)
        plt.colorbar(shrink=0.5)
    else:
        print('Dimensions do not agree')

    plt.suptitle('Updraft selection: ' + type_ + ' vs. PDF (z=' + str(iz*dz) + 'm)')
    plt.savefig(os.path.join(path, 'figures_Labels', file_name))
    # plt.show()
    plt.close()
    return


def plot_labels_comparison_vert(qt, ql, w, thetali, labels_pdf, labels_tr, type_, time, y0, krange, kmin, kmax, dz, path, file_name):
    print('')
    print('PLOTTING vertical: '+ type_)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
    n_pdf = labels_pdf.shape
    n_tr = labels_tr.shape
    print('')
    # print('pdf model: ', labels_pdf.shape)
    # print('tracer model: ', labels_tr.shape)
    # print('')

    # kmax = np.int(krange[-1])rxrange = np.arange(0,qt.shape[0],1)
    plt.figure(figsize=(30,8))

    plt.subplot(2, 6, 1)
    plt.imshow(w[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'w')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')

    plt.subplot(2, 6, 2)
    plt.imshow(thetali[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$\theta_l$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 3)
    plt.imshow(qt[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$q_t$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 4)
    plt.imshow(ql[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$q_l$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 5)
    plt.imshow(labels_pdf[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.title('PDF Labels')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 7)
    plt.imshow(w[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'w')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 8)
    plt.imshow(thetali[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$\theta_l$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 9)
    plt.imshow(qt[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$q_t$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 10)
    plt.imshow(ql[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(r'$q_l$')
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    plt.subplot(2, 6, 11)
    plt.imshow(labels_tr[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.title('Tracer Labels: '+type_, fontsize=10)
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'z (m)')

    if nx_pdf == nx_tr:
        arr = np.zeros(shape=labels_pdf[:,y0,kmin:kmax].shape, dtype=np.int16)
        n_env = 0
        n_updr = 0
        n_i = 0
        n_j = 0
        j = y0
        for i in range(nx_pdf):
            for k in range((kmax-kmin)):
                iz = np.int(krange[k])
                if labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 0:
                    arr[i,k] = 0
                    n_env += 1.
                    # print('label 0')
                elif labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 1:
                    arr[i, k] = 1
                    n_updr += 1.
                    # print('label 1')
                elif labels_pdf[i,j,k] == 1 and labels_tr[i,j,iz] == 0:
                    arr[i,k] = 2
                    # print('label 2')
                    n_i += 1.
                elif labels_pdf[i, j, k] == 0 and labels_tr[i, j, iz] == 1:
                    arr[i, k] = 3
                    # print('label 3')
                    n_j += 1.
        n_env /= (nx_pdf*ny_pdf)
        n_updr /= (nx_pdf * ny_pdf)
        n_i /= (nx_pdf * ny_pdf)
        n_j /= (nx_pdf * ny_pdf)
        print('total (z='+str(iz*dz)+'): ', n_env, n_updr, n_i, n_j, n_env+n_updr+n_i+n_j)

        plt.subplot(2,6,6)
        plt.imshow(arr[:,:].T, clim=(0,3), interpolation="nearest")
        plt.title('env both: '+str(np.round(n_env*100,1))+'%, updr both: '+str(np.round(n_updr*100,1))+'%', fontsize=10)
        set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')
        plt.colorbar(shrink=0.5)
        aux = np.zeros((8,4))
        aux[2:4] = 1
        aux[4:6] = 2
        aux[6:8] = 3
        plt.subplot(2, 6, 12)
        plt.imshow(aux.T, interpolation="nearest")
        plt.title('blue: good, yellow: PDF only, red: tracer only', fontsize=10)
        plt.colorbar(shrink=0.5)
    else:
        print('Dimensions do not agree')

    plt.suptitle('Updraft selection: ' + type_ + ' vs. PDF (y=' + str(y0 * dy) + 'm)')
    plt.savefig(os.path.join(path, 'figures_Labels', file_name))
    # plt.show()
    plt.close()
    return

#
# def plot_variable_hor(var, var_name, labels_pdf, labels_tr, type_, time, y0, k, krange, dz, path, file_name):
#     iz = np.int(krange[k])
#
#     plt.figure(12,6)
#     plt.subplot(1,2,1)
#     plt.imshow(var[:, :, iz])
#     plt.colorbar(shrink=0.5)
#     plt.contour(lables_pdf[:,:,k], [0.9], colors='w')
#     plt.title(var_name)
#
#     return


def plot_variable_vert(var, var_name, labels_pdf, labels_tr, type_, time, y0, k, krange, kmin, kmax, dz, path,
                      file_name):
    xrange = np.arange(0, var.shape[0], 1)

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.imshow(var[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(var_name)
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')

    plt.subplot(1,2,2)
    plt.imshow(var[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(var_name)
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')

    plt.suptitle('w: ' + type_ + ' vs. PDF (y=' + str(y0 * dy) + 'm)')
    plt.savefig(os.path.join(path, 'figures_Labels', 'im_' + file_name))
    plt.close()



    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.contourf(var[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_pdf[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(var_name)
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')

    plt.subplot(1, 2, 2)
    plt.contourf(var[:, y0, kmin:kmax].T)
    plt.colorbar(shrink=0.5)
    plt.contour(labels_tr[:, y0, kmin:kmax].T, [0.9], colors='w')
    plt.title(var_name)
    set_ticks(plt.gca(), xrange[:] * dx, krange[kmin:kmax] * dz, 'x (m)', 'k (m)')

    plt.suptitle('w: ' + type_ + ' vs. PDF (y=' + str(y0 * dy) + 'm)')
    plt.savefig(os.path.join(path, 'figures_Labels', 'cont_' + file_name))
    plt.close()

    return


def plot_labels(qt, ql, w, labels_, time, krange, dz, type_, path):
    for k in range(len(krange)):
        iz = krange[k]
        if type_ == 'PDF':
            k_ = k
        else:
            k_ = iz

        plt.figure(figsize=(12, 8))

        plt.subplot(2, 3, 1)
        plt.imshow(labels_[:, :, k_])
        plt.colorbar(shrink=0.6)
        plt.title(type_ +' Labels')

        plt.subplot(2, 3, 4)
        plt.imshow(qt[:, :, iz])
        plt.colorbar(shrink=0.6)
        # if np.amax(np.abs(labels_[:,:,k])) != 0.0:
        #     plt.contour(labels_[:, :, k],colors='w', levels=[0.9])
        plt.contour(labels_[:, :, k_] , [0.9],colors='w')
        plt.title(r'$q_t$ (cont: pdf label)')

        plt.subplot(2, 3, 5)
        plt.imshow(ql[:, :, iz])
        plt.colorbar(shrink=0.6)
        plt.contour(labels_[:, :, k_],[0.9], colors='w')
        plt.title(r'$q_l$')

        plt.subplot(2, 3, 6)
        plt.imshow(w[:, :, iz])
        plt.colorbar(shrink=0.6)
        plt.contour(labels_[:, :, k_], [0.9],colors='w')
        plt.title('w')

        plt.savefig(os.path.join(path, 'Updrafts_figures', 'Labels_' + type_ + '_t' + str(time) + '_z' + str(iz*dz) + 'm.pdf'))
        # # plt.show()
        plt.close()
    return


# ------------------------------------------------
def set_ticks(ax, var_x, var_y, xtitle, ytitle):
    # plt.gca().invert_yaxis()
    ax.invert_yaxis()
    global dz
    labels_x = ax.get_xticks()
    labels_y = ax.get_yticks()
    nx = len(var_x)
    ny = len(var_y)
    lab_x = [var_x[0]]
    lab_y = [var_y[0]]

    for i in range(1, labels_x.shape[0]):
        if labels_x[i] < nx:
            lab_x = np.append(lab_x, int(var_x[int(labels_x[i])]))
    lab_x = lab_x.astype(int)
    for i in range(1, labels_y.shape[0]):
        if labels_y[i] < ny:
            lab_y = np.append(lab_y, np.round(var_y[int(labels_y[i])], 3))
    ax.set_xticklabels(lab_x, fontsize=9)
    ax.set_yticklabels(lab_y, fontsize=9)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    return


def read_in_fields(t, path_):
    print('')
    print('--- read in 3D Fields ---')
    filename = str(t)+'.nc'
    path = os.path.join(path_, filename)
    root = nc.Dataset(path, 'r')
    qt = root.groups['fields'].variables['qt'][:,:,:]
    ql = root.groups['fields'].variables['ql'][:, :, :]
    w = root.groups['fields'].variables['w'][:, :, :]
    thetali = root.groups['fields'].variables['thetali'][:, :, :]
    root.close()
    return qt, ql, w, thetali


if __name__ == '__main__':
    main()
