import os
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
# from sklearn.preprocessing import StandardScaler
import time
# import traceback


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 10
plt.rcParams['lines.linewidth'] = 0.8




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





def plot_labels_test(qt, ql, w, labels_, time, krange, dz, type_, path):
    print('Plotting Test')
    print('data size', qt.shape, 'labels size', labels_.shape)
    for k in range(len(krange)):
        iz = krange[k]
        if type_ == 'PDF':
            k_ = k
        else:
            k_ = iz
        # print('data size', qt.shape, 'labels size', labels_.shape, k_, iz)

        plt.figure(figsize=(12, 8))

        plt.subplot(2, 3, 1)
        plt.imshow(labels_[:, :, k_])
        plt.colorbar(shrink=0.6)
        plt.title(type_ +' Labels')

        plt.subplot(2, 3, 4)
        # plt.imshow(qt[:, :, iz])
        plt.contour(ql[:, :, iz], [1e-5], colors = 'k')
        # plt.colorbar(shrink=0.6)
        # if np.amax(np.abs(labels_[:,:,k])) != 0.0:
        #     plt.contour(labels_[:, :, k],colors='w', levels=[0.9])
        # plt.contour(labels_[:, :, k_] , [0.9],colors='w')
        plt.title(r'$q_l$ cont')

        plt.subplot(2, 3, 5)
        # plt.imshow(ql[:, :, iz])
        # plt.colorbar(shrink=0.6)
        plt.contour(labels_[:, :, k_], [0.9], colors='k')
        plt.title(r'labels cont')

        plt.subplot(2, 3, 6)
        plt.imshow(ql[:, :, iz])
        plt.colorbar(shrink=0.6)
        plt.contour(labels_[:, :, k_], [0.9],colors='w')
        plt.title('ql (cont=labels)')

        plt.savefig(os.path.join(path, 'Updrafts_figures', 'test_Labels_' + type_ + '_t' + str(time) + '_z' + str(iz*dz) + 'm.pdf'))
        # # plt.show()
        plt.close()
    return



def plot_labels_comparison(qt, ql, labels_pdf, labels_tr, type_, time, krange, dz, path):
    print('')
    print('PLOTTING: '+ type_)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
    n_pdf = labels_pdf.shape
    n_tr = labels_tr.shape
    print('')
    print('pdf model: ', labels_pdf.shape)
    print('tracer model: ', labels_tr.shape)
    print('')
    # x_pdf = np.arange(0,nx_pdf)
    # y_pdf = np.arange(0,ny_pdf)
    # X, Y = np.meshgrid(x_, y_)
    #
    #
    # # (1) make samples
    # n_sample = 300
    # x_ = np.linspace(np.amin(data[:, 0]), np.amax(data[:, 0]), n_sample)
    # y_ = np.linspace(np.amin(data[:, 1]), np.amax(data[:, 1]), n_sample)
    # XX_ = np.ndarray(shape=(n_sample ** nvar_, nvar_))
    # X_ = np.ndarray(shape=(n_sample, n_sample))
    # Y_ = np.ndarray(shape=(n_sample, n_sample))
    # delta_i = n_sample
    # for i in range(n_sample):
    #     for j in range(n_sample):
    #         shift = i * delta_i + j
    #         XX_[shift, 0] = x_[i]
    #         XX_[shift, 1] = y_[j]
    #         X_[i, j] = x_[j]
    #         Y_[i, j] = y_[i]
    # print('')

    for k in range(len(krange)):
        iz = krange[k]

        plt.figure(figsize=(12,8))

        plt.subplot(2, 3, 1)
        plt.imshow(qt[:,:,iz])
        plt.colorbar(shrink=0.5)
        plt.title(r'$q_t$')

        plt.subplot(2, 3, 4)
        plt.imshow(ql[:,:,iz])
        plt.colorbar(shrink=0.5)
        plt.title(r'$q_l$')

        plt.subplot(2, 3, 2)
        # plt.scatter(labels_pdf[:, :, 0], labels_pdf[:, :, 1], s=5, alpha=0.2)
        plt.imshow(labels_pdf[:,:,k])
        print('')
        # print('.........')
        # plt.contour(ql[:,:,iz], [0.0], colors = 'w')
        plt.colorbar(shrink=0.5)
        plt.title('PDF Labels')

        plt.subplot(2, 3, 3)
        plt.imshow(labels_tr[:,:,iz])
        plt.colorbar(shrink=0.5)
        plt.title('Tracer Labels: '+type_)
        # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)

        if nx_pdf == nx_tr and ny_pdf == ny_tr:
            arr = np.zeros(shape=labels_pdf.shape, dtype=np.int16)
            n_env = 0
            n_updr = 0
            for i in range(nx_pdf):
                for j in range(ny_pdf):
                    if labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 0:
                        arr[i,j,k] = 3
                        n_env += 1
                        # print('label 0')
                    elif labels_pdf[i, j, k] == labels_tr[i, j, iz] and labels_pdf[i,j,k] == 1:
                        arr[i, j, k] = 2
                        n_updr += 1
                        # print('label 1')
                    elif labels_pdf[i,j,k] == 1 and labels_tr[i,j,iz] == 0:
                        arr[i,j,k] = 1
                        # print('label 2')
                    elif labels_pdf[i, j, k] == 0 and labels_tr[i, j, iz] == 1:
                        arr[i, j, k] = 0
                        # print('label 3')
            n_env /= (nx_pdf*ny_pdf)
            n_updr /= (nx_pdf * ny_pdf)
            plt.subplot(2,3,5)
            plt.imshow(arr[:,:,k], interpolation="nearest")
            plt.title('0: '+str(n_env*100)+'%, 1: '+str(n_updr*100)+'%', fontsize=12)
            plt.colorbar(shrink=0.5)
            plt.subplot(2, 3, 6)
            aux = np.zeros((8,8))
            aux[2:4] = 1
            aux[4:6] = 2
            aux[6:8] = 3
            plt.imshow(aux, interpolation="nearest")
            plt.title('blue: good, yellow: PDF only, red: tracer only', fontsize=10)
            plt.colorbar(shrink=0.5)
        else:
            print('Dimensions do not agree')



        plt.savefig(os.path.join(path, 'Updrafts_figures', type_ +'_t'+str(time)+'_z'+str(iz*dz)+'m.pdf'))
        # plt.show()
        plt.close()
    return