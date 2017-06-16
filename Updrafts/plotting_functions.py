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




def plot_parameters_pdf(labels_, means, covars, weights, krange, zrange, time_, dz, path):
    for k in range(len(krange)):
        iz = krange[k]
        k_ = k
        # if type_ == 'PDF':
        #     k_ = k
        # else:
        #     k_ = iz

        plt.figure(figsize=(12, 8))

        plt.subplot(2, 3, 1)
        plt.imshow(labels_[:, :, k_])
        plt.colorbar(shrink=0.6)
        plt.title('Labels')

        # plt.subplot(2, 3, 4)
        # plt.imshow(qt[:, :, iz])
        # plt.colorbar(shrink=0.6)
        # # if np.amax(np.abs(labels_[:,:,k])) != 0.0:
        # #     plt.contour(labels_[:, :, k],colors='w', levels=[0.9])
        # plt.contour(labels_[:, :, k_] , [0.9],colors='w')
        # plt.title(r'$q_t$ (cont: pdf label)')
        #
        # plt.subplot(2, 3, 5)
        # plt.imshow(ql[:, :, iz])
        # plt.colorbar(shrink=0.6)
        # plt.contour(labels_[:, :, k_],[0.9], colors='w')
        # plt.title(r'$q_l$')
        #
        # plt.subplot(2, 3, 6)
        # plt.imshow(w[:, :, iz])
        # plt.colorbar(shrink=0.6)
        # plt.contour(labels_[:, :, k_], [0.9],colors='w')
        # plt.title('w')

        plt.subplot(2,3,2)
        plt.plot(means[:, 0, 1], zrange, label='comp 1')
        plt.plot(means[:, 1, 1], zrange, label='comp 2')
        plt.legend()
        plt.title('<qt>')
        plt.subplot(2,3,3)
        plt.plot(means[:, 0, 0], zrange, label='comp 1')
        plt.plot(means[:, 1, 0], zrange, label='comp 2')
        plt.legend()
        plt.title(r'<$\theta_l$>')
        # plt.plot(weights[:, 0], krange, label='comp 1')
        # plt.plot(weights[:, 1], krange, label='comp 2')

        plt.savefig(os.path.join(path, 'Updrafts_figures', 'PDF_parameters' + '_t' + str(time_) + '_z' + str(iz*dz) + 'm.pdf'))
        # # plt.show()
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



# def plot_labels_test(qt, ql, w, labels_, time, krange, dz, type_, path):
#     print('Plotting Test')
#     for k in range(len(krange)):
#         iz = krange[k]
#         if type_ == 'PDF':
#             k_ = k
#         else:
#             k_ = iz
#         # print('data size', qt.shape, 'labels size', labels_.shape, k_, iz)
#
#         plt.figure(figsize=(12, 8))
#
#         plt.subplot(2, 3, 1)
#         plt.imshow(labels_[:, :, k_])
#         plt.colorbar(shrink=0.6)
#         plt.title(type_ +' Labels')
#
#         plt.subplot(2, 3, 4)
#         # plt.imshow(qt[:, :, iz])
#         plt.contour(ql[:, :, iz], [1e-5], colors = 'k')
#         # plt.colorbar(shrink=0.6)
#         # if np.amax(np.abs(labels_[:,:,k])) != 0.0:
#         #     plt.contour(labels_[:, :, k],colors='w', levels=[0.9])
#         # plt.contour(labels_[:, :, k_] , [0.9],colors='w')
#         plt.title(r'$q_l$ cont')
#
#         plt.subplot(2, 3, 5)
#         # plt.imshow(ql[:, :, iz])
#         # plt.colorbar(shrink=0.6)
#         plt.contour(labels_[:, :, k_], [0.9], colors='k')
#         plt.title(r'labels cont')
#
#         plt.subplot(2, 3, 6)
#         plt.imshow(ql[:, :, iz])
#         plt.colorbar(shrink=0.6)
#         plt.contour(labels_[:, :, k_], [0.9],colors='w')
#         plt.title('ql (cont=labels)')
#
#         plt.savefig(os.path.join(path, 'Updrafts_figures', 'test_Labels_' + type_ + '_t' + str(time) + '_z' + str(iz*dz) + 'm.pdf'))
#         # # plt.show()
#         plt.close()
#     return



def plot_labels_comparison(qt, ql, labels_pdf, labels_tr, type_, time, krange, dz, path):
    print('')
    print('PLOTTING: '+ type_)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
    n_pdf = labels_pdf.shape
    n_tr = labels_tr.shape
    print('')
    # print('pdf model: ', labels_pdf.shape)
    # print('tracer model: ', labels_tr.shape)
    # print('')

    for k in range(len(krange)):
        iz = krange[k]

        plt.figure(figsize=(12,6))

        plt.subplot(2, 4, 1)
        plt.imshow(qt[:,:,iz])
        plt.colorbar(shrink=0.5)
        plt.contour(labels_pdf[:,:,k], [0.9], colors='w')
        plt.title(r'$q_t$')

        plt.subplot(2, 4, 2)
        plt.imshow(ql[:,:,iz])
        plt.colorbar(shrink=0.5)
        plt.contour(labels_pdf[:, :, k], [0.9], colors='w')
        plt.title(r'$q_l$')

        plt.subplot(2, 4, 3)
        plt.imshow(labels_pdf[:,:,k])
        plt.colorbar(shrink=0.5)
        plt.title('PDF Labels')

        plt.subplot(2, 4, 5)
        plt.imshow(qt[:, :, iz])
        plt.colorbar(shrink=0.5)
        plt.contour(labels_tr[:, :, iz], [0.9], colors='w')
        plt.title(r'$q_t$')

        plt.subplot(2, 4, 6)
        plt.imshow(ql[:, :, iz])
        plt.colorbar(shrink=0.5)
        plt.contour(labels_tr[:, :, iz], [0.9], colors='w')
        plt.title(r'$q_l$')

        plt.subplot(2, 4, 7)
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
            print('total (z='+str(iz*dz)+'): ', n_env, n_updr, n_i, n_j, n_env+n_updr+n_i+n_j)

            plt.subplot(2,4,4)
            plt.imshow(arr[:,:,k], clim=(0,3), interpolation="nearest")
            plt.title('env both: '+str(np.round(n_env*100,1))+'%, updr both: '+str(np.round(n_updr*100,1))+'%', fontsize=10)
            plt.colorbar(shrink=0.5)
            aux = np.zeros((8,8))
            aux[2:4] = 1
            aux[4:6] = 2
            aux[6:8] = 3
            plt.subplot(2, 4, 8)
            plt.imshow(aux, interpolation="nearest")
            plt.title('blue: good, yellow: PDF only, red: tracer only', fontsize=10)
            plt.colorbar(shrink=0.5)
        else:
            print('Dimensions do not agree')

        plt.savefig(os.path.join(path, 'Updrafts_figures', type_ +'_t'+str(time)+'_z'+str(iz*dz)+'m.pdf'))
        # plt.show()
        plt.close()
    return