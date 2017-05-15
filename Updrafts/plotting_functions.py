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
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 3




def plot_labels_pdf(qt, ql, w, labels_pdf, time, krange, dz, path):
    for k in range(len(krange)):
        iz = krange[k]

        plt.figure(figsize=(12, 3))

        plt.subplot(2, 3, 1)
        plt.imshow(labels_pdf[:, :, k])
        plt.colorbar()
        plt.title('PDF Labels')

        plt.subplot(2, 3, 4)
        plt.imshow(qt[:, :, iz])
        plt.contour(labels_pdf[:, :, k])
        plt.colorbar()
        plt.title(r'$q_t$')

        plt.subplot(2, 3, 5)
        plt.imshow(ql[:, :, iz])
        plt.contour(labels_pdf[:, :, k])
        plt.colorbar()
        plt.title(r'$q_l$')

        plt.subplot(2, 3, 6)
        plt.imshow(w[:, :, iz])
        plt.contour(labels_pdf[:, :, k])
        plt.colorbar()
        plt.title('w')

        plt.savefig(os.path.join(path, 'Updrafts_figures', 'Labels_PDF' + '_t' + str(time) + '_z' + str(iz) + 'm.pdf'))
        # plt.show()
        plt.close()
    return



# def plot_labels(data, labels_pdf, labels_tr, path):
def plot_labels(qt, ql, labels_pdf, labels_tr, type, time, krange, dz, path):
    print('')
    print('PLOTTING: '+ type)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
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

        plt.figure(figsize=(12,3))

        plt.subplot(1, 4, 1)
        plt.imshow(qt[:,:,iz])
        plt.colorbar()
        plt.title(r'$q_t$')

        plt.subplot(1, 4, 2)
        plt.imshow(ql[:,:,iz])
        plt.colorbar()
        plt.title(r'$q_l$')

        plt.subplot(1, 4, 3)
        # plt.scatter(labels_pdf[:, :, 0], labels_pdf[:, :, 1], s=5, alpha=0.2)
        # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)
        # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)
        plt.imshow(labels_pdf[:,:,k])
        plt.colorbar()
        plt.title('PDF Labels')

        plt.subplot(1, 4, 4)
        plt.imshow(labels_tr[:,:,iz])
        plt.colorbar()
        plt.title('Tracer Labels: '+type)
        # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)


        plt.savefig(os.path.join(path, 'Updrafts_figures', type+'_t'+str(time)+'_z'+str(iz)+'m.pdf'))
        # plt.show()
        plt.close()
    return