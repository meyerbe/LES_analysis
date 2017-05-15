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

# def plot_labels(data, labels_pdf, labels_tr, path):
def plot_labels(labels_pdf, labels_tr, type, path):
    print('')
    print('PLOTTING: '+type)
    [nx_pdf, ny_pdf, nz_pdf] = labels_pdf.shape
    [nx_tr, ny_tr, nz_tr] = labels_tr.shape
    print('')
    print('!!!')
    print(labels_pdf.shape, labels_tr.shape)
    # x_pdf = np.arange(0,nx_pdf)
    # y_pdf = np.arange(0,ny_pdf)
    # X, Y = np.meshgrid(x_, y_)
    #
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

    plt.figure(figsize=(12,3))

    plt.subplot(1,3,1)
    # plt.scatter(labels_pdf[:, :, 0], labels_pdf[:, :, 1], s=5, alpha=0.2)
    # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)
    # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)
    plt.imshow(labels_pdf[:,:,0])
    plt.colorbar()
    plt.title('PDF Labels')

    plt.subplot(1,3,2)
    plt.imshow(labels_tr[:,:,20])
    plt.colorbar()
    plt.title('Tracer Labels: '+type)
    # plt.scatter(labels_pdf[:, :, 0], labels_tr[:, :, 1], s=5, alpha=0.2)


    plt.savefig(os.path.join(path, 'Updrafts_figures', type+'.pdf'))
    # plt.show()
    plt.close()
    return