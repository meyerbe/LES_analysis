
import pylab as plt
import numpy as np
import matplotlib.mlab as mlab
from matplotlib.colors import LogNorm
import matplotlib.cm as cm

label_size = 12
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 12
# plt.rcParams['title.fontsize'] = 15


def gaussian(x,x0,sigma):
  return 1/np.sqrt(2*np.pi*sigma)*np.exp(-np.power((x - x0)/sigma, 2.)/2.)

def main():
    x = 1
    x0 = 0
    sigma = 1
    a = gaussian(x, x0, sigma)

    x = np.linspace(-10,20,1e3)
    x0 = 0
    sigma1 = 1
    b1 = gaussian(x, x0, sigma1)
    x0 = 5
    sigma2 = 2.5
    b2 = gaussian(x, x0, sigma2)
    x = np.linspace(-30, 30, 1e3)
    x0 = -2
    sigma3 = 10
    b3 = gaussian(x, x0, sigma3)



    cmap1 = cm.get_cmap('winter')
    cmap2 = cm.get_cmap('spring')
    cmap3 = cm.get_cmap('gray')
    cmap4 = cm.get_cmap('jet')

    # plt.figure()
    # plt.plot(b1, linewidth=3)
    # plt.savefig('./Gaussian_1.pdf')
    #
    # plt.figure()
    # plt.plot(b2, linewidth=3)
    # plt.savefig('./Gaussian_2.pdf')
    #
    # plt.figure()
    # plt.plot(b1, 'b', linewidth=3)
    # plt.plot(b2, 'b', linewidth=3)
    # plt.plot(b1+b2, 'k--', linewidth=3)
    # plt.savefig('./Gaussian_ncomp2_a.pdf')
    # plt.figure()
    # plt.plot(b1, 'b', linewidth=3)
    # plt.plot(b2, 'b', linewidth=3)
    # plt.plot(1.0/2.0*(b1 + b2), 'k--', linewidth=3)
    # plt.savefig('./Gaussian_ncomp2_b.pdf')
    #
    # plt.figure()
    # plt.plot(x, b1, 'b', linewidth=3)
    # plt.plot(x, b2, 'b', linewidth=3)
    # plt.plot(x, b3, 'b', linewidth=3)
    # plt.plot(x, b1 + b2 + b3, 'k--', linewidth=3)
    # plt.savefig('./Gaussian_ncomp3_a.pdf')
    #
    # plt.figure()
    # plt.plot(b1, 'b', linewidth=3)
    # plt.plot(b2, 'b', linewidth=3)
    # plt.plot(b3, 'b', linewidth=3)
    # plt.plot(1.0/3.0*(b1 + b2 + b3), 'k--', linewidth=3)
    # plt.savefig('./Gaussian_ncomp3_b.pdf')

    plt.figure()
    col1 = '0.0'
    col2 = '0.3'
    col3 = '0.6'
    plt.plot(b1, 'b', linewidth=1)
    plt.plot(b2, color=col1, linewidth=3, label='single Gaussian')
    plt.plot(b3, 'b', linewidth=1)
    plt.plot(1.0 / 2.0 * (b1 + b3), color=col2, linewidth=3, label='double Gaussian')
    plt.plot(1.0 / 3.0 * (b1 + b2 + b3), color=col3, linewidth=3, label='triple Gaussian')
    # plt.plot(1.0 / 3.0 * (b1 + b2 + b3), 'r', linewidth=3)
    plt.legend()
    plt.savefig('./Gaussian_ncomp2_ncomp3_a.pdf')

    plt.figure()
    plt.plot(b1)


    return



if __name__ == "__main__":
    main()