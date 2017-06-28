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
import copy


label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 6
plt.rcParams['figure.titlesize'] = 18
plt.rcParams['lines.linewidth'] = 0.8


def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("time_field")
    parser.add_argument("--type", nargs='+', type=str)
    args = parser.parse_args()
    path = args.path
    path_out_ = os.path.join(path, 'bulk_vs_single')
    case_name = args.casename
    time_field = args.time_field
    print('')
    print('time: ' + time_field)
    print('')
    print('path out: ' + path_out_)

    if args.type:
        type_list = args.type
    else:
        type_list = ['Couvreux']

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
    krange = np.arange(0,nz)
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
    path_files = os.path.join(path, 'fields', time_field + '.nc')
    print(path_files)
    thetal_flag = False
    buoyancy_flag = False
    # root = nc.Dataset(path_files, 'r')
    # qt_ = root.groups['fields'].variables['qt'][:, :, :]
    # ql_ = root.groups['fields'].variables['ql'][:, :, :]
    # w_ = root.groups['fields'].variables['w'][:, :, :]
    # try:
    #     buoyancy_ = root.groups['fields'].variables['buoyancy'][:, :, :]
    #     buoyancy_flag = True
    # except:
    #     buoyancy_flag = False
    #     print('!! buoyancy not in fields')
    # try:
    #     thetal_ = root.groups['fields'].variables['thetali'][:, :, :]
    #     thetal_flag = True
    # except:
    #     # thetal = np.zeros(shape=qt_.shape)
    #     # for i in range(nx):
    #     #     for j in range(ny):
    #     #         for k in range(nz):
    #     #             aux = thetali_c(p_ref[iz + k_], T_[i, j, iz + k_], qt_[i, j, iz + k_], ql_[i, j, iz + k_], qi_, Lv)
    #     #             theta_l = np.append(theta_l, aux)
    #     thetal_flag = False
    #     print('!! thetal not in fields')
    # root.close()



    # (2) read in Labels
    for type_ in type_list:
        print('')
        print('--- ' + type_ + ' ---')
        print('')
        if type_ == 'PDF':
            path_labels = os.path.join(path, 'Updrafts', 'Labeling_t' + str(time_field) + '.nc')
            print(path_labels)
            # root = nc.Dataset(path_labels, 'r')
            # labels = root.groups['fields'].variables['labels'][:, :, :]
        else:
            path_labels = os.path.join(path, 'tracer_fields')
            labels = read_in_updrafts_colleen(type_, time_field, path_labels)
            # read_in_updrafts_colleen(type_, time_field, path_labels)



    #     # (3) select the datapoints that are environemnt / updraft
        ini = time.time()
        qt = np.zeros(shape=(nx * ny, nk))
        qt_env = np.zeros(shape=(nx * ny, nk))
        qt_up = np.zeros(shape=(nx * ny, nk))
        ql = np.zeros(shape=(nx * ny, nk))
        ql_up = np.zeros(shape=(nx * ny, nk))
        ql_env = np.zeros(shape=(nx * ny, nk))
        w = np.zeros(shape=(nx * ny, nk))
        w_env = np.zeros(shape=(nx * ny, nk))
        w_up = np.zeros(shape=(nx * ny, nk))
        if thetal_flag:
            thetal = np.zeros(shape=(nx * ny, nk))
            thetal_env = np.zeros(shape=(nx * ny, nk))
            thetal_up = np.zeros(shape=(nx * ny, nk))
        if buoyancy_flag:
            buoy = np.zeros(shape=(nx * ny, nk))
            buoy_env = np.zeros(shape=(nx * ny, nk))
            buoy_up = np.zeros(shape=(nx * ny, nk))

        n_up = np.zeros(nk, dtype=np.int32)
        n_env = np.zeros(nk, dtype=np.int32)


        # labels_updrafts = np.zeros(shape=labels.shape, dtype=np.int)#copy.copy(labels)
        # labels_updraft_roots = np.zeros(shape=labels.shape, dtype=np.int)  # copy.copy(labels)
        # labels_updrafts[0,:,:] = labels[0,:,:]
        # labels_updrafts[:,0,:] = labels[:,0,:]
        # labels_updrafts[:,:,0] = labels[:,:,0]
        # labels_updrafts_ = copy.copy(labels_updrafts)
        # labels_updrafts_2 = copy.copy(labels_updrafts)
        # labels_updrafts_3 = copy.copy(labels_updrafts)
        #
        # label_root = 0
        # label = 1
        # root_list = []
        #
        # m0 = 40
        # n0 = 40
        # k0 = 5
        # print(np.nonzero(labels[0:m0,0:n0,k0]))
        # print('labels', labels[0:m0,0:n0,k0])
        # print('labels updrafts', labels_updrafts[0:m0, 0:n0, k0])

    #     # for i in range(1,nx-1):
    #     for i in range(1, m0):
    #         ishift = i * ny
    #         # for j in range(1,ny-1):
    #         for j in range(1, n0):
    #             ij = ishift + j
    #             # for k in range(1,nk-1):
    #             for k in range(5, 6):
    #                 iz = krange[k]
    #                 # iz = k
    #                 qt[ij, k] = qt_[i, j, iz]
    #                 ql[ij, k] = ql_[i, j, iz]
    #                 w[ij, k] = w_[i, j, iz]
    #                 if thetal_flag:
    #                     thetal[ij, k] = thetal_[i, j, iz]
    #                 if buoyancy_flag:
    #                     buoy[ij, k] = buoyancy_[i, j, iz]
    #
    #                 # (a) decompose fields
    #                 if labels[i, j, iz] == 0:
    #                     n_env[k] += 1
    #                     qt_env[ij,k] = qt_[i,j,iz]
    #                     ql_env[ij,k] = ql_[i,j,iz]
    #                     w_env = w_[i,j,iz]
    #                     if thetal_flag:
    #                         th_env = thetal_[i,j,iz]
    #                     if buoyancy_flag:
    #                         buoy_env = buoyancy_[i,j,iz]
    #                 elif labels[i, j, iz] == 1:
    #                     print('')
    #                     print('updraft: ', i,j,iz)
    #                     n_up[k] += 1
    #                     qt_up[ij, k] = qt_[i, j, iz]
    #                     ql_up[ij, k] = ql_[i, j, iz]
    #                     w_up = w_[i, j, iz]
    #                     if thetal_flag:
    #                         th_up = thetal_[i, j, iz]
    #                     if buoyancy_flag:
    #                         buoy_up = buoyancy_[i, j, iz]
    #
    #                     # (b) find root coordinates
    #                     if labels[i-1,j,iz] == 0 and labels[i,j-1,iz] == 0 and labels[i,j,iz-1] == 0:
    #                         label_root += 1
    #                         root = (i,j,iz, label_root)
    #                         root_list.append(root)
    #                         # print('root', type(root), root_list, type(root_list), type(root_list[-1]))
    #
    #
    #                     # (c) find connected updraft points
    #                     #        (i) assign correct label to point
    #                     l = labels[i, j, iz]
    #                     # if labels[i-1, j, iz] == 0 and labels[i, j - 1, iz] == 0 and labels[i, j, iz - 1] == 0:
    #                     # if labels[i-1, j, iz] == 0 and labels[i, j-1, iz] and labels[i, j+1, iz] == 0 and labels[i+1, j, iz] == 0:
    #                     if labels[i-1, j, iz] == 0 and labels[i, j-1, iz] == 0:
    #                         label += 1
    #                         l = label
    #                         l2 = label
    #                         if labels[i + 1, j, iz] == 1:
    #                             print('l=10', i, j, iz)
    #                             l2 = 10
    #                         elif labels[i, j + 1, iz] == 1:
    #                             print('l=20', i, j, iz)
    #                             l2 = 20
    #                     else:
    #                         if labels[i-1, j, iz] == 1:
    #                             l = labels_updrafts[i-1, j, iz]
    #                             l2 = labels_updrafts[i - 1, j, iz]
    #                             # l = 5
    #                             # print(i, j, iz, k0,'2', l)
    #                         elif labels[i,j-1, iz] == 1:
    #                             l = labels_updrafts[i, j-1, iz]
    #                             l2 = labels_updrafts[i, j - 1, iz]
    #                         # elif labels[i, j, iz-1] == 1:
    #                         #     l = labels_updrafts[i, j, iz-1]
    #
    #
    #
    #                     labels_updrafts[i,j,iz] = l
    #                     labels_updrafts_2[i, j, iz] = l2
    #                     if labels[i-1, j, iz] == 10:
    #                         labels_updrafts_2[i-1,j,iz] = l
    #                     elif labels[i, j-1, iz] == 20:
    #                         labels_updrafts_2[i,j-1, iz] = l
    #
    #                     print('updraft: ', i, j, iz, l)
    #
    #                     # print('labels', labels[0:10, 0:10, k0])
    #                     # print('labels updrafts', labels_updrafts[0:m0, 0:n0, k0])
    #
    #
    #
    #                     # if labels[i-1:i+2, j, iz].all() == 0 and labels[i, j - 1:j + 2, iz].all() == 0:
    #                     #     label += 1
    #                     #     l = label
    #                     #     print(i, j, iz, k0, '1', l)
    #                     # elif labels[i - 1, j, iz] == 1:
    #                     #     l = labels_updrafts_[i - 1, j, iz]
    #                     #     print(i, j, iz, k0, '2', l)
    #                     # elif labels[i, j - 1, iz] == 1:
    #                     #     l = labels_updrafts_[i, j - 1, iz]
    #                     #     print(i, j, iz, k0, '3', l)
    #                     # # elif labels[i, j, iz-1] == 1:
    #                     # #     l = labels_updrafts_[i, j, iz-1]
    #                     # #     print(i,j,iz,k0,'4', l)
    #                     # labels_updrafts_[i, j, iz] = l
    #                     #
    #                     # print('labels updrafts_', labels_updrafts_[0:m0, 0:n0, k0])
    #
    #     for i in range(m0-1,-1,-1):
    #     # for i in range(nx-2,-1,-1):
    #         for j in range(n0-1,-1,-1):
    #         # for j in range(ny-2,-1,-1):
    #             for k in range(6-1,5-1,-1):
    #             # for k in range(nk-2,-1,-1):
    #                 l = labels_updrafts[i, j, k]
    #                 if labels[i+1,j,iz] == 0 and labels[i,j-1,iz] == 0:
    #                     l = labels_updrafts[i,j,k]
    #                 elif labels[i+1,j,iz] != 0:
    #                     l = labels_updrafts[i+1,j,iz]
    #                 elif labels[i, j+1, iz] != 0:
    #                     l = labels_updrafts[i, j+1, iz]
    #                 # elif labels[i, j, iz+1] != 0:
    #                 #     l = labels_updrafts[i1, j, iz+1]
    #                 labels_updrafts_3[i,j,k] = l
    #
    #     print('roots:', label, len(root_list))
    #     print('')
    #
    #     # # (b) loop through all roots and label updraft points from there
    #     # # for i in len(root_list):
    #     # for n, val in enumerate(root_list):
    #     #     print('root ' + str(n) +': ', val)
    #     #     label = val[-1]
    #     #     x0 = val[0]
    #     #     y0 = val[1]
    #     #     z0 = val[2]
    #     #     i = x0-1
    #     #     j = y0-1
    #     #     z = z0-1
    #     #     # # print(x0,y0,z0)
    #     #     labels_updraft_roots[x0, y0, z0] = np.int(label)
    #     #     # # while(i<nx and j<ny and k<nz):
    #     #     # while(i<nx-1):
    #     #     #     i += 1
    #     #     #     while(j<ny-1):
    #     #     #         j += 1
    #     #     #         while(k<nz-1):
    #     #     #             k += 1
    #     #     #             print(i,j,k)
    #     #     #             if labels[i,j,k] == 1:
    #     #     #                 print('yes', label, i, j, k)
    #     #     #                 labels_updrafts[i,j,k] = np.int(label)
    #     #     #     i += 1
    #     #     #     j += 1
    #     #     #     k += 1
    #     #
    #     #     # for i in range(x0,nx):
    #     #     #     for j in range(y0,ny):
    #     #     #         # for k in range(z0,nz):
    #     #     #         k = z0
    #     #     #         while(k<nz):
    #     #     #             print(i, j, k)
    #     #     #             if labels[i, j, k] == 1:
    #     #     #                 print('yes', label, i, j, k)
    #     #     #                 labels_updrafts[i, j, k] = np.int(label)
    #     #     #                 k+=1
    #     #     #             else:
    #     #     #                 k = nz
    #     #
    #     # for i in range(nx):
    #     #     for j in range(ny):
    #     #         for k in range(nz):
    #     #             pass
    #     #
    #     #
    #     #
    #     #             # k = z0
    #     #             # while (k < nz):
    #     #             #     print(i, j, k)
    #     #             #     if labels[i, j, k] == 1:
    #     #             #         print('yes', label, i, j, k)
    #     #             #         labels_updrafts[i, j, k] = np.int(label)
    #     #             #         k += 1
    #     #             #     else:
    #     #             #         k = nz
    #     #
    #     #             # for i in range(x0,nx):
    #     #             #     for j in range(y0,ny):
    #     #             #         for k in range(z0,nz):
    #     #             #
    #     #             #             pass
    #     #             # l = len(val)
    #     #             # print(type(val), l, type(l))
    #     #             # print(val[0], val[1], val[2], label)
    #
    #
    #
    #     plt.figure(figsize=(20,5))
    #     plt.subplot(1,4,1)
    #     plt.imshow(labels[0:m0,0:n0,k0],interpolation="nearest")
    #     plt.colorbar(shrink=0.75)
    #     plt.title('labels')
    #     plt.xlabel('y')
    #     plt.ylabel('x')
    #     plt.subplot(1, 4, 2)
    #     plt.imshow(labels_updrafts[0:m0,0:n0,k0],interpolation="nearest")
    #     plt.colorbar(shrink=0.75)
    #     plt.title('labels updrafts')
    #     plt.xlabel('y')
    #     plt.ylabel('x')
    #     plt.subplot(1, 4, 3)
    #     plt.imshow(labels_updrafts_2[0:m0,0:n0,k0],interpolation="nearest")
    #     plt.colorbar(shrink=0.75)
    #     plt.title('labels updrafts 2')
    #     plt.xlabel('y')
    #     plt.ylabel('x')
    #     plt.subplot(1, 4, 4)
    #     plt.imshow(labels_updrafts_3[0:m0, 0:n0, k0], interpolation="nearest")
    #     plt.colorbar(shrink=0.75)
    #     plt.title('labels updrafts 3')
    #     plt.xlabel('y')
    #     plt.ylabel('x')
    #     plt.show()
    #     plt.close()


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

    # print('')
    # path = os.path.join(path_, root + 'environment.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())

    path = os.path.join(path_, root + 'time_'+ str(t) + '_Grid.pkl')
    print(path)
    print('')
    print('')
    f = open(path, 'r')
    print(type(f))
    labels = pickle.load(f)
    # # labels = pickle.load(open(path))
    # return labels
    return


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

if __name__ == '__main__':
    main()