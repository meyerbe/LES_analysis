#!/bin/sh

#  Vis.py
#
#
#  Created by Meyer  Bettina on 01/04/16.
#
from h5py import File
import h5py
import numpy as np
import matplotlib.cm as cm
import pylab as plt

import glob
import json as  simplejson


import pickle
from matplotlib.backends.backend_pdf import PdfPages


#from Namelist import Namelist
#nml = NamelistDCBL()
#nml = Namelist()



#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
    global fullpath_out, file_name
    global t, dt, dx, dz
    case_name = 'Bomex'

    path_list = ['../bomex/161130_test/n24/2_full_old_EV12/']
    # path = '../bomex/161130_test/n24/2_QL_old_EV12/'
    for path in path_list:
        fullpath_out = path + 'vis/'
        print('fullpath_out', fullpath_out)
        scheme = 'QL, 2nd, TKE, CFL = 0.1'

        nml = simplejson.loads(open(path + case_name + '.in').read())
        # namelist_files = glob.glob(path +'*.in')
        # print(namelist_files)
        # for namelist in namelist_files:
        #     nml = simplejson.loads(open(namelist).read())
        # dt = nml['stats_io']['frequency']
        dt = nml['visualization']['frequency']
        dx = nml['grid']['dx']
        dz = nml['grid']['dz']
        print('dt:', dt, dz)

        T = np.linspace(0,2400,5)
        print('T',T)
        for t in T:
            if t < 10:
                file_name = np.str(1000000) + np.str(np.int(t))
            elif t < 100:
                file_name = np.str(100000) + np.str(np.int(t))
            elif t < 1000:
                file_name = np.str(10000) + np.str(np.int(t))
            elif t < 10000:
                file_name = np.str(1000) + np.str(np.int(t))
            else:
                file_name = np.str(100) + np.str(np.int(t))
            print('name:', file_name)
            fullpath_in = fullpath_out + file_name  + '.pkl'
            print('fullpath_in',fullpath_in)
            # ----------------------------------------------
            # input = open(fullpath_in)
            # pickle.load(fullpath_in,'r')
            # import pickle
            try:
                with open(fullpath_in, 'rb') as f:
                    restart_data = pickle.load(f)

                f = open(fullpath_in)
                data = pickle.load(f)

                var_name = 'phi'
                var = data[var_name]
                levels = np.linspace(-1e-9, 1+1e-9, 250)
                # plot_data_levels(var, var_name, levels)
                plot_data(var, var_name)

                var_name = 'w'
                f = open(fullpath_in)
                data = pickle.load(f)
                # print('data:', data[var_name].shape)
                var = data[var_name]
                levels = np.linspace(-10,10,100)
                plot_data(var,var_name)

                var_name = 'potential_temperature'
                f = open(fullpath_in)
                data = pickle.load(f)
                # print('data:', data[var_name].shape)
                var = data[var_name]
                levels = np.linspace(-10, 10, 100)
                plot_data(var, var_name)

                var_name = 's'
                f = open(fullpath_in)
                data = pickle.load(f)
                # print('data:', data[var_name].shape)
                var = data[var_name]
                levels = np.linspace(-1e-6, 1, 100)
                # print('levels:', levels)
                plot_data(var, var_name)

                var_name = 'ql'
                f = open(fullpath_in)
                data = pickle.load(f)
                # print('data:', data[var_name].shape)
                var = data[var_name]
                levels = np.linspace(-1e-6, 1, 100)
                # print('levels:', levels)
                plot_data(var, var_name)
            except:
                continue





    #data = (1.4,42)
    #output = open('data.pkl', 'w')
    #pickle.dump(data, output)
    #output.close()


    #list_of_names = ['u', 'v', 'w', 'specific_entropy']#['u', 'v', 'w']
    #list_of_names = ['u']
    #for name in list_of_names:
    #    data = read_in(name, 'fields', fullpath_in)       # type == 'profiles', 'fields'

    print('Ende')


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
def read_in(variable_name, group_name, fullpath_in):
    f = File(fullpath_in)

    # Get access to the profiles group
    profiles_group = f[group_name]
    # Get access to the variable dataset
    variable_dataset = profiles_group[variable_name]
    # Get the current shape of the dataset
    variable_dataset_shape = variable_dataset.shape

    variable = np.ndarray(shape=variable_dataset_shape)
    for t in range(variable_dataset_shape[0]):
        if group_name == "timeseries":
            variable[t] = variable_dataset[t]
        elif group_name == "profiles":
            variable[t, :] = variable_dataset[t, :]
        elif group_name == "correlations":
            variable[t, :] = variable_dataset[t, :]
        elif group_name == "fields":
            variable[t] = variable_dataset[t]

    f.close()
    return variable


# ----------------------------------------------------------------------
def plot_data(data, var_name):
    print(data.shape)
    plt.figure()
    # plt.contourf(data.T, levels = level)
    plt.contourf(data.T)
    # plt.show()
    plt.colorbar()
    plt.title(var_name + ', (t=' + np.str(t) + 's)')
    plt.xlabel('x (dx=' + np.str(dx) + 'm)')
    plt.ylabel('height z (dz=' + np.str(dz) + 'm)')
    # plt.savefig(fullpath_out + 'phi/' + var_name + '_' + file_name + '.png')
    plt.savefig(fullpath_out + var_name + '_' + file_name + '.png')
    plt.title(var_name)


def plot_data_levels(data, var_name, level):
    print(data.shape)
    plt.figure()
    plt.contourf(data.T, levels=level)
    # plt.contourf(data.T)
    # plt.show()
    plt.colorbar()
    # plt.savefig(fullpath_out + 'phi/' + var_name + '_' + file_name + '.png')
    plt.savefig(fullpath_out + var_name + '_a' + file_name + '.pdf')
    plt.title(var_name)


def output(fullpath_in, fullpath_out, file_name, var_name):
    print('out')

# ----------------------------------------------------------------------

if __name__ == "__main__":
    main()