#!/bin/sh

#  Vis.py
#
#
#  Created by Meyer  Bettina on 01/04/16.
#

import argparse
import os
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

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 35
plt.rcParams['lines.linewidth'] = 3

#----------------------------------------------------------------------
#----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("--var_name")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    if args.var_name:
        var_list = [args.var_name]
    else:
        var_list = ['w', 's', 'potential_temperature', 'temperature', 'ql', 'qt', 'u', 'v']
    print(var_list)
    print(' ')

    global fullpath_out, file_name
    global t, dt, dx, dz
    # case_name = 'Bomex'

    # path_list = ['../bomex/161130_test/n24/2_full_old_EV12/']
    # path = '../bomex/161130_test/n24/2_QL_old_EV12/'
    path_list = [path]
    for path in path_list:
        fullpath_out = os.path.join(path,'vis/')
        print('fullpath_out', fullpath_out)
        scheme = 'QL, 2nd, TKE, CFL = 0.1'

        nml = simplejson.loads(open(os.path.join(path,case_name + '.in')).read())
        # namelist_files = glob.glob(path +'*.in')
        # print(namelist_files)
        # for namelist in namelist_files:
        #     nml = simplejson.loads(open(namelist).read())
        dt = nml['visualization']['frequency']
        dx = nml['grid']['dx']
        dz = nml['grid']['dz']
        print('vis dt:' + str(dt) +', dz: '+ str(dz))

        # files = os.listdir(os.path.join(path, 'vis/*.pkl'))       # type = list
        files = os.listdir(os.path.join(path, 'vis'))  # type = list
        print('visualisation files: ', str(files))
        print('')


        # T = np.linspace(0,2400,5)
        # print('T',T)
        # for t in T:
        #     if t < 10:
        #         file_name = np.str(1000000) + np.str(np.int(t))
        #     elif t < 100:
        #         file_name = np.str(100000) + np.str(np.int(t))
        #     elif t < 1000:
        #         file_name = np.str(10000) + np.str(np.int(t))
        #     elif t < 10000:
        #         file_name = np.str(1000) + np.str(np.int(t))
        #     else:
        #         file_name = np.str(100) + np.str(np.int(t))
        #     print('name:', file_name)
        #     fullpath_in = fullpath_out + file_name  + '.pkl'
        #     print('fullpath_in: ' + fullpath_in)

        for file_name in files:
            if file_name[-4:] == '.pkl':
                t = np.int(file_name[0:-4])
                fullpath_in = fullpath_out + file_name

                f = open(fullpath_in)
                data = pickle.load(f)
                for var_name in var_list:
                    print(var_name+ ', fullpath_in: ' + fullpath_in)
                    try:
                        var = data[var_name]
                    except:
                        print("No variable " + var_name)
                        print('')
                        continue

                        # if var_name == 'phi':
                        #     levels = np.linspace(-1e-9, 1 + 1e-9, 250)
                        # elif var_name == 'w':
                        #     levels = np.linspace(-10, 10, 100)
                        # elif var_name == 'potential_temperature':
                        #     levels = np.linspace(-10, 10, 100)
                        # elif var_name == 's':
                        #     levels = np.linspace(-1e-6, 1, 100)
                        # elif var_name == 'ql':
                        #     levels = np.linspace(-1e-6, 1, 100)
                        # elif var_name == 'qt':
                        #     levels = np.linspace(0.0, 0.02, 100)
                        # else:
                        #     levels = np.linspace(-10, 10, 100)
                    levels = np.linspace(np.amin(var), np.amax(var), 100)

                    plot_data_levels(var, var_name, t, levels)
                    plot_data(var, var_name, t)



            # # ----------------------------------------------
            # # input = open(fullpath_in)
            # # pickle.load(fullpath_in,'r')
            # # import pickle
            # try:
            #     # with open(fullpath_in, 'rb') as f:
            #     #     restart_data = pickle.load(f)
            #
            #     f = open(fullpath_in)
            #     data = pickle.load(f)
            #
            #     var_name = 'phi'
            #     var = data[var_name]
            #     levels = np.linspace(-1e-9, 1+1e-9, 250)
            #     # plot_data_levels(var, var_name, levels)
            #     plot_data(var, var_name)
            #
            #     var_name = 'w'
            #     f = open(fullpath_in)
            #     data = pickle.load(f)
            #     # print('data:', data[var_name].shape)
            #     var = data[var_name]
            #     levels = np.linspace(-10,10,100)
            #     plot_data(var,var_name)
            #
            #     var_name = 'potential_temperature'
            #     f = open(fullpath_in)
            #     data = pickle.load(f)
            #     # print('data:', data[var_name].shape)
            #     var = data[var_name]
            #     levels = np.linspace(-10, 10, 100)
            #     plot_data(var, var_name)
            #
            #     var_name = 's'
            #     f = open(fullpath_in)
            #     data = pickle.load(f)
            #     # print('data:', data[var_name].shape)
            #     var = data[var_name]
            #     levels = np.linspace(-1e-6, 1, 100)
            #     # print('levels:', levels)
            #     plot_data(var, var_name)
            #
            #     var_name = 'ql'
            #     f = open(fullpath_in)
            #     data = pickle.load(f)
            #     # print('data:', data[var_name].shape)
            #     var = data[var_name]
            #     levels = np.linspace(-1e-6, 1, 100)
            #     # print('levels:', levels)
            #     plot_data(var, var_name)
            #
            #     # var_name = 'qt'
            #     # print(var_name)
            #     # # print('data:', data[var_name].shape)
            #     # var = data[var_name]
            #     # levels = np.linspace(0.0,0.02, 100)
            #     # # print('levels:', levels)
            #     # plot_data(var, var_name)
            #     # plot_data_levels(var,var_name,levels)
            # except:
            #     print('no file found')
            #     print('')
            #     continue


    # #list_of_names = ['u', 'v', 'w', 'specific_entropy']#['u', 'v', 'w']
    # #list_of_names = ['u']
    # #for name in list_of_names:
    # #    data = read_in(name, 'fields', fullpath_in)       # type == 'profiles', 'fields'

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
def plot_data(data, var_name, t):
    # print(data.shape)
    plt.figure()
    if var_name == 'w':
        plt.contourf(data.T,cmap = cm.bwr)
    else:
        plt.contourf(data.T)
    # plt.show()
    plt.colorbar()
    plt.title(var_name + ', (t=' + np.str(t) + 's)')
    plt.xlabel('x (dx=' + np.str(dx) + 'm)')
    plt.ylabel('height z (dz=' + np.str(dz) + 'm)')
    plt.savefig(os.path.join(fullpath_out, 'pdf-pics', var_name + '_' + str(t) + '.pdf'))
    plt.savefig(fullpath_out + var_name + '_' + str(t) + '.png')
    plt.close()
    return


def plot_data_levels(data, var_name, t, levels_):
    # print(data.shape)
    plt.figure()
    if var_name == 'w':
        # print('bwr')
        # plt.contourf(data.T, levels=levels_, cmap=plt.cm.bwr)
        plt.contourf(data.T, cmap=cm.bwr, levels=levels_)
    else:
        plt.contourf(data.T, levels=levels_, cmap=cm.viridis)
    # plt.contourf(data.T)
    # plt.show()
    plt.title(var_name + ', (t=' + np.str(t) + 's)')
    plt.xlabel('x (dx=' + np.str(dx) + 'm)')
    plt.ylabel('height z (dz=' + np.str(dz) + 'm)')
    plt.colorbar()
    plt.savefig(os.path.join(fullpath_out, 'pdf-pics', 'levels_' + var_name + '_' + str(t) + '.pdf'))
    plt.savefig(fullpath_out + 'levels_' + var_name + '_' + str(t) + '_levels.png')
    plt.close()
    return



def plot_data_levels_cont(var_field,var_name_field,levels_field,var_cont,var_name_cont,levels_cont):
    print('cont print')
    plt.figure()
    if var_name_field == 'w':
        ax1 = plt.contourf(var_field.T, cmap=cm.bwr, levels=levels_field)
    else:
        ax1 = plt.contourf(var_field.T, levels=levels_field)
    ax2 = plt.contour(var_cont.T, levels=levels_cont,colors='k')
    ax3 = plt.contour(var_cont.T, levels=[0.018,0.02], colors='b',linewidths=1)
    plt.clabel(ax2,inline=1)
    plt.colorbar(ax1,shrink=0.8)
    plt.colorbar(ax3,shrink=0.8)
    plt.title(var_name_field + ', (t=' + np.str(t) + 's)')
    plt.xlabel('x (dx=' + np.str(dx) + 'm)')
    plt.ylabel('height z (dz=' + np.str(dz) + 'm)')
    # plt.savefig(fullpath_out + 'phi/' + var_name + '_' + file_name + '.png')
    print('hoi')
    path = fullpath_out + var_name_field + '_cont_' + file_name + '.png'
    print('hohioi')
    print('out: ', fullpath_out + var_name_field + '_cont_' + file_name + '.png')
    plt.savefig(fullpath_out + var_name_field + '_cont_' + file_name + '.png')
    plt.close()


def output(fullpath_in, fullpath_out, file_name, var_name):
    print('out')



# ----------------------------------------------------------------------

if __name__ == "__main__":
    main()