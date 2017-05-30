# (0) chose time
# (1) import Fields
#     --> load field data for a advected variable phi and all velocities u,v,w
# (2) import Statistical File (nc-file)
#     (mean Profiles --> horizontal domain mean)
#     --> load mean-profile for advected variable phi and all velocities u,v,w
#     --> chose array[var,z] at time t
# (3) phi' = phi - mean[phi], u' = ... etc.
import netCDF4 as nc
import argparse
import os, sys
import numpy as np
import json as  simplejson
import matplotlib.cm as cm
import pylab as plt


label_size = 18
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = label_size
plt.rcParams['figure.titlesize'] = 70
plt.rcParams['lines.linewidth'] = 2


def main():
    global case_name
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("--var_name")
    parser.add_argument("--cont_name")
    parser.add_argument("--files")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    if args.var_name:
        var_list = [args.var_name]
    else:
        if case_name == 'TRMM_LBA':
            var_list = ['ql', 'qt', 'qr', 'w', 's', 'temperature', 'u', 'v']
        elif case_name == 'DYCOMS_RF01':
            var_list = ['w', 's', 'thetali', 'temperature', 'ql', 'qt', 'u', 'v']
            var_list = ['ql', 'qt', 'w', 's', 'temperature', 'u', 'v']
        else:
            var_list = ['ql', 'qt', 'w', 's', 'thetali', 'temperature', 'u', 'v']

    if args.cont_name:
        cont_list = [args.cont_name]
    else:
        cont_list = ['qt', 'ql']
    # var_list_corr = ['wphi','uphi', 'vphi', 'uphi_div', 'vphi_div', 'wphi_div']
    print('var list:', var_list)
    print('contours list:', cont_list)
    print('')



    # -----------
    global fullpath_out
    fullpath_out = os.path.join(path,'fields_figures/')
    print('fullpath_out: ', fullpath_out)
    # T = [1800,3600,5400]

    path_fields = os.path.join(path,'fields/')
    files = os.listdir(path_fields)  # type = list
    # files = [files[1]]
    print('Found the following fields: ', str(files), type(files), type(files[0]), files[0], len(files))
    print('')

    # -----------


    # (0) import Namelist --> to chose right mean profile, fitting with time
    global nx0, ny0, nz0
    global time, nml
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    dt = nml['stats_io']['frequency']
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    print('dt:', dt, 'dz:', dz)

    # (1) Read in Test field
    field = read_in_netcdf_fields('w',os.path.join(path_fields,files[0]))

    # (2) read in grid dimensions
    #       ni:                   field dimensions
    #       nx0, ny0, nz0:        coordinates of center
    global n, ntot
    ni_ = np.zeros((3,))
    n = np.zeros((3), dtype=np.int16)
    # n = n.astype(int)
    ni_[0] = nml['grid']['nx']
    ni_[1] = nml['grid']['ny']
    ni_[2] = nml['grid']['nz']
    for i in range(3):
        n[i] = field.shape[i]
        if n[i] != ni_[i]:
            print('Dimensions do not fit!')
            sys.exit()
    print(n[0], n[1], n[2])
    nx0 = np.int(n[0]/2)
    ny0 = np.int(n[1]/2)
    nz0 = np.int(n[2]/2)

    # ntot = n[0]*n[1]*n[2]
    print('x0,y0,z0', nx0, ny0, nz0, n[:])



    # (3) Set Levels and Times
    zrange = np.asarray(np.linspace(1, n[2] - 1, 10), dtype=np.int16)
    # zrange = np.asarray([10])
    yrange = np.asarray(np.linspace(1, n[1] - 1, 10), dtype=np.int16)
    # yrange = np.asarray([10])
    zrange, files = set_zrange(case_name)
    if args.files:
        files = ['']
        files[0] = args.files + '.nc'




    # (4) Visualize Fields
    print('')
    print('--- plot crosssections with contours ---')
    for file_name in files:
        if file_name[-3:] == '.nc' and file_name[0] != '0':
            print(file_name)
            tt = np.int(file_name[0:-3])
            for var_name in var_list:
                print('')
                print(var_name+ ', path: ' + os.path.join(path_fields, file_name))
                try:
                    field_data = read_in_netcdf_fields(var_name, os.path.join(path_fields, file_name))
                except:
                    print('!! Problem reading in '+ var_name + ' in '+os.path.join(path_fields, file_name))
                    print " "
                    sys.exit()
                    continue

                levels = np.linspace(np.amin(field), np.amax(field), 100)

                # (a) Visualize Fields & Fields plus 1 Contour
                for cont_name in cont_list:
                    print('cont: ', cont_name+ ', path: ' + os.path.join(path_fields, file_name))
                    # cont_name = 'phi'
                    # cont_name = 'qt'
                    if cont_name != var_name and cont_name != ' ':
                        cont_data = read_in_netcdf_fields(cont_name, os.path.join(path_fields, file_name))
                    print('Contour variable: ' + cont_name + ' in ' + os.path.join(path_fields, file_name))
                    for k in zrange:
                        plot_name = var_name + '_t'+ np.str(tt) + '_z' + np.str(np.int(k * dz)) + 'm'
                        plot_field(var_name, field_data[:, :, k], k*dz, plot_name, tt, 'hor')
                        if cont_name != var_name and cont_name != ' ':
                            plot_name = var_name + '_' + cont_name + '-cont_t' + np.str(tt) + '_z' + np.str(np.int(k*dz)) + 'm'
                            plot_field_cont(var_name, field_data[:, :, k], cont_name,cont_data[:, :, k], k*dz, plot_name, tt, 'hor')
                    for k in yrange:
                        plot_name = var_name + '_t'+ np.str(tt) + '_y' + np.str(np.int(k * dy)) + 'm'
                        plot_field(var_name, field_data[:, k, :], k*dy, plot_name, tt, 'vert')
                        if cont_name != var_name and cont_name != ' ':
                            plot_name = var_name + '_' + cont_name + '-cont_t' + np.str(tt) + '_y' + np.str(np.int(k*dy)) + 'm'
                            plot_field_cont(var_name, field_data[:, k, :], cont_name,
                                            cont_data[:, k, :], k*dz, plot_name, tt, 'vert')

                # (b) Visualize Fields plus 2 Contours
                # cont_name1 = 'qt'
                # cont_name2 = 'w'
                # cont_data1 = read_in_netcdf_fields(cont_name1, os.path.join(path_fields, file_name))
                # cont_data2 = read_in_netcdf_fields(cont_name2, os.path.join(path_fields, file_name))
                # for k in zrange:
                #     plot_name = var_name + '_' + cont_name1 + '-cont_' + cont_name2 + '-cont_t' \
                #                 + np.str(tt) + '_z' + np.str(np.int(k * dz)) + 'm'
                #     plot_field_cont1_cont2(var_name, field_data[:,:,k],
                #                                cont_name1,cont_data1[:,:,k],cont_name2,cont_data2[:,:,k], plot_name)

        else:
            print('!!!', file_name)


    return




# ----------------------------------
def set_zrange(case_name):
    if case_name[0:8] == 'ZGILS_S6':
        # ZGILS 6
        files = ['1382400.nc']
        # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
        # files = files[0:25:2]
        krange = np.asarray([25, 25, 40, 50, 60, 65], dtype=np.int32)
    elif case_name[0:9] == 'ZGILS_S12':
        # ZGILS S12
        # files = ['432000.nc']
        files = ['86400.nc']
        # files = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
        krange = np.asarray([35,40,45])
    elif case_name == 'DYCOMS_RF01':
        # DYCOMS RF01 large
        krange = np.asarray([140, 150, 160, 166, 180])
        # files = ['3600.nc']
        # DYCOMS RF01
        krange = np.asarray([140, 150, 160, 166, 180])
        # files = ['10800.nc', '12600.nc', '14400.nc']
        files = ['10800.nc']
    elif case_name == 'DYCOMS_RF02':
        # DYCOMS RF02
        # krange = np.asarray([120, 170])
        krange = np.asarray([120, 140, 150, 160, 170, 200])
        # files = ['18000.nc', '19800.nc', '21600.nc']
        # files = ['10800.nc']
        files = ['3600.nc']
    elif case_name == 'Bomex':
        # Bomex large, kyle
        # krange = np.asarray([27, 91])
        ## Bomex 170314_weno7 (dz=40)
        # krange = np.asarray([10, 12, 15, 18, 20, 22, 25, 40, 50])
        # files = ['21600.nc']
        ## Bomex (dz=20)
        krange = np.asarray([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 125], dtype=np.int32)
        files = ['21600.nc']
        files = ['18000.nc']
        files = ['19800.nc']
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
        krange = np.asarray([10, 15, 20, 30, 40, 50, 60])

    return krange, files





# ----------------------------------
def plot_field(field_name, field_data, level, file_name, tt, type):
    print('plot field: ', field_name)
    global nml
    dx = nml['grid']['dx']
    dy = nml['grid']['dy']
    dz = nml['grid']['dz']
    plt.figure(figsize=(12,10))
    if field_name == 'w':
        ax1 = plt.contourf(field_data.T, cmap=cm.bwr)
    else:
        ax1 = plt.contourf(field_data.T, cmap=cm.viridis)
    plt.colorbar(ax1)

    max_field = np.amax(field_data)
    plt.title(field_name + ', max:' + "{0:.2f}".format(max_field), fontsize=35)
    if type == 'hor':
        # plt.title(field_name + ', z='+ str(level) + 'm (max:' + "{0:.2f}".format(max_field), fontsize=28)
        plt.title(field_name + ', (z=' + str(level) + 'm, t=' + str(tt) + 's)', fontsize=28)
        plt.xlabel(r'x ($\Delta $x='+str(dx)+'m)')
        plt.ylabel(r'y ($\Delta $y='+str(dy)+'m)')
    elif type == 'vert':
        # plt.title(field_name + ', y=' + str(level) + 'm  (max:' + "{0:.2f}".format(max_field), fontsize=28)
        plt.title(field_name + ', y=' + str(level) + 'm, t=' + str(tt) + 's)', fontsize=28)
        plt.xlabel('x')
        plt.ylabel(r'z ($\Delta $z=' + str(dz) + 'm)')

    print('saving: ', fullpath_out + file_name+ '.png')
    plt.savefig(fullpath_out + file_name + '.png')
    # # plt.show()
    plt.close()
    return




def plot_field_cont(field_name, field_data,cont_name,cont_data, level, file_name, tt, type):
    print('plot field & cont: ', field_name, cont_name)
    plt.figure(figsize=(15,10))
    if field_name == 'w':
        ax1 = plt.contourf(field_data.T, cmap=cm.bwr)
    else:
        ax1 = plt.contourf(field_data.T, cmap=cm.viridis)
    # cont = np.linspace(1.0,1.1,11)
    max = np.amax(cont_data)
    min = np.amin(cont_data)
    if np.amax(np.abs(cont_data)) > 0.0 and max != min:
        cont = np.linspace(min, max, 11)
        ax2 = plt.contour(cont_data.T, cont, cmap=cm.Greys)
        plt.colorbar(ax2)
    plt.colorbar(ax1)

    max_field = np.amax(field_data)
    max_data = np.amax(cont_data)
    if type == 'hor':
        plt.title(field_name + ', z=' + str(level) + 'm, t=' + str(tt) + 's (contours: ' + cont_name + ', max: ' + "{0:.2f}".format(
            max_data) + ')', fontsize=28)
        plt.xlabel('x')
        plt.ylabel('y')
    elif type == 'vert':
        # plt.title(field_name + ', y=' + str(level) + 'm  , (contours: ' + cont_name + ', max: ' + "{0:.2f}".format(
        #     max_data) + ')', fontsize=35)
        plt.title(field_name + ', y=' + str(level) + 'm, t=' + str(tt) + 's (contours: ' + cont_name + ', max: ' + "{0:.2f}".format(
            max_data) + ')', fontsize=28)
        plt.xlabel('x')
        plt.ylabel('z')

    # print('saving: ', fullpath_out + file_name + '.png')
    # print('')
    plt.savefig(fullpath_out + file_name + '.png')
    # # plt.show()
    plt.close()
    return



def plot_field_cont1_cont2(field_name, field_data, cont_name1, cont_data1, cont_name2, cont_data2, file_name):
    # print('plot corr/cont: ', field_name, field_data.shape)
    plt.figure(figsize=(19,10))
    if field_name == 'w':
        levels=np.linspace(-6,6,250)
        ax1 = plt.contourf(field_data.T, cmap=cm.bwr, levels=levels)
    else:
        ax1 = plt.contourf(field_data.T)
    cont1 = np.linspace(0.93,1.01,9)
    ax2a = plt.contour(cont_data1.T, levels=cont1)
    cont2 = [-3.0,-2.0,-1.0,1.0,2.0,3.0]
    ax2b = plt.contour(cont_data2.T, levels=cont2, colors='k', linewidths=0.6)
    plt.colorbar(ax2a,shrink=0.75)
    plt.colorbar(ax2b,shrink=0.75)
    plt.colorbar(ax1)
    max_field = np.amax(field_data)
    max_data1 = np.amax(cont_data1)
    max_data2 = np.amax(cont_data2)
    plt.title(field_name+', max:'+"{0:.2f}".format(max_field)+', (contours: '+cont_name1+', max: '+"{0:.2f}".format(max_data1)
              +', '+cont_name2+', max: '+"{0:.2f}".format(max_data2)+')')
    plt.xlabel('x')
    plt.ylabel('z')
    # plt.savefig(fullpath_out + file_name + '.png')
    # print('')
    # print('saving: ', fullpath_out + file_name + '.png')
    # print('')
    plt.close()




# ----------------------------------
# def plot_data_vertical(data, var_name, file_name):
#     print('plot vertical: ', var_name, data.shape)
#     plt.figure()
#     ax1 = plt.contourf(data.T)
#     if var_name == 'phi':
#         cont = np.linspace(1.0,1.1,11)
#         ax2 = plt.contour(data.T, levels = cont)
#         plt.colorbar(ax2)
#     # plt.show()
#     plt.colorbar(ax1)
#     max = np.amax(data)
#     plt.title(var_name + ', max:' + "{0:.2f}".format(np.amax(data)), fontsize=12)
#     plt.xlabel('x')
#     plt.ylabel('z')
#     plt.savefig(fullpath_out + file_name + '.png')
#     plt.close()


# def plot_data_vertical_levels(data, var_name, level):
#     print(data.shape)
#     plt.figure()
#     plt.contourf(data.T, levels = level)
#     if var_name == 'phi':
#         cont = np.linspace(1.0,1.1,11)
#         ax2 = plt.contour(data.T, levels = cont)
#         plt.colorbar(ax2)
#     plt.colorbar(ax1)
#     plt.title(var_name + ', max:' + "{0:.2f}".format(np.amax(data)), fontsize=12)
#     plt.xlabel('x')
#     plt.ylabel('z')
#     plt.savefig(fullpath_out + file_name + '.png')
#     plt.close()


# def plot_data_horizontal(data, var_name):
#     print(data.shape)
#     plt.figure()
#     plt.contourf(data.T)
#     # plt.show()
#     plt.colorbar()
#     max = np.amax(data)
#     plt.title(var_name + ', max:' + "{0:.2f}".format(np.amax(data)))
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.savefig(fullpath_out + file_name + '.png')
#     plt.close()
#
#
# def plot_data_horizontal_levels(data, var_name, level):
#     print(data.shape)
#     plt.figure()
#     plt.contourf(data.T, levels = level)
#     # plt.show()
#     plt.colorbar()
#     max = np.amax(data)
#     plt.title(var_name + ', max:' + "{0:.2f}".format(np.amax(data)))
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.savefig(fullpath_out + file_name + '.png')
#     plt.close()


# ----------------------------------
def read_in_netcdf_fields(variable_name, fullpath_in):
    print('.....', fullpath_in, variable_name)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups['fields'].variables[variable_name][:, :, :]
    rootgrp.close()
    return var



def read_in_netcdf_profile_all(variable_name, group_name, fullpath_in):
    #    print(fullpath_in)
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]
    shape = var.shape
    data = np.ndarray(shape = var.shape)
    if group_name != 'profiles':
        var = rootgrp.groups[group_name].variables[variable_name]
    for t in range(shape[0]):
        if group_name == "profiles":
            data[t,:] = var[t, :]
            nkr = rootgrp.groups['profiles'].variables['z'].shape[0]
        if group_name == "correlations":
            data[t,:] = var[t, :]
        if group_name == "timeseries":
            data[t] = var[t]
    rootgrp.close()
    return data


def read_in_netcdf_profile(variable_name, group_name, fullpath_in):
    rootgrp = nc.Dataset(fullpath_in, 'r')
    var = rootgrp.groups[group_name].variables[variable_name]
    shape = var.shape
    #print('read_in_profile: ', time, var.shape, nt, type(nt))
    
    if group_name == "profiles":
        data = np.ndarray((shape[1],))
        data[:] = var[nt, :]
    if group_name == "correlations":
        data = np.ndarray((shape[1],))
        data[:] = var[nt, :]
    if group_name == "timeseries":
        data = var[nt]
    rootgrp.close()
    return data


# ----------------------------------



if __name__ == "__main__":
    main()