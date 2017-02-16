# author: Bettina Meyer
import argparse
import sys, os

import numpy as np
import pylab as plt
from math import fabs

label_size = 6
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 6
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['lines.linewidth'] = 2
# Legend
plt.rcParams['legend.loc'] = 'best'
plt.rcParams['legend.fontsize'] = 'small'
plt.rcParams['legend.framealpha'] = 0.8      # legend patch transparency
#legend.frameon       : True     # if True, draw the legend on a background patch
#legend.facecolor     : inherit  # inherit from axes.facecolor; or color spec
#legend.edgecolor     : 0.8      # background patch boundary color
#legend.fancybox      : True     # if True, use a rounded box for the
                                 # legend background, else a rectangle
#legend.numpoints     : 1        # the number of marker points in the legend line
#legend.scatterpoints : 1        # number of scatter points
#legend.markerscale   : 1.0      # the relative size of legend markers vs. original

# Dimensions as fraction of fontsize:
#legend.borderpad     : 0.4      # border whitespace
#legend.labelspacing  : 0.5      # the vertical space between the legend entries
#legend.handlelength  : 2.0      # the length of the legend lines
#legend.handleheight  : 0.7      # the height of the legend handle
#legend.handletextpad : 0.8      # the space between the legend line and legend text
#legend.borderaxespad : 0.5      # the border between the axes and legend edge
#legend.columnspacing : 2.0      # column separation

sys.path.append("..")
from io_read_in_files import read_in_nml, read_in_netcdf
from CC_thermodynamics import sat_adj_fromentropy, sat_adj_fromthetali
from thermodynamic_functions import theta_li

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    # ______________________
    global nx, dx
    case_name = args.casename
    dx, nx, dt = read_in_nml(args.path, case_name)
    global path
    path = args.path
    print('path:', path)
    # ______________________
    files = os.listdir(os.path.join(args.path, 'fields'))
    N = len(files)
    print('Found the following directories', files, N)
    print('')
    # ______________________
    global time
    time = np.zeros((1))
    for d in files:
        time = np.sort(np.append(time, np.int(d[0:-3])))
    print('time: ', time)
    print('')
    # ______________________
    ''' read in reference state '''
    if case_name == 'ZGILS6':
        fullpath_in_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub_D15.nc')
    else:
        fullpath_in_ref = os.path.join(path, 'Stats.'+case_name+'.nc')
    print(fullpath_in_ref)
    try:
        p_ref = read_in_netcdf('p0_half', 'reference', fullpath_in_ref)
    except:
        p_ref = read_in_netcdf('p0', 'reference', fullpath_in_ref)
    z_ref = read_in_netcdf('z', 'reference', fullpath_in_ref)
    z_half_ref = read_in_netcdf('z_half', 'reference', fullpath_in_ref)
    print('')
    print('nz, dz: ', nx[2], dx[2])
    print('height: ', z_ref, z_ref.shape)
    print(z_half_ref)
    # plot_profiles(fullpath_in_ref, z_ref)

    # ______________________
    ''' zrange '''
    global zrange
    zrange = [15,18,20]
    # ______________________
    # ______________________

    # ______________________

    files_ = [files[0]]
    for d in files_:
        print('')
        print(files, d, files_)
        nc_file_name = str(d)
        fullpath_in = os.path.join(path, 'fields', nc_file_name)

        '''read in fields'''
        s_ = read_in_netcdf('s', 'fields', fullpath_in)
        qt_ = read_in_netcdf('qt', 'fields', fullpath_in)
        T_ = read_in_netcdf('temperature', 'fields', fullpath_in)
        ql_ = read_in_netcdf('ql', 'fields', fullpath_in)
        alpha_ = np.zeros(shape=ql_.shape)
        T_comp = np.zeros(shape=ql_.shape)
        ql_comp = np.zeros(shape=ql_.shape)
        alpha = np.zeros(shape=ql_.shape)
        theta_l = np.zeros(shape=ql_.shape)
        T_comp_thl = np.zeros(shape=ql_.shape)
        ql_comp_thl = np.zeros(shape=ql_.shape)
        alpha_thl = np.zeros(shape=ql_.shape)
        for i in range(nx[0]):
            for j in range(nx[1]):
                for k in range(nx[2]):
                    if ql_[i,j,k] != 0.0:
                        alpha_[i,j,k] = 1

        max_T_sat = 0.0
        max_T_unsat = 0.0
        max_ql = 0.0
        min_ql = 0.0

        max_T_sat_thl = 0.0
        max_T_unsat_thl = 0.0
        max_ql_thl = 0.0
        min_ql_thl = 0.0
        # for k in range(nx[2]):
        for k in zrange:
            for i in range(nx[0]):
                for j in range(nx[1]):
                    # T_comp, ql_comp = sat_adj(p, s[i,j], qt[i,j])
                    T_comp[i, j, k], ql_comp[i, j, k], alpha[i, j, k] = sat_adj_fromentropy(p_ref[k], s_[i, j, k],qt_[i, j, k])

                    theta_l[i, j, k] = theta_li(p_ref[k], T_[i, j, k], qt_[i, j, k], ql_[i, j, k], 0)
                    T_comp_thl[i, j, k], ql_comp_thl[i, j, k], alpha_thl[i, j, k] = sat_adj_fromthetali(p_ref[k], theta_l[i, j, k],qt_[i, j, k])
                    if np.isnan(T_comp_thl[i,j,k]):
                        print('T_comp_thl is nan')
                        sys.exit()

                    if (ql_comp[i,j,k] - ql_[i, j, k]) > max_ql:
                        max_ql = (ql_comp[i,j,k]- ql_[i, j, k])
                    elif (ql_comp[i, j, k] - ql_[i, j, k]) < min_ql:
                        min_ql = (ql_comp[i, j, k] - ql_[i, j, k])

                    if ql_[i,j,k] > 0.0 and alpha[i,j,k] > 0.0:
                        if np.abs(T_comp[i,j,k] - T_[i, j, k]) > max_T_sat:
                            max_T_sat = np.abs(T_comp[i,j,k] - T_[i, j, k])
                    elif alpha[i,j,k] == 0.0:
                        if np.abs(T_comp[i,j,k] - T_[i, j, k]) > max_T_unsat:
                            max_T_unsat = np.abs(T_comp[i,j,k] - T_[i, j, k])
                            print('unsat max T: ', max_T_unsat)

                    if (ql_comp_thl[i,j,k] - ql_[i,j,k]) > max_ql_thl:
                        max_ql_thl = ql_comp_thl[i,j,k] - ql_[i,j,k]
                    elif (ql_comp_thl[i,j,k] - ql_[i,j,k]) < max_ql_thl:
                        min_ql_thl = ql_comp_thl[i,j,k] - ql_[i,j,k]

                    if ql_[i,j,k] > 0.0 and alpha[i,j,k] > 0:
                        print('sat: ', np.abs(T_comp_thl[i,j,k] - T_[i,j,k]))
                        if np.abs(T_comp_thl[i,j,k] - T_[i,j,k]) > max_T_sat_thl:
                            max_T_sat_thl = np.abs(T_comp_thl[i,j,k] - T_[i,j,k])
                    elif alpha[i,j,k] == 0:
                        print('unsat', np.abs(T_comp_thl[i,j,k]-T_[i,j,k]) )
                        if np.abs(T_comp_thl[i,j,k]-T_[i,j,k]) > max_T_unsat_thl:
                            max_T_unsat_thl = np.abs(T_comp_thl[i,j,k]-T_[i,j,k])


        print('')
        print('From Entropy:')
        print('max T sat: ', max_T_sat)             # max_T_sat = 0.096
        print('max T unsat: ', max_T_unsat)         # max_T_unsat = 0.05
        print('max ql:', max_ql)                    # max_ql = 4.4e-5
        print('min ql:', min_ql)                    # min_ql = -6.7e-5
        print('')
        print('From Thetali:')
        print('max T sat: ', max_T_sat_thl)         # max_T_sat = 77
        print('max T unsat: ', max_T_unsat_thl)     # max_T_unsat = 0.05
        print('max ql:', max_ql_thl)                # max_ql = 0.56 --> ql_comp >> ql_data
        print('min ql:', min_ql_thl)                # min_ql = -5e-6
        print('')


        plot_snapshots(ql_, ql_comp, alpha_, alpha, 'ql')
        plot_snapshots(T_, T_comp, alpha_, alpha, 'T')


    return




# ________________________________________________________________________________________
# ________________________________________________________________________________________
# ________________________________________________________________________________________
''' plot mean profiles for visual reference'''
def plot_profiles(fullpath_in_ref, z_ref):
    global path
    ''' profiles '''
    qt_prof = read_in_netcdf('qt_mean', 'profiles', fullpath_in_ref)
    qt_max = read_in_netcdf('qt_max', 'profiles', fullpath_in_ref)
    qt_min = read_in_netcdf('qt_min', 'profiles', fullpath_in_ref)
    ql_prof = read_in_netcdf('ql_mean', 'profiles', fullpath_in_ref)
    ql_max = read_in_netcdf('ql_max', 'profiles', fullpath_in_ref)
    ql_min = read_in_netcdf('ql_min', 'profiles', fullpath_in_ref)
    s_prof = read_in_netcdf('s_mean', 'profiles', fullpath_in_ref)
    s_max = read_in_netcdf('s_max', 'profiles', fullpath_in_ref)
    s_min = read_in_netcdf('s_min', 'profiles', fullpath_in_ref)
    time_prof = read_in_netcdf('t', 'timeseries', fullpath_in_ref)
    dt_prof = time_prof[2] - time_prof[1]
    n_prof = np.int(time[1] / dt_prof)
    print('time prof: ', n_prof, time_prof[n_prof], time[1])

    fig = plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.plot(1e2*ql_prof[n_prof,:], z_ref, label='<ql>*100')
    plt.plot(ql_max[n_prof, :], z_ref, 'k--', linewidth=1)
    plt.plot(ql_min[n_prof, :], z_ref, 'k--', linewidth=1)
    plt.plot([0,0],[0,z_ref[-1]])
    plt.legend()
    plt.xlabel('ql')
    plt.subplot(1, 3, 2)
    plt.plot(qt_prof[n_prof, :], z_ref, label='<qt>')
    plt.plot(qt_max[n_prof, :], z_ref, 'k--', linewidth=1, label='max: '+str(np.amax(qt_max[n_prof])))
    plt.plot(qt_min[n_prof, :], z_ref, 'k--', linewidth=1, label='min: '+str(np.amin(qt_min[n_prof])))
    plt.plot([0, 0], [0, z_ref[-1]])
    plt.legend()
    plt.xlabel('qt')
    plt.subplot(1, 3, 3)
    plt.plot(s_prof[n_prof, :], z_ref)
    plt.plot(s_max[n_prof,:], z_ref, 'k--', linewidth=1)
    plt.plot(s_min[n_prof, :], z_ref, 'k--', linewidth=1)
    plt.plot([0, 0], [0, z_ref[-1]])
    plt.xlabel('s')
    plt.xlim([np.amin(s_prof[n_prof,:]),np.amax(s_max[n_prof,:])])
    plt.savefig(os.path.join(path,'mean_profiles.png'))
    return

def plot_snapshots(field_data_, field_comp, alpha_, alpha, var_name):
    global path
    plt.figure(figsize=(12,15))
    plt.subplot(6, 3, 1)
    plt.imshow(field_data_[0:30, 0:30, 15],interpolation="nearest")
    plt.title(var_name+' (data), nz =' + str(dx[2] * 15))
    plt.colorbar()
    plt.subplot(6, 3, 4)
    plt.imshow(field_comp[0:30, 0:30, 15], interpolation="nearest")
    plt.title(var_name+ ' (comp)', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 7)
    plt.imshow(field_data_[0:30, 0:30, 15] - field_comp[0:30, 0:30, 15], interpolation="nearest")
    plt.title(var_name+' difference', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 10)
    plt.imshow(alpha_[0:30, 0:30, 15],interpolation="nearest")
    plt.colorbar()
    plt.title('alpha (data)')
    plt.subplot(6, 3, 13)
    plt.imshow(alpha[0:30, 0:30, 15],interpolation="nearest")
    plt.colorbar()
    plt.title('alpha comp')
    plt.subplot(6, 3, 16)
    plt.imshow(alpha[0:30, 0:30, 15]-alpha_[0:30, 0:30, 15],interpolation="nearest")
    plt.colorbar()
    plt.title('difference alpha: data - comp',fontsize=8)


    plt.subplot(6, 3, 2)
    plt.imshow(field_data_[0:30, 0:30, 18],interpolation="nearest")
    plt.colorbar()
    plt.title('data, nz =' + str(dx[2] * 18))
    plt.subplot(6, 3, 5)
    plt.imshow(field_comp[0:30, 0:30, 18], interpolation="nearest")
    plt.title('comp', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 8)
    plt.imshow(field_data_[0:30, 0:30, 18] - field_comp[0:30, 0:30, 18], interpolation="nearest")
    plt.title(var_name + ' difference', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 11)
    plt.imshow(alpha_[0:30, 0:30, 18],interpolation="nearest")
    plt.colorbar()
    plt.title('alpha (data)')
    plt.subplot(6, 3, 14)
    plt.imshow(alpha[0:30, 0:30, 18],interpolation="nearest")
    plt.colorbar()
    plt.title('alpha_comp')
    plt.subplot(6, 3, 17)
    plt.imshow(alpha[0:30, 0:30, 18] - alpha_[0:30, 0:30, 18], interpolation="nearest")
    plt.colorbar()
    plt.title('difference alpha: data - comp', fontsize=8)


    plt.subplot(6, 3, 3)
    plt.imshow(field_data_[0:30, 0:30, 20], interpolation="nearest")
    plt.colorbar()
    plt.title('data, nz =' + str(dx[2] * 20))
    plt.subplot(6, 3, 6)
    plt.imshow(field_comp[0:30, 0:30, 20], interpolation="nearest")
    plt.title('comp', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 9)
    plt.imshow(field_data_[0:30, 0:30, 20] - field_comp[0:30, 0:30, 20], interpolation="nearest")
    plt.title(var_name + ' difference', fontsize=8)
    plt.colorbar()
    plt.subplot(6, 3, 12)
    plt.imshow(alpha_[0:30, 0:30, 20], interpolation="nearest")
    plt.colorbar()
    plt.title('alpha (data)')
    plt.subplot(6, 3, 15)
    plt.imshow(alpha[0:30, 0:30, 20], interpolation="nearest")
    plt.colorbar()
    plt.title('alpha comp')
    plt.subplot(6, 3, 18)
    plt.imshow(alpha[0:30, 0:30, 20] - alpha_[0:30, 0:30, 20], interpolation="nearest")
    plt.colorbar()
    plt.title('difference alpha: data - comp', fontsize=8)



    plt.savefig(os.path.join(path,'snapshot_'+var_name+'.png'))



if __name__ == "__main__":
    main()