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


'''
Test Types:
input data (p_ref, s_[i,j,k], qt_[i,j,k] are float64
s_dry, s_vap, s_cond: all float64
'''

sys.path.append("..")
sys.path.append("../Thermo/")
sys.path.append("../CloudClosure/")
# pyximport CC_thermodynamics
import CC_thermodynamics
import CC_thermodynamics_c
from io_read_in_files import read_in_nml, read_in_netcdf
from CC_thermodynamics import sat_adj_fromentropy, sat_adj_fromthetali
from CC_thermodynamics_c import sat_adj_fromentropy_c
from thermodynamic_functions import theta_li

def main():
    print('test CC cython')
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
        fullpath_in_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
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
    # ______________________
    ''' zrange '''
    global zrange
    zrange = [15,18,20]
    # ______________________
    # ______________________
    # Initialize Lookup Table for Clausius Clapeyron
    # CC = CC_thermodynamics.ClausiusClapeyron()
    CC = CC_thermodynamics_c.ClausiusClapeyron_c()
    # CC = CC_thermodynamics.ClausiusClapeyron_pycles()
    CC.initialize()

    LH = CC_thermodynamics_c.LatentHeat()
    # if microphysics == 'dry':
    #     LH.Lambda_fp = lambda_constant
    #     LH.L_fp = latent_heat_constant
    # elif microphysics == 'sa':
    #     LH.Lambda_fp = lambda_constant
    #     LH.L_fp = latent_heat_variable
    # L_fp = LH.L_fp
    # Lambda_fp = LH.Lambda_fp$

    microphysics = 'sa'
    # ______________________
    files_ = [files[1]]
    # files_ = files
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
        print('')
        print('types: ', type(p_ref[0]), type(s_[0,0,0]), type(qt_[0,0,0]))     #
        print('shapes: ', s_.shape)
        print('zrange: ', zrange)
        print('')
        # for k in zrange:
        for k in [15]:
            # for i in range(nx[0]):
            #     for j in range(nx[1]):
            for i in range(10):
                for j in range(10):
                    print('ijk', i, j, k, p_ref[k], s_[i,j,k], qt_[i,j,k])
                    # T_comp, ql_comp = sat_adj(p, s[i,j], qt[i,j])
                    # T_comp[i, j, k], ql_comp[i, j, k], alpha[i, j, k] = sat_adj_fromentropy(p_ref[k], s_[i, j, k],qt_[i, j, k])
                    # T_comp[i, j, k], ql_comp[i, j, k], alpha[i, j, k] = sat_adj_fromentropy_double(p_ref[k], s_[i, j, k],
                    #                                                                         qt_[i, j, k])
                    T_comp[i, j, k], ql_comp[i, j, k], alpha[i, j, k] = sat_adj_fromentropy_c(p_ref[k], s_[i, j, k],
                                                                                            qt_[i, j, k], microphysics)

                    # theta_l[i, j, k] = theta_li(p_ref[k], T_[i, j, k], qt_[i, j, k], ql_[i, j, k], 0)
                    # T_comp_thl[i, j, k], ql_comp_thl[i, j, k], alpha_thl[i, j, k] = sat_adj_fromthetali(p_ref[k], theta_l[i, j, k],qt_[i, j, k])






        #             if np.isnan(T_comp_thl[i,j,k]):
        #                 print('T_comp_thl is nan')
        #                 sys.exit()
        #
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

    #                 if (ql_comp_thl[i,j,k] - ql_[i,j,k]) > max_ql_thl:
    #                     max_ql_thl = ql_comp_thl[i,j,k] - ql_[i,j,k]
    #                 elif (ql_comp_thl[i,j,k] - ql_[i,j,k]) < min_ql_thl:
    #                     min_ql_thl = ql_comp_thl[i,j,k] - ql_[i,j,k]
    #
    #                 if ql_[i,j,k] > 0.0 and alpha[i,j,k] > 0:
    #                     print('sat: ', np.abs(T_comp_thl[i,j,k] - T_[i,j,k]))
    #                     if np.abs(T_comp_thl[i,j,k] - T_[i,j,k]) > max_T_sat_thl:
    #                         max_T_sat_thl = np.abs(T_comp_thl[i,j,k] - T_[i,j,k])
    #                 elif alpha[i,j,k] == 0:
    #                     print('unsat', np.abs(T_comp_thl[i,j,k]-T_[i,j,k]) )
    #                     if np.abs(T_comp_thl[i,j,k]-T_[i,j,k]) > max_T_unsat_thl:
    #                         max_T_unsat_thl = np.abs(T_comp_thl[i,j,k]-T_[i,j,k])


        print('')
        print('From Entropy:')
        print('max T sat: ', max_T_sat)             # max_T_sat = 0.096
        print('max T unsat: ', max_T_unsat)         # max_T_unsat = 0.05
        print('max ql:', max_ql)                    # max_ql = 4.4e-5
        print('min ql:', min_ql)                    # min_ql = -6.7e-5
        print('')
    #     print('From Thetali:')
    #     print('max T sat: ', max_T_sat_thl)         # max_T_sat = 0.12
    #     print('max T unsat: ', max_T_unsat_thl)     # max_T_unsat = 0.013
    #     print('max ql:', max_ql_thl)                # max_ql = 3.31e-5
    #     print('min ql:', min_ql_thl)                # min_ql = 0.0
    #     print('')
    #
    #
    # #     plot_snapshots(ql_, ql_comp, alpha_, alpha, 'ql')
    # #     plot_snapshots(T_, T_comp, alpha_, alpha, 'T')
    # #     plot_snapshots(theta_l, theta_l, alpha_, alpha, 'thetali')


    return


if __name__ == "__main__":
    main()