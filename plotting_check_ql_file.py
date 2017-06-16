import os, sys
import argparse
import json as simplejson
import numpy as np
import pylab as plt
import netCDF4 as nc

sys.path.append("..")
from io_read_in_files import read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 8


def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    parser.add_argument("--files")
    parser.add_argument("--files_cum")
    args = parser.parse_args()
    path = args.path
    global case_name
    case_name = args.casename
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dz, dt_stats
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']
    dt_stats = nml['stats_io']['frequency']
    print('Grid: dz='+str(dz))
    global max_height
    print('')

    files = os.listdir(os.path.join(path, 'fields'))
    print('All Files: ', files)
    t = np.int(files[0][0:-3])



    # path_fields = os.path.join(path, 'fields', str(t)+'.nc')
    # ql = read_in_netcdf('ql', 'fields', path_fields)
    print('path_ref', path_ref)
    # time_stats = read_in_netcdf('z', 'reference', path_ref)
    # print('time stats: ', time_stats[0:5], time_stats[-3:])
    time_stats = nc.Dataset(path_ref, 'r').groups['profiles'].variables['t'][:]
    print('time stats: ', time_stats[0:5], time_stats[-3:])
    zrange_stats1 = read_in_netcdf('z', 'reference', path_ref)
    zrange_stats = nc.Dataset(path_ref, 'r').groups['profiles'].variables['z_half'][:]
    zrange_stats2 = nc.Dataset(path_ref, 'r').groups['profiles'].variables['z'][:]

    # if str(files_[0])[-3:] == '.nc':
    if str(files[0])[-3:] == '.nc':
        path_z = os.path.join(path, 'fields', files[0])
    else:
        path_z = os.path.join(path, 'fields', str(files[0])+'.nc')
    print(path_z)
    try:
        zrange_field = nc.Dataset(path_z, 'r').groups['fields'].variables['z'][:]
        zrange = zrange_field
    except:
        zrange = zrange_stats
    print('')
    # print('zrange stats: ', zrange_stats1[0:5], zrange_stats1[20:23])
    print('zrange stats z_half:  ', zrange_stats[0:5], zrange_stats[20:23])
    if case_name == 'TRMM_LBA':
        zrange_stats_zp = nc.Dataset(path_ref, 'r').groups['reference'].variables['zp'][:]
        zrange_stats_zp_half = nc.Dataset(path_ref, 'r').groups['reference'].variables['zp_half'][:]
        print('zrange stats zp:      ', zrange_stats_zp[0:3], zrange_stats_zp[-2:])
        print('zrange stats zp_half: ', zrange_stats_zp_half[0:3], zrange_stats_zp_half[-2:])
        print('zrange field:         ', zrange_field[0:3], zrange_field[-2:])
        zrange_stats = zrange_stats

    day = np.int(24 * 3600)
    hour = np.int(3600)
    levels, files_, files_cum, time_prof = set_levels(case_name, files, zrange_stats, dz)
    print('Selected Files: ', files_)

    max_height = 120
    plot_mean_profile('thetali', time_prof, zrange_stats, max_height, path_ref, path, False, 4)
    plot_mean_profile('thetali', time_prof, zrange_stats, max_height, path_ref, path, True, 4)
    plot_mean_profile('cloud_fraction', time_prof, zrange, max_height, path_ref, path, True)
    # plot_mean_profile('fraction_core', time_prof, zrange, max_height, path_ref, path, True)
    # plot_mean_profile('fraction_cloud', time_prof, zrange, max_height, path_ref, path, True)
    plot_mean_profile('qt', time_prof, zrange, max_height, path_ref, path, True)
    plot_mean_profile('ql', time_prof, zrange, max_height, path_ref, path, True)
    # plot_ql_n_all(files_, zrange, path, prof=False)
    plot_mean('qt', files_, zrange, levels, path)
    plot_mean('ql', files_, zrange, levels, path)
    plot_mean('s', files_, zrange, levels, path)
    # try:
    #     plot_mean('thetali', files_, zrange, levels, path)
    #     plot_mean_var('thetali', files_, zrange, path)
    # except:
    #     print('thetali not in variables')
    plot_mean_levels('qt', files_cum, zrange, path)
    plot_mean_levels('ql', files_cum, zrange, path)

    plot_mean_var('ql', files_, zrange, path)
    plot_mean_var('qt', files_, zrange, path)
    plot_mean_var('s', files_, zrange, path)


    plot_mean_cumulated('ql', files_cum, zrange, levels, path)
    plot_mean_cumulated_BL('ql', files_cum, zrange, levels, path, 125)
    plot_mean_cumulated('s', files_cum, zrange, levels, path)
    plot_mean_cumulated_BL('ql', files_cum, zrange, levels, path, max_height)
    plot_mean_cumulated_BL('s', files_cum, zrange, levels, path, max_height)

    plot_max_var('ql', zrange, path)

    return



# ----------------------------------
def set_levels(case_name, files, zrange_stats, dz):
    global max_height
    day = np.int(24 * 3600)
    hour = np.int(3600)

    files_ = files
    files_cum = files
    time_prof = [0]

    if case_name[0:8] == 'ZGILS_S6':
        # files_ = ['1382400.nc']
        # files_ = [1317600, 1339200, 1360800, 1382400]  # ZGILS 6
        # files = files[0:25:2]
        files_ = files[0:44:10]
        files_cum = files[0:44:2]
        levels = dz*np.asarray([25, 25, 40, 50, 60, 65], dtype=np.int32)
    elif case_name[0:9] == 'ZGILS_S12':
        # ZGILS S12
        files_ = [0, 1 * day, 2 * day, 3 * day, 4 * day, 5 * day, 6 * day, 7 * day, 8 * day]
        files_cum = ['172800.nc', '259200.nc', '345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
        # files_cum = ['345600.nc', '432000.nc', '518400.nc', '604800.nc', '691200.nc']
        # levels = dz*np.asarray([35,40,45])
        levels = np.asarray([700,800,900,1000])
        time_prof = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5*hour]
    elif case_name == 'DYCOMS_RF01':
        # DYCOMS RF01 large
        # files_ = ['3600.nc']
        files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour]
        files_cum = ['10800.nc', '12600.nc', '14400.nc']
        time_prof = [1 * hour, 2 * hour, 3 * hour, 4 * hour]
        # DYCOMS RF01
        files_ = ['10800.nc']
        # files = ['10800.nc', '12600.nc', '14400.nc']
        levels = dz*np.asarray([140, 150, 160, 166, 180], dtype=np.int32)
    elif case_name == 'DYCOMS_RF02':
        files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
        files_ = [0, 1 * hour, 2 * hour]
        # files_ = ['7200.nc']
        # files_cum = ['10800.nc', '12600.nc', '14400.nc']
        files_cum = ['18000.nc', '19800.nc', '21600.nc']
        levels = dz*np.asarray([120, 140, 150, 160, 170, 200], dtype=np.int32)
    elif case_name == 'Bomex':
        ## Bomex 170314_weno7 (dz=40)
        # files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
        # # files_ = ['18000.nc', '19800.nc', '21600.nc']
        # files_cum = [5*hour, np.int(5.5*hour), 6*hour]
        # levels = [400, 500, 600, 720, 800, 900, 1000]
        ## Bomex (dz=20)
        files_ = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6 * hour]
        files_cum = [5 * hour, np.int(5.5 * hour), 6 * hour]
        levels = [500, 600, 720, 800, 1000, 1200, 1500, 2000]
        time_prof = [0, 1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 6*hour]
        ## Bomex test
        # levels = [400, 500, 600, 720, 800, 900, 1000]
        # files_ = ['21600.nc']
        # files_cum = files_

    elif case_name == 'TRMM_LBA':
        files_ = ['1003600.nc', '1007200.nc', '1010800.nc', '1014400.nc', '1018000.nc']
        # files_cum = ['1014400.nc', '1016200.nc', '1018000.nc']
        files_cum = ['1016200.nc', '1018000.nc', '1019800.nc']
        k_levels = np.asarray([10, 20, 30, 40, 50, 60, 68, 75, 85, 95, 105, 127], dtype=np.int32)
        levels = np.zeros(shape=k_levels.shape)
        for k in range(k_levels.shape[0]):
            levels[k] = zrange_stats[k_levels[k]]

        # levels_ = np.asarray([10, 20, 30, 40, 50, 60, 75, 85, 95, 105, 127], dtype=np.int32)
        # levels = np.asarray([1e3, 2e3, 3e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3, 12e3, 14e3, 16e3, 18e3])
        time_prof = [1 * hour, 2 * hour, 3 * hour, 4 * hour, 5 * hour, 5.5*hour]
        max_height = 130
        # Note:
        # - 'z' from fields is equivalent to 'z_half' from the Stats-File, group 'profiles'
        #           and 'zp_half' from Stats-file, gropu 'reference'
        # - Stats-file: 'z' from 'reference' is not identical with 'z' from 'profiles'

    elif case_name == 'Rico':
        files_ = [0, 2 * hour, 4 * hour, 6 * hour, 8 * hour, 10 * hour, 12 * hour,
                        14 * hour, 16 * hour, 18 * hour, 20 * hour, 22 * hour, 24 * hour]
        files_cum = [23 * hour, np.int(23.5 * hour), 24 * hour]

    print('levels', levels)
    return levels, files_, files_cum, time_prof


# ----------------------------------

def plot_mean_profile(var_name, time_range, zrange, max_height, path_ref, path, BL=False, location=1):
    print('-- plot mean from profile: '+ var_name + ' --')
    global dt_stats
    # path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    time_ = read_in_netcdf('t', 'timeseries', path_ref)
    if var_name == 'cloud_fraction' or var_name == 'fraction_core' or var_name == 'fraction_cloud':
        var = read_in_netcdf(var_name, 'profiles', path_ref)
    else:
        var_name = var_name + '_mean'
        var = read_in_netcdf(var_name, 'profiles', path_ref)


    print('')
    print('var shape: ', var.shape, 'time: ', time_.shape, time_range, ' zrange:', zrange.shape, 'dt_stats: ', dt_stats)
    print('')
    plt.figure(figsize=(9,6))
    cm1 = plt.cm.get_cmap('bone')
    # for t in range(time_.shape[0]):
    #     if time_[t]>=3600.0 and np.mod(time[t], 100*dt_stats) == 0.0:
    #         print('timetimetime', time_[t], 20*dt_stats)
    #         plt.plot(var[t,:], zrange, label='t='+str(time[t]))
    count_color = 0
    t_ini = 0
    count_t = 0
    for t in time_range:
        for t_ in range(t_ini, time_.shape[0]):
            if np.abs(time_[t_] - time_range[count_t]) < dt_stats:
                lab = set_tlabel(time_[t_])
                if BL:
                    plt.plot(var[t_, 0:max_height], zrange[0:max_height], color=cm1(np.double(count_color)/len(time_range)), label=lab)
                else:
                    plt.plot(var[t_, :], zrange, color=cm1(np.double(count_color)/len(time_range)), label=lab)
                t_ini = t_+1
                count_color += 1
                continue
        count_t += 1

    plt.legend(loc=location)
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean ' + var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    if BL:
        plt.savefig(
            os.path.join(path, 'figs_stats', var_name + '_fromprofile_BL.pdf'))
    else:
        plt.savefig(
            os.path.join(path, 'figs_stats', var_name + '_fromprofile.pdf'))
    # plt.show()
    # plt.close()

    return


def plot_mean_cumulated_BL(var_name, files_cum, zrange, levels, path, max_height):
    print('')
    print('-- plot mean cumulated BL: ' + var_name + ' --')
    global case_name
    global nz
    mean_all = np.zeros(nz)
    # path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    # time = read_in_netcdf('t', 'timeseries', path_references)

    print('files_cum', files_cum, len(files_cum))

    fig1 = plt.figure(figsize=(9, 6))
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    min = 9999.9
    max = -9999.9
    for t in files_cum:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.double(t)
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)
        mean_all += var_mean_field

        min_ = np.amin(var_mean_field[0:max_height])
        max_ = np.amax(var_mean_field[0:max_height])
        if min_ < min:
            min = min_
        if max_ > max:
            max = max_

        lab = set_tlabel(it)
        plt.plot(var_mean_field[0:max_height], zrange[0:max_height], '--', color=cm1(count_color / len(files_cum)), label='t=' + lab + ' (from field)')
        count_color += 1.0
    mean_all /= len(files_cum)
    plt.plot(mean_all[0:max_height], zrange[0:max_height], 'k', label='time mean')

    print('plotting levels', levels)
    for i in levels:
        if i < zrange[max_height]:
            # k = levels[i]
            # print(i, zrange[i], levels[i])
            # print('plotting levels: ', k)
            # k = zrange[i]
            # plt.plot([min,max],[k,k], linewidth=0.5, color='0.5', label=str(np.int(i))+'m (k=' + str(i) + ')' )i]
            plt.plot([min,max],[i,i], linewidth=0.5, color='0.5', label=str(np.int(i))+'m' )

    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean ' + var_name + ' (' + case_name + ', n*nx*ny=' + str(nx * ny * len(files_cum)) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_cum_BL.pdf'))
    # plt.show()
    plt.close()

    return

def plot_mean_cumulated(var_name, files_cum, zrange, levels, path):
    print('-- plot mean cumulated: ' + var_name + ' --')
    print(levels)
    global case_name
    global nz, max_height
    mean_all = np.zeros(nz)
    # path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    # time = read_in_netcdf('t', 'timeseries', path_references)

    print('')
    print('files_cum', files_cum, len(files_cum))

    fig1 = plt.figure(figsize=(9, 6))
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0


    for t in files_cum:
        if str(t)[-2:] == 'nc':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.double(t)
        var_field = read_in_netcdf(var_name, 'fields', path_fields)

        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)
        mean_all += var_mean_field

        lab = set_tlabel(it)
        plt.plot(var_mean_field[:], zrange, '--', color=cm1(count_color / len(files_cum)),
                     label='t=' + lab + ' (from field)')
        count_color += 1.0
    mini = np.amin(var_mean_field[:])
    maxi = np.amax(var_mean_field[:])
    for l in levels:
        plt.plot([mini, maxi], [l, l], color='0.75', linewidth=0.8, label=str(l) + 'm')

    mean_all /= len(files_cum)
    plt.plot(mean_all[:], zrange, 'k', label='time mean')
    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean ' + var_name + ' (' + case_name + ', n*nx*ny=' + str(nx * ny * len(files_cum)) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_cum.pdf'))
    # plt.show()
    plt.close()

    return


def plot_mean_var(var_name, files_, zrange, path):
    print('-- plot mean var: ' + var_name + ' --')
    global case_name
    global nx, ny, nz, dz, dt_stats
    day = 24 * 3600
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(14, 7))
    for t in files_:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = t

        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_fields = np.mean(np.mean(var, axis=0), axis=0)
        var2_mean = np.mean(np.mean(var*var, axis=0), axis=0)
        var_variance = var2_mean - var_mean_fields
        lab = set_tlabel(it)
        plt.subplot(1, 2, 1)
        plt.plot(var_mean_fields, zrange, color = cm2(count_color/len(files_)), label='t=' + lab)
        plt.subplot(1, 2, 2)
        plt.plot(var_variance, zrange, color = cm2(count_color/len(files_)), label='t=' + lab)
        count_color += 1.0

    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('mean ' + var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name+' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('Var['+var_name+']')
    plt.ylabel('height z [m]')
    plt.title('Var['+var_name+'] (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name+'_mean_var_fromfield.png'))
    # plt.show()
    plt.close()
    return


def plot_max_var(var_name, zrange, path):
    print('-- plot max var: ' + var_name + ' --')
    global case_name
    global nx, ny, nz, dz, dt_stats
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_max = read_in_netcdf(var_name + '_max', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)

    plt.figure(figsize=(14, 7))

    for it in range(var_max.shape[0]):
        tt = np.int((np.double(it) - time[0]) / np.double(dt_stats))
        plt.plot(var_max[tt, :], zrange, 'k--', label='t=' + str(tt * dt_stats) + 's')

    plt.legend()
    plt.xlabel('max' + var_name)
    plt.ylabel('height z [m]')
    plt.title('max' + var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.title('max[' + var_name + '] (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_max_fromprofile.png'))
    # plt.show()
    plt.close()
    return


def plot_ql_n_all(files_, zrange, path, prof=False):
    print('-- plot ql n all --')
    print(files_)
    global case_name
    global nx, ny, nz, dz, dt_stats
    day = 24 * 3600
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')

    count_color = 0.0
    plt.figure(figsize=(14,7))
    for t in files_:
        print('')
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
            if it >= 1000000:
                it = np.int(t[0:-3]) - 1000000
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.int(t)
            if it >= 1000000:
                it -= 1000000

        lab = set_tlabel(it)
        print('it: ', it, lab)

        ql = read_in_netcdf('ql', 'fields', path_fields)
        zql = []
        nql = []
        nql_all = []
        for z in range(nz):
            n = np.count_nonzero(ql[:, :, z])
            nql_all.append(n)
            if n > 0:
                zql.append(z)
                nql.append(n)

        ql_mean_fields = np.mean(np.mean(ql,axis=0),axis=0)
        ax = plt.subplot(1,2,1)
        plt.plot(nql_all, zrange, color = cm2(count_color/len(files_)), label='t=' + lab)
        if np.any(nql_all) > 0.0:
            try:
                ax.set_xscale('log')
            except:
                print('no log-scaling possible (nql_all=' + str(nql_all) + ')')
                pass
        plt.subplot(1, 2, 2)
        plt.plot(ql_mean_fields, zrange, color=cm2(count_color/len(files_)), label='t=' + lab)
        if prof:
            path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
            ql_mean = read_in_netcdf('ql_mean', 'profiles', path_references)
            time = read_in_netcdf('t', 'timeseries', path_references)
            tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
            plt.plot(ql_mean[tt,:], zrange,  'k--', label='t=' + lab)
        count_color += 1.0

    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('# non-zero ql')
    plt.ylabel('height z [m]')
    plt.title('# non-zero ql (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('mean ql')
    plt.ylabel('height z [m]')
    plt.title('mean ql (' + case_name + ', nx*ny=' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_number_fromfield.png'))
    # plt.show()
    plt.close()
    return



def plot_mean(var_name, files_, zrange, levels, path, prof = False):
    print('plotting mean')
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    # print(time)
    # cm1 = plt.cm.get_cmap('viridis')
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(9, 6))
    # mini = np.min([np.amin(var_mean[:, :]), np.amin(var_mean_field)])
    # maxi = np.max([np.amax(var_mean[:, :]), np.amax(var_mean_field)])
    mini = np.amin(var_mean[:, :])
    maxi = np.amax(var_mean[:, :])
    for l in levels:
        plt.plot([mini, maxi], [l, l], color='0.75', linewidth=0.8, label=str(l) + 'm')

    for t in files_:
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            it = np.double(t)
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)

        t_label = set_tlabel(it)
        plt.plot(var_mean_field[:], zrange, color=cm1(count_color/len(files_)), label=t_label)
        if prof:
            tt = np.int((np.double(it) - time[0]) / np.double(dt_stats))
            plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color/len(files_)), label='t=' + str(tt * dt_stats) + 's (from Stats)')
        count_color += 1.0

    plt.legend()
    plt.xlabel('mean '+var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield.png'))
    # plt.show()
    plt.close()
    return


def plot_mean_levels(var_name, files_, zrange, path, profile=False):
    print('')
    print('plotting mean levels')
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    var_mean = read_in_netcdf(var_name+'_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    # print(time)
    # cm1 = plt.cm.get_cmap('viridis')
    cm1 = plt.cm.get_cmap('bone')
    cm2 = plt.cm.get_cmap('winter')
    count_color = 0.0

    plt.figure(figsize=(9, 6))
    # print('files'), files_
    for t in files_:
        print('t', t)
        if str(t)[-1] == 'c':
            path_fields = os.path.join(path, 'fields', str(t))
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-3])
            else:
                it = np.int(t[0:-3])
        else:
            path_fields = os.path.join(path, 'fields', str(t) + '.nc')
            if case_name == 'TRMM_LBA':
                it = np.int(t[3:-1])
            else:
                it = np.double(t)
        # path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        var_field = read_in_netcdf(var_name, 'fields', path_fields)
        var_mean_field = np.mean(np.mean(var_field, axis=0), axis=0)

        tt = np.int((np.double(it)-time[0]) / np.double(dt_stats))
        # print('tt', tt, 'dt_stats', dt_stats, 'it', it, 'time[0]', time[0])
        try:
            mini = np.min([np.amin(var_mean[tt, :]), np.amin(var_mean_field)])
            maxi = np.max([np.amax(var_mean[tt, :]), np.amax(var_mean_field)])
        except:
            mini = np.amin(var_mean_field)
            maxi = np.amax(var_mean_field)

        if var_name == 'ql':
            location = 1
            ql = np.ndarray(shape=(0))
            z_ = np.ndarray(shape=(0))
            k_ = np.ndarray(shape=(0), dtype=np.int)
            for k in range(zrange.shape[0]):
                if var_mean_field[k] > 0.0:
                    ql = np.append(ql, var_mean_field[k])
                    z_ = np.append(z_, zrange[k])
                    k_ = np.append(k_, k)
            # if count_color == len(files_)-1:
            if count_color == 0:
                for l in z_:
                    if np.mod(l, 50) == 0:
                        plt.plot([mini, maxi], [l, l], color='0.5', linewidth=1.5, label=str(l) + 'm')
                    else:
                        plt.plot([mini, maxi], [l, l], color='0.5', linewidth=0.5)
            lab = set_tlabel(it)
            plt.plot(ql[:], z_, color=cm1(count_color / len(files_)),
                         label=lab)
            # plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color / len(files_)),
            #              label='t=' + str(tt * dt_stats) + 's (from Stats)')
        else:
            location = 3
            if count_color == 0.0:
                for l in zrange:
                    if np.mod(l,100) == 0:
                        plt.plot([mini, maxi], [l, l], color='0.5', linewidth=1.0, label=str(l)+'m')
                    elif np.mod(l,10*dz) == 0:
                        plt.plot([mini, maxi], [l, l], color='0.2', linewidth=0.2)
            lab = set_tlabel(it)
            plt.plot(var_mean_field[:], zrange, color=cm1(count_color/len(files_)), label=lab)
            if profile:
                plt.plot(var_mean[tt, :], zrange, '--', color=cm2(count_color/len(files_)), label='t=' + str(tt * dt_stats) + 's (from Stats)')
        count_color += 1.0


    plt.legend(loc=location, fontsize=6)
    plt.xlabel('mean '+var_name)
    plt.ylabel('height z [m]')
    plt.title('mean '+var_name + ' (' + case_name + ', nx*ny=' + str(nx * ny) + ', dz='+str(dz)+')')

    plt.savefig(
        os.path.join(path, 'figs_stats', var_name + '_mean_fromfield_levels.pdf'))
    # plt.show()
    plt.close()
    return


# _______
def set_tlabel(tt):
    if tt >= 1000000:
        tt -= 1000000
    if tt < 3600:
        lab = str(tt) + 's'
    elif tt < 3600 * 24:
        h = np.double(np.floor(tt / 1800)) / 2
        s = np.mod(tt, h*3600)
        if s > 0:
            lab = str(np.round(h, 1)) + 'h ' + str(np.round(s, 0)) + 's'
        else:
            lab = str(np.round(h, 1)) + 'h'
        print('set_tlabel: ', 'tt', tt, h, s, lab)
    else:
        d = tt / day
        h = np.mod(tt, day) / 3600
        s = np.mod(np.mod(tt, day), 3600)
        if s > 0:
            lab = str(np.round(d, 0)) + 'days ' + str(np.round(h, 0)) + 'h ' + str(np.round(s, 0)) + 's'
        elif h > 0:
            lab = str(np.round(d, 0)) + 'days ' + str(np.round(h, 0)) + 'h '
        else:
            lab = str(np.round(d, 0)) + 'days'
    return lab

if __name__ == "__main__":
    main()