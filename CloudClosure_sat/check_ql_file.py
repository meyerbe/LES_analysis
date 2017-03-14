import os, sys
import argparse
import json as simplejson
import numpy as np
import pylab as plt

sys.path.append("..")
from io_read_in_files import read_in_netcdf

label_size = 8
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['lines.linewidth'] = 2


def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
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

    files = os.listdir(os.path.join(path, 'fields'))
    print(files)
    t = np.int(files[0][0:-3])

    # files_ = [14400, 21600]
    # t = files_[0]
    files_ = [files[0][0:-3]]

    path_fields = os.path.join(path, 'fields', str(t)+'.nc')
    ql = read_in_netcdf('ql', 'fields', path_fields)
    zrange = read_in_netcdf('z', 'reference', path_ref)

    plot_ql_n_all(files_, zrange, path, prof=False)
    # plot_mean(files_, zrange, path)
    plot_qt(files_, zrange, path)
    return


def plot_qt(files_, zrange, path):
    global case_name
    global nx, ny, nz, dz, dt_stats

    plt.figure(figsize=(14, 7))
    for t in files_:
        path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        qt = read_in_netcdf('qt', 'fields', path_fields)
        qt_mean_fields = np.mean(np.mean(qt, axis=0), axis=0)
        qt_qt_mean = np.mean(np.mean(qt*qt, axis=0), axis=0)
        qt_var = qt_qt_mean - qt_mean_fields
        plt.subplot(1, 2, 1)
        plt.plot(qt_mean_fields, zrange, label='t=' + str(t) + 's')
        plt.subplot(1, 2, 2)
        plt.plot(qt_var, zrange, label='t=' + str(t) + 's')


    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('mean qt')
    plt.ylabel('height z [m]')
    plt.title('mean qt (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('Var[qt]')
    plt.ylabel('height z [m]')
    plt.title('Var[qt] (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.savefig(
        os.path.join(path, 'figs_stats', 'qt_mean_fromfield.png'))
    # plt.show()
    plt.close()
    return



def plot_ql_n_all(files_, zrange, path, prof=True):
    global case_name
    global nx, ny, nz, dz, dt_stats

    plt.figure(figsize=(14,7))
    for t in files_:
        path_fields = os.path.join(path, 'fields', str(t) + '.nc')
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
        plt.plot(nql_all, zrange, label='t=' + str(t) + 's')
        ax.set_xscale('log')
        plt.subplot(1, 2, 2)
        plt.plot(ql_mean_fields, zrange, label='t=' + str(t) + 's')
        if prof:
            path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
            ql_mean = read_in_netcdf('ql_mean', 'profiles', path_references)
            it = np.int(t / dt_stats)
            plt.plot(ql_mean[it,:], zrange, 'k--', label='t=' + str(t) + 's')

    plt.subplot(1, 2, 1)
    plt.legend()
    plt.xlabel('# non-zero ql')
    plt.ylabel('height z [m]')
    plt.title('# non-zero ql (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')
    plt.subplot(1, 2, 2)
    plt.legend()
    plt.xlabel('mean ql')
    plt.ylabel('height z [m]')
    plt.title('mean ql (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_number_fromfield.png'))
    # plt.show()
    plt.close()
    return



def plot_mean(files_, zrange, path):
    global case_name
    global nx, ny, nz, dz, dt_stats
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    ql_mean = read_in_netcdf('ql_mean', 'profiles', path_references)
    time = read_in_netcdf('t', 'timeseries', path_references)
    print(time)

    plt.figure(figsize=(10, 7))
    print('files'), files_
    for t in files_:
        path_fields = os.path.join(path, 'fields', str(t) + '.nc')
        ql = read_in_netcdf('ql', 'fields', path_fields)
        ql_mean_fields = np.mean(np.mean(ql, axis=0), axis=0)

        it = np.int(np.double(t) / np.double(dt_stats))
        print('it', it, 'dt_stats', dt_stats)
        # print('type t', type(t), t, type(dt_stats), dt_stats, type(it))
        # print('dt stats: ', t, dt_stats, it, ql_mean.shape)
        plt.plot(ql_mean[it, :], zrange, label='t=' + str(it * dt_stats) + 's (from Stats)')
        plt.plot(ql_mean_fields[:], zrange, label='t=' + str(it * dt_stats) + 's (from field)')
    plt.legend()
    plt.xlabel('mean ql')
    plt.ylabel('height z [m]')
    plt.title('mean ql (' + case_name + ', nx*ny=1' + str(nx * ny) + ')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_mean_fromfield.png'))
    # plt.show()
    plt.close()
    return



if __name__ == "__main__":
    main()