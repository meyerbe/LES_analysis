import os, sys
import argparse
import json as simplejson
import numpy as np
import pylab as plt

sys.path.append("..")
from io_read_in_files import read_in_netcdf

plt.rcParams['lines.linewidth'] = 2

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    global nx, ny, nz, dz
    nx = nml['grid']['nx']
    ny = nml['grid']['ny']
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']

    files = os.listdir(os.path.join(path, 'fields'))
    t = np.int(files[0][0:-3])
    path_fields = os.path.join(path, 'fields', files[0])
    ql = read_in_netcdf('ql', 'fields', path_fields)


    zql = []
    nql = []
    nql_all = []
    for z in range(nz):
        n = np.count_nonzero(ql[:,:,z])
        nql_all.append(n)
        if n > 0:
            zql.append(z)
            nql.append(n)

    plot_n(nql_all, path, case_name, t)

    return


def plot_n(nql_all, path, case_name, t):
    global nx, ny, nz, dz
    path_references = os.path.join(path, 'Stats.' + case_name + '.nc')
    zrange = read_in_netcdf('z', 'reference', path_references)

    plt.figure(figsize=(10,7))
    plt.plot(nql_all, zrange, label = 't='+str(t)+'s')
    plt.legend()
    plt.xlabel('# non-zero ql')
    plt.ylabel('height z [m]')
    plt.title('# non-zero ql ('+case_name+', nx*ny=1'+str(nx*ny)+')')

    plt.savefig(
        os.path.join(path, 'figs_stats', 'ql_number.png'))
    # plt.show()
    return


if __name__ == "__main__":
    main()