import os
import argparse
import json as simplejson
import numpy as np

import Updrafts
# from Updrafts import Updrafts

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    path = args.path
    # path_tr = '/Volumes/Data/ClimatePhysics/LES/updrafts_colleen/'
    path_tr = os.path.join(path, 'tracer_fields')

    case_name = args.casename
    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    dz = nml['grid']['dz']
    nz = nml['grid']['nz']
    print('dz = '+str(dz))

    files_3d = os.listdir(os.path.join(path, 'fields'))
    # dz_range = 1
    ncomp_range = [2]

    # Bomex test
    files_3d = ['21600.nc']
    # krange = np.asarray([20, 30, 40, 50, 60, 70, 80, 90, 100], dtype=np.int32)
    # krange = np.arange(0, 120, 1, dtype=np.int32)
    krange = np.arange(0, nz, 1, dtype=np.int32)
    print('krange: ', type(krange), type(krange[0]))
    # krange = np.asarray([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100, 125], dtype=np.int32)

    N = len(files_3d)
    print('Found the following directories', files_3d, N)
    print('zrange: ' + str(krange * dz))
    print('dz: ' + str(dz))
    print('ncomp: ' + str(ncomp_range))
    print('')




    # Updrafts.CloudClosure()
    Up = Updrafts.Updrafts()
    Up.initialize(krange, path, case_name)
    # Up.predict_pdf(files_3d, path, ncomp_range, dz_range, krange, nml)
    # Up.update(files_3d, ncomp_range, dz_range, krange, nml, path, path_tr)
    Up.update(files_3d, ncomp_range, krange, nml, path, path_tr)

    return




if __name__ == '__main__':
    main()