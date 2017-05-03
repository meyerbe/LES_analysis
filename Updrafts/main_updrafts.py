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
    case_name = args.casename

    nml = simplejson.loads(open(os.path.join(path, case_name + '.in')).read())
    dz = nml['grid']['dz']

    files = os.listdir(os.path.join(path, 'fields'))
    # ncomp_range = [1, 2, 3, 4, 5, 6, 7, 8]
    dz_range = 1
    ncomp_range = [2]

    # Bomex test
    files = ['21600.nc']
    # krange = np.asarray([10, 17, 20, 25, 50])
    krange = np.asarray([20, 25], dtype=np.int32)
    # krange = np.asarray([20, 25])
    # krange = np.asarray([18,30,38])
    print('...', type(krange), type(krange[0]))

    N = len(files)
    print('Found the following directories', files, N)
    print('zrange: ' + str(krange * dz))
    print('dz: ' + str(dz))
    print('ncomp: ' + str(ncomp_range))
    print('')


    # Updrafts.CloudClosure()
    Up = Updrafts.Updrafts()
    Up.initialize(krange, path, case_name)
    Up.predict_pdf(files, path, ncomp_range, dz_range, krange, nml)
    return




if __name__ == '__main__':
    main()