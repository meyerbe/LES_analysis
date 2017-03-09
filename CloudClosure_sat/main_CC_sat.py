import sys, os
import argparse
import json as simplejson
import numpy as np

import CloudClosure

# from CloudClosure import do_everything_with_pycles
# from CloudClosure import verification_CC

def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("path")
    parser.add_argument("casename")
    args = parser.parse_args()
    path = args.path
    case_name = args.casename
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')

    # (0) Namelist File
    # path = '../test_ZGILS6/'
    # path_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
    # case_name = 'ZGILS_S6_1xCO2_SST_FixSub'

    # path = '../test/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)

    path = '../test_bomex_n1024/'
    path_ref = os.path.join(path, 'Stats.Bomex.nc')
    case_name = 'Bomex'

    path = '../test_bomex/'
    case_name = 'Bomex'
    path_ref = os.path.join(path, 'Stats.'+case_name+'.nc')


    nml = simplejson.loads(open(os.path.join(path, case_name+'.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']

    ClCl = CloudClosure.CloudClosure()
    ClCl.initialize(path, path_ref, case_name)
    # ClCl.verification_CC(path, path_ref)

    ncomp = 2
    krange = np.arange(6, 20, 8)
    print('zrange', krange*dz)
    print('')
    ClCl.predict_pdf(path, path_ref, ncomp, krange, nml)

    return



if __name__ == "__main__":
    main()