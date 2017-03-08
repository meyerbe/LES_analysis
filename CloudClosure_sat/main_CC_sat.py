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
    path_ref = os.path.join(path, 'Stats.Bomex.nc')

    # (0) Namelist File
    nml = simplejson.loads(open(os.path.join(path, 'Bomex.in')).read())
    nz = nml['grid']['nz']
    dz = nml['grid']['dz']

    # path = '../test_ZGILS6/'
    # path_ref = os.path.join(path, 'Stats.ZGILS_S6_1xCO2_SST_FixSub.nc')
    # do_everything(path, path_ref)
    # path = '../test/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)

    # path = '../test_bomex_n1024/'
    # path_ref = os.path.join(path, 'Stats.Bomex.nc')
    # do_everything(path, path_ref)
    # do_everything_with_pycles(path, path_ref)

    ClCl = CloudClosure.CloudClosure()
    ClCl.initialize(path, path_ref)
    ClCl.verification_CC(path, path_ref)

    print('--- PDF Prediction ---')
    print('')
    ncomp = 2
    zrange = np.arange(6, 20, 8)
    print('zrange', zrange*dz)
    ClCl.predict_pdf(path, path_ref, ncomp)

    return



if __name__ == "__main__":
    main()