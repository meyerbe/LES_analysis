import sys, os
import argparse

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
    # initialize_CC()
    ClCl.verification_CC(path, path_ref)

    return



if __name__ == "__main__":
    main()