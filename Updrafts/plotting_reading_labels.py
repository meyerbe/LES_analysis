import os
import numpy as np
import pylab as plt
import argparse
import netCDF4 as nc

label_size = 10
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['lines.linewidth'] = 0.8




def main():
    parser = argparse.ArgumentParser(prog='PyCLES')
    parser.add_argument("--path")
    parser.add_argument("--casename")
    args = parser.parse_args()
    path = '/Volumes/Data/ClimatePhysics/LES/updrafts_colleen/'
    case_name = 'Bomex'
    if args.path:
        path = args.path
    if args.casename:
        case_name = args.casename

    time = 21600
    path_pdf = os.path.join(path, 'Updrafts')
    labels_pdf = read_in_pdf_labels(time, path_pdf)
    type = 'Coherent'
    path_colleen = os.path.join(path, 'tracer_fields')
    # labels_tracers = read_in_updrafts_colleen(type, time, path_colleen)
    # print('labels_tr', labels_tracers.shape)
    return


def read_in_pdf_labels(t, path_):
    print('--- read in PDF Labels ---')
    filename = 'Labeling_t'+str(t)+'.nc'
    # filename = 'Labeling_t' + '.nc'
    path = os.path.join(path_, filename)
    print(path)
    print('')

    root = nc.Dataset(path, 'r')
    # labels = root.group['fields'].variables['labels'][:,:,:]
    labels = root.group['fields']


    return 0



def read_in_updrafts_colleen(type, t, path_):
    print('')
    print('--- Updraft Colleen: read in ---')
    import pickle
    print('')


    path = path_
    files = os.listdir(path)
    print(path_)
    print('time: ', t)
    print('')

    if type == 'Cloud':
        root = 'Bomex_Cloud_'
    elif type == 'Coherent':
        root = 'Bomex_Coherent_'
    elif type == 'Couvreux':
        root = 'Bomex_Couvreux_'
    elif type == 'Core':
        root = 'Bomex_Core_'
    print(root)

    # print('')
    # path = os.path.join(path_, root + 'updraft.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())
    #
    # print('')
    # path = os.path.join(path_, root + 'environment.pkl')
    # data = pickle.load(open(path))
    # print(path + ': ', data.keys())

    print('')
    path = os.path.join(path_, root + 'time_'+ str(t) + '_Grid.pkl')
    print(path)
    print('')
    labels = pickle.load(open(path))
    # print(type(labels))
    # print(path + ': ', labels.shape())

    return labels




if __name__ == '__main__':
    main()
