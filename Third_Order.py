''' read in and plot 3rd order terms from LES output '''

__author__ = 'bettinameyer'

from h5py import File
import numpy as np
import matplotlib.cm as cm
import pylab as plt
#from Namelist import Namelist
#nml = Namelist()

# _________________________________________________________________________________
def read_in_hdf5(variable_name, group_name, fullpath_in):
    f = File(fullpath_in)

    #Get access to the profiles group
    profiles_group = f[group_name]
    #Get access to the variable dataset
    variable_dataset = profiles_group[variable_name]
    #Get the current shape of the dataset
    variable_dataset_shape = variable_dataset.shape

    variable = np.ndarray(shape = variable_dataset_shape)
    for t in range(variable_dataset_shape[0]):
        if group_name == "profiles":
            variable[t,:] = variable_dataset[t, :]
        if group_name == "correlations":
            variable[t,:] = variable_dataset[t, :]
        if group_name == "timeseries":
            variable[t] = variable_dataset[t]

    f.close()
    return variable




# _________________________________________________________________________________
# _________________________________________________________________________________

T = '7200'
tt = 2
dt_profile = 3600


data_full = 'dcbl_full.statistics.hdf5'
data_ql = 'dcbl_ql.statistics.hdf5'

path_full = 'dcbl/150918_2nd/full_SGSfull/'
path_ql = 'dcbl/150918_2nd/ql_SGSfull/'
path_out = 'dcbl/150918_2nd/'
print(path_full)
print(path_ql)

time_full = read_in_hdf5('time', "timeseries", path_full + data_full)
time_ql = read_in_hdf5('time', "timeseries", path_ql + data_ql)
print('time full',time_full[0:5]/3600, 't_max: ', time_full[-1]/3600, np.shape(time_full))
print('time ql',time_ql[0:5]/3600, 't_max: ', time_ql[-1]/3600, np.shape(time_ql))


wss_ql = read_in_hdf5('wss',"profiles",path_ql+data_ql)
wws_ql = read_in_hdf5('wws',"profiles",path_ql+data_ql)
wthth_ql = read_in_hdf5('wthth',"profiles",path_ql+data_ql)
wwth_ql = read_in_hdf5('wwth',"profiles",path_ql+data_ql)
www_ql = read_in_hdf5('w_resolved_third_central',"profiles",path_ql+data_ql)
wwu_ql = read_in_hdf5('wwu',"profiles",path_ql+data_ql)
#wwv_ql = read_in_hdf5('wwv',"profiles",path_ql+data_ql)


#nx = wss_ql.shape[0]
#ny = wss_ql.shape[1]
nz = wss_ql.shape[1]
##
nz = nz/2

dz = 25
zvector = dz*np.linspace(0,nz-1,nz)



plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
for t in [0,3,6,9]:
    plt.plot(wss_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title('wss')
plt.ylabel('height (dz = ' + np.str(dz) + '), [km]')
plt.subplot(1,2,2)
for t in [0,3,6,9]:
    plt.plot(wws_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title('wws')
plt.savefig('entropy_3rd.png')
plt.savefig(path_out + 'entropy_3rd.png')
plt.close()


plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
for t in [0,3,6,9]:
    plt.plot(wthth_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title(r'w$\theta\theta$')
plt.ylabel('height (dz = ' + np.str(dz) + '), [km]')
plt.subplot(1,2,2)
for t in [0,3,6,9]:
    plt.plot(wws_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title(r'ww$\theta$')
plt.savefig('pottemp_3rd.png')
plt.savefig(path_out + 'pottemp_3rd.png')
plt.close()

plt.figure(1,figsize=(20,10))
plt.subplot(1,2,1)
for t in [0,3,6,9]:
    plt.plot(www_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title(r'www')
plt.ylabel('height (dz = ' + np.str(dz) + '), [km]')
plt.subplot(1,2,2)
for t in [0,3,6,9]:
    plt.plot(wwu_ql[t,0:nz],zvector[0:nz],linewidth=2,label = 't = ' + np.str(t))
plt.legend()
plt.title(r'wwu')
plt.savefig('mom_3rd.png')
plt.savefig(path_out + 'mom_3rd.png')
plt.close()


