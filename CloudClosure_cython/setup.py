from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import sys
import platform
import subprocess as sp
import os.path



# Now get include paths from relevant python modules
# include_path = [mpi4py.get_include()]

include_path = [np.get_include()]
include_path += ['./Csrc']
# include_path += ['../Thermo/']

if sys.platform == 'darwin':
    #Compile flags for MacOSX
    library_dirs = []
    libraries = []
    extensions = []
    extra_compile_args = []
    extra_compile_args += ['-O3', '-march=native', '-Wno-unused', '-Wno-#warnings','-fPIC']
    # extra_objects=['./RRTMG/rrtmg_build/rrtmg_combined.o']
    extra_objects = []
    netcdf_include = '/opt/local/include'
    netcdf_lib = '/opt/local/lib'
    f_compiler = 'gfortran'
elif 'euler' in platform.node():
    #Compile flags for euler @ ETHZ
    library_dirs = ['/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib/']
    libraries = []
    libraries.append('mpi')
    libraries.append('gfortran')
    extensions = []
    extra_compile_args=[]
    extra_compile_args+=['-std=c99', '-O3', '-march=native', '-Wno-unused',
                         '-Wno-#warnings', '-Wno-maybe-uninitialized', '-Wno-cpp', '-Wno-array-bounds','-fPIC']
    # extra_objects=['./RRTMG/rrtmg_build/rrtmg_combined.o']
    extra_objects = []
    netcdf_include = '/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/include'
    netcdf_lib = '/cluster/apps/netcdf/4.3.1/x86_64/gcc_4.8.2/openmpi_1.6.5/lib'
    f_compiler = 'gfortran'

else:
    print('Unknown system platform: ' + sys.platform  + 'or unknown system name: ' + platform.node())
    sys.exit()

_ext = Extension('CC_thermodynamics_c', ['CC_thermodynamics_c.pyx'], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)

# _ext = Extension('PyCLES_Thermodynamics', ['PyCLES_Thermodynamics.pyx'], include_dirs=include_path,
#                  extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
#                  runtime_library_dirs=library_dirs)
# extensions.append(_ext)


_ext = Extension('thermodynamic_functions_c', ['thermodynamic_functions_c.pyx'], include_dirs=include_path,
                 extra_compile_args=extra_compile_args, libraries=libraries, library_dirs=library_dirs,
                 runtime_library_dirs=library_dirs)
extensions.append(_ext)


setup(
    ext_modules=cythonize(extensions, verbose=1, include_path=include_path)
)
