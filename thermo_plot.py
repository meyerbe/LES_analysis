''' Auxiliary thermodynamic functions from LES '''
import pylab as plt
import numpy as np
from parameters import *
from math import *


from thermo_aux import pv_c



''' Vapor pressure '''
qv = np.linspace(0,0.5,11)
qt_ = np.linspace(0.0,0.8,5)
#qt = 0.3
p0 = 1e5
plt.figure()
for qt in qt_:
    plt.plot(qv,pv_c(p0,qt,qv),label='qt='+str(qt))
plt.legend()
plt.xlabel(r'$q_v$')
plt.ylabel('vapor pressure '+r'$p_v$')
plt.title('vapor pressure')
plt.savefig('figures/vapor_pressure.png')

''' Liquid potential temperature '''
