import matlab.engine
from arfit import arfit
import numpy as np


eng = matlab.engine.start_matlab()

# tf = eng.isprime(37)
# print(tf)

# eng.test_m(nargout=0)
# eng.edit('test_m',nargout=0)
ret = eng.test_m(1.,5.)
print(ret)


# Creates a list containing 5 lists, each of 8 items, all set to 0
w, h = 3, 10
Matrix = [[0 for x in range(w)] for y in range(h)]
print Matrix
# # v = np.ones((10,2))
# v = [[1,1,1,1,1,1,1,1,1,1,1,1,1],[1,1,1,1,1,1,1,1,1,1,1,1]]
pmin = 1
pmax = 2
# # arfit(v,pmin,pmax)
eng.arfit(Matrix,pmin,pmax)


