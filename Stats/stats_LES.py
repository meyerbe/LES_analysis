import numpy as np
import statsmodels.api as sm

import statsmodels.tsa.api

def test():
    # test data
    mdata = sm.datasets.macrodata.load().data
    mdata = mdata[['realgdp','realcons','realinv']]
    names = mdata.dtype.names
    data = mdata.view((float,3))
    data = np.diff(np.log(data), axis=0)


    model = statsmodels.tsa.api.VAR(data)
    model2 = statsmodels.tsa.vector_ar.var_model.VAR(data)
    # model = VAR(data, names=names)
    res = model.fit(2)
    print res.summary()

# real data
def var_modeling(phi,pmax):
    print('')
    print('----- var modeling -----')
    print phi.shape
    model = statsmodels.tsa.vector_ar.var_model.VAR(phi)
    res = model.fit(pmax)
    print res.summary()
    print('----------')
    return