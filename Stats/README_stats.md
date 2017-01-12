**HOW TO USE:**

**(1) run EM_PDF_univar.py (or EM_PDF_bivar.py), which will output nc-files with PDF parameters
python EM_PDF path-to-fields-repository
(2)**



**(1a) PDF Fitting:**

`EM_PDF.py`: fit Gaussian mixed model (GMM) independently at every level and save parameters (means, covariances, weights) in nc-files
    using scikit package: sklearn.mixture.GaussianMixture

`EM_PDF_univar.py`: fit UNIVARIATE Gaussian mixed model (GMM) independently at every level and save parameters (means, covariances, weights) in nc-files
    using scikit package: sklearn.mixture.GaussianMixture

`EM_PDF_bivar.py`: fit BIVARIATE Gaussian mixed model (GMM) independently at every level and save parameters (means, covariances, weights) in nc-files
    using scikit package: sklearn.mixture.GaussianMixture

**(1b) Plotting PDF parameters:**

`EM_PDF_plot.py`: plot how PDF means evolve over vertical levels and time

`EM_PDF_bivar_plot.py`: plot how PDF means evolve over vertical levels and time for bivariate PDF

**(2) Autoregression:**

`EM_PDF_stochastic.py`: read in files with GMM parameters and feed to VAR model
- using package:  `statsmodels` (former `scikits.statsmodels`)
- _Installation:_ https://pypi.python.org/pypi/statsmodels
- _Documentation:_ http://www.statsmodels.org/stable/
- _IO:_ can load pickle, .dta-files; can be extended by using pandas.io (Excel, CSV, HDF5)





**INFORMATION**

**Statsmodels package** (former scikits.statsmodels package)
- linear regression
- time series analysis:
    - univariate time series analysis (AR, ARIMA)
    - vector autregressive models (VAR and structural VAR)
    - Datasets for examples and testing
- statistical tests

Fitting a model in statsmodels typically involves 3 easy steps:
(1) Use the model class to describe the model
(2) Fit the model using a class method (--> use dir(res))
(3) Inspect the results using a summary method

_**OLS:**_
yt = μ +a1*yt−1 +...+ak"yt−p +εt +b1*εt−1 +...+bq*εt−q
(1) model = sm.OLS(y,X)
(2) res = model.fit(n)  # n is the lag order
(3) print(res.summary())
other attributes:
res.plot()

_**Vector Autorgressive (VAR) models**_
--> model relationships among multiple time series
Y_t = A_1*Y_t−1 + ... + A_p*Y_t−p + ε_t.
Y_t: K-dimensional vector
A_i: square matrices
ε_t: typically assumed to be normally distributed and uncorrelated over time

- Lag order selection: how myna lags Y_t-i to include


**Datasets**
A dataset is a dictionary-like object that holds all the data and some metadata about the data.
This data is stored in the .data member, which is a n_samples, n_features array.




