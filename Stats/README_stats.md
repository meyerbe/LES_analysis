
Statsmodels package (former scikits.statsmodels package)
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


OLS:
yt = μ +a1*yt−1 +...+ak"yt−p +εt +b1*εt−1 +...+bq*εt−q
(1) model = sm.OLS(y,X)
(2) res = model.fit(n)  # n is the lag order
(3) print(res.summary())
other attributes:
res.plot()


Vector Autorgressive (VAR) models

--> model relationships among multiple time series


Y_t = A_1*Y_t−1 + ... + A_p*Y_t−p + ε_t.
Y_t: K-dimensional vector
A_i: square matrices
ε_t: typically assumed to be normally distributed and uncorrelated over time

- Lag order selection: how myna lags Y_t-i to include

------------------------------------------------------------------------------------
------------------------------------------------------------------------------------
Datasets
A dataset is a dictionary-like object that holds all the data and some metadata about the data.
This data is stored in the .data member, which is a n_samples, n_features array.



EM_PDF.py: fit Gaussian mixed model (GMM) and save parameters (means, covariances, weights) in nc-files
    using scikit package: sklearn.mixture.GaussianMixture

EM_PDF_stochastic.py: read in files with GMM parameters and feed to VAR model
    using package:  statsmodels (former scikits.statsmodels)
                    Installation: https://pypi.python.org/pypi/statsmodels
                    Documentation: http://www.statsmodels.org/stable/
                    IO: can load pickle, .dta-files; can be extended by using pandas.io (Excel, CSV, HDF5)


