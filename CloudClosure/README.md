



# Cloud Closure

**Description:**

`CloudClosure.py`, `CloudClosure_s.py`:
1. from LES fields compute bivariate Gaussian (one component) PDF f(th_l, q_t),
using likelihood maximization in sklearn.GaussianMixture
(`CloudClosure_s.py`: provides the option to use different 'temperature variables' for bivariate PDF)
2. Compare to Kernel Density Estimation
3. estimate q_l from this PDF (Monte Carlo calculation) and compare to the 'real' q_l from the LES data
(using the satuartion adjustment scheme from PyCLES)



**HOW TO USE:**

`CloudClosure.py`:

_input:_ reads in nc-files with 3D field output, computes (1) the bivariate Gaussian from maximum likelihood estimation and (2) the empirical PDF from Kernel Density Estimation (KDE)

_output:_ PDF parameters for Gaussian PDF




**Kernel Density Estimation**
The principle behind nearest _neighbor methods_ is to find a predefined number of training samples closest in distance to the new point, and predict the label from these.

The implemented Kernel Density algorithm uses the Ball Tree or KD Tree for efficient queries

- *Kernel Density Estimation*: `from sklearn.neighbors import KernelDensity`
`class sklearn.neighbors.KernelDensity(bandwidth=1.0, algorithm='auto', kernel='gaussian', metric='euclidean', atol=0, rtol=0, breadth_first=True, leaf_size=40, metric_params=None)`
- parameters:
    - kernel = 'gaussian' (‘gaussian’|’tophat’|’epanechnikov’|’exponential’|’linear’|’cosine’; default is 'gaussian')
    - bandwith = number (controlls "smoothness" of estimation)
- methods:
    - `KernelDensity.fit(data)`
    - `KernelDensity.fit(data).score_samples()`:  Evaluate the density model on the data.
    - `KernelDensity.fit(data).score(data)`:          total log probability
    - `KernelDensity.fit(data).sample(data)`: generate random variables from model




