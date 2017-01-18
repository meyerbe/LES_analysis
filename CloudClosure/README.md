**HOW TO USE:**



# Cloud Closure
`CloudClosure.py`:
The principle behind nearest _neighbor methods_ is to find a predefined number of training samples closest in distance to the new point, and predict the label from these.

The implemented Kernel Density algorithm uses the Ball Tree or KD Tree for efficient queries

**Kernel Density Estimation**
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




