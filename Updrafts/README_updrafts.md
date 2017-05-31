#Concept:
###1. read in
read in **3D output files from LES**
###2. PDF Model
1. **GMM**: compute double-Gaussian PDF (GMM algorithm) from 3D output file per level
2. **Labeling**: find labels to which cluster each data point belongs) for each data point
Attribute: `clf.predict(X[, y])`, `X.shape = (# data, # features)`
    * rearrange labels: returns array with 0/1 of size #data >> return into 2D array
    * sort PDF components, s.t. 0/1 always correspond to the same component (environment vs. updrafts)

Option: read in fields and PDF parameters from Cloud Closure and use clf.predict(field_data) for producing labels

###3. **Tracer Model**
read in Colleens output data (pickles-data with 0/1 encoding of updraft or environment)

###4. Comparison
Compare both clustering algorithms:
- __plotting:__ coloring scheme
blue = environment in both
green = updraft in both
red = differences (could be further distinguised into which algorithm classifies it as an updraft point)
- compute __conditional statistics__, e.g. cloud fraction and compare



# Algorithm:
### Updrafts.initialize

**1.** `Updrafts.initialize:`
- case_name
- ref.pressure, ref.height

**2.** `Updrafts.update:`

 **(a)**  read in files at time d (path from args.parse): qt_, s_, ql_, T_ = (nx,ny,nz)

**(b) PDF Model:**
- initialize CC, Lv-Lookup Table
- 3D Data collected in 2D array
- normalise Data
- GMM model --> PDF (clf)
- compute Labels array
- sort PDF: sort parameters and labels according to <qt>
- rearrange Data into 2D array

return: labels = (nx, ny, nz)

**(c) Tracer model:**

for type in list:
- read in 3D array = (nx, ny, nz)
- plotting(...)
