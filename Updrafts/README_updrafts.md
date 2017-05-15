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


