# LES_analysis

parameters: LES parameters

# IO

read_in_netcdf4: example for reading in NetCDF4 files




# Offline Calculations

budget_TKE: offline TKE budget computation from LES output (old LES; hdf5 output)
budget_Th: offline potential temperature budget (potential energy)

thermodynamics_1.py: offline thermodynamic computations (Clausius Clapeyron by Magnus Formula; saturation adjustment) as in LES

compatibility: check compatibility of scalar tracer convergence according to Chr. Schär



# Fluctuations and Correlations

EddyField_output: compute Eddy Fields (fluctuations) from LES field output and save as new NC-files

Correlations: compute Correlations from LES output fields and output as NC-files



# Visualization
Vis: plot from pickle

Vis_fields: 
    - plot from LES output fields and profiles
    - compute Correlations
    - compute 

