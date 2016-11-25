???????????
- entropies.py: sv_c() -- correct to set it to zero if pv=0?

- pv_star_1 > p0 --> qv_star_1<0 (bcs. in calculation qv_star_c: pd<0) --> qt>qv_star_1 (saturated even for qt=0) but pv=0 
    --> changed formulation for (saturated && pv==0) to s_1 = s (initial entropy)
- if continuing saturated thermodynamics with qt=0, pv=0:  
    --> pv_star_2 >> pv_star_1 >> p0 --> qv_star_2 << qv_star_1


- s_1 negative for 





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

EddyField_output: 
    - compute Eddy Fields (fluctuations) from LES field output
    - save as new NC-files

Correlations: 
    - compute Correlations from Eddy Fields
    - output as NC-files 



# Visualization
Vis: plot from pickle

Vis_fields: 
    - plot from LES output fields and profiles
    - compute Correlations
    - compute 



# Thermodynamic
