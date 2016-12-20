???????????
- entropies.py: sv_c() -- correct to set it to zero if pv=0?

- pv_star_1 > p0 --> qv_star_1<0 (bcs. in calculation qv_star_c: pd<0) --> qt>qv_star_1 (saturated even for qt=0) but pv=0 
    --> changed formulation for (saturated && pv==0) to s_1 = s (initial entropy)
- if continuing saturated thermodynamics with qt=0, pv=0:  
    --> pv_star_2 >> pv_star_1 >> p0 --> qv_star_2 << qv_star_1


- s_1 negative for 





# LES Files

parameters: LES parameters



# Statistics
EM_PDF: fit a Gaussian PDF by using the 'expectation-maximation' (EM) algorithm
        --> output means, covariances and relative weights
EM_PDF_univariate: fit a uni-variate Gaussian PDF by using the 'expectation-maximation' (EM) algorithm
        --> output means, covariances and relative weights
EM_PDF_plot: read in output from EM_PDF.py / EM_PDF_univariate.py & plot PDFs from Gaussian
EM_PDF_stochastic: Fit an auto-regression (AR) model to the EM PDF parameters
    - arfit_py.py
    - arqr_py.py: QR decomposition


# IO
read_in_netcdf4: example for reading in NetCDF4 files




# Offline Calculations
budget_TKE: offline TKE budget computation from LES output (old LES; hdf5 output)
budget_Th: offline potential temperature budget (potential energy)

compatibility: check compatibility of scalar tracer convergence according to Chr. Schär





# Fluctuations and Correlations
EddyField_output: 
    - compute Eddy Fields (fluctuations) from LES field output
    - save as new NC-files

Correlations: 
    - compute Correlations from Eddy Fields
    - output as NC-files 



# Visualization
Vis: plot from Visualization outputs (pickle files)

Vis_fields: 
    - plot from LES output fields and profiles
    - compute Correlations
    - compute 



# Thermodynamics
thermodynamics_1: offline saturation adjustment with (Clausius Clapeyron by Magnus Formula; saturation adjustment) as in LES;
                    normal loop --> gives nan in T
thermodynamics_2: same as thermodynamics_1, but with correction on first iteration loop (qv_star_1 > 0)
thermodynamics_3: same as thermodynamics_1, but with additional corrections on all iteration loop (qv_star_i > 0, or equivalently T<373K or pd>0)
thermodynamics_4: same as thermodynamics_1, but with correction in T_2 --> if pv_star_2>p0: T_2 = 371K
    - T_2 <= 371K:
            if pv_star_2[i,j] >= p0:
                T_2[i,j] = 371.0 K
                ...
    - T_2 <= 350K:
            if T_2[i,j] >= 350:
                T_2[i,j] = 350.0 K
                ...

CC_Magnus: saturation adjustment 'step by step', Clausius Clapeyron with Magnus formula and from PyCLES
CC_Magnus_mod: like CC_Magnus, but with modified T_2 formulation (if T_2>373.3K: set T_2=373.3K)

thermo_plot: simple file for plotting thermodynamic relations (e.g. defined in thermo_aux)




