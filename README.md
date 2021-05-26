# lfca
Low-frequency component analysis (LFCA) is a method that transforms the leading empirical orthogonal functions (EOFs) 
of a data set in order to identify a pattern with the maximum ratio of low-frequency to total variance (based on 
application of a lowpass filter). The resulting low-frequency patterns (LFPs) and low-frequency components (LFCs) 
isolate low-frequency climate variability and are useful in diagnosing the corresponding mechanisms. This method is 
presented in Wills et al. (2018, GRL). 

Here, we provide an example LFCA script in Matlab and Python. The scripts run_lfca_example.m or run_lfca_example.py run an example that creates Fig. 4 of the associated paper (Wills et al. 2018). By changing the value of truncation in the script to 3, the script will create Fig. 1 of this paper instead. The method is contained within lfca.m or the function lfca within signal_processing.py. 

Reference for Method:
Wills, R.C., T. Schneider, J.M. Wallace, D.S. Battisti, and D.L. Hartmann, 2018: Disentangling global warming, multidecadal variability, and El Niño in Pacific temperatures. Geophysical Research Letters, 45, doi:10.1002/2017GL076327.

The data file for the example includes data from the NOAA Extended Reconstructed Sea-Surface 
Temperature data set (Smith et al. 2008) and can be downloaded here: https://atmos.uw.edu/~rcwills/data/ERSST_1900_2016.mat. All use of this data should reference the appropriate publication:

Smith, T.M., R.W. Reynolds, T.C. Peterson, and J. Lawrimore, 2008: Improvements to NOAA’s historical merged land–ocean surface temperature analysis (1880–2006). Journal of Climate, 21 (10), 2283–2296.

The data file for the bivariate example includes data from the AVISO monthly-mean sea-surface height anomaly data set (https://www.aviso.altimetry.fr/en/index.php?id=1526) and can be downloaded here: https://atmos.uw.edu/~rcwills/data/AVISO_1993_2016.mat. All use of this data should use the citation: "This product was processed by SSALTO/DUACS and distributed by AVISO+ (https://www.aviso.altimetry.fr) with support from CNES”

