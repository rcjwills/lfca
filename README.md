# lfca
Low-frequency component analysis (LFCA) is a method that transforms the leading empirical orthogonal functions (EOFs) 
of a data set in order to identify a pattern with the maximum ratio of low-frequency to total variance (based on 
application of a lowpass filter). The resulting low-frequency patterns (LFPs) and low-frequency components (LFCs) 
isolate low-frequency climate variability and are useful in diagnosing the corresponding mechanisms. This method is 
presented in Wills et al. (2018, GRL). 

Here, we provide an example LFCA script in Matlab. The script run_lfca_example.m runs an example that creates Fig. 4 
of the associated paper. By changing the value of truncation in the script to 3, the script will create Fig. 1 of the 
paper instead. The method is contained within lfca.m. Included data is from the NOAA Extended Reconstructed Sea-Surface 
Temperature data set (Smith et al. 2008).

Reference for Method:
Wills, R.C., T. Schneider, J.M. Wallace, D.S. Battisti, and D.L. Hartmann, 2018: Disentangling global warming, multidecadal variability, and El Niño in Pacific temperatures. Geophysical Research Letters, 45, doi:10.1002/2017GL076327. [PDF] [SI] [Official version]

Reference for Data Used:
Smith, T.M., R.W. Reynolds, T.C. Peterson, and J. Lawrimore, 2008: Improvements to NOAA’s historical merged land–ocean surface temperature analysis (1880–2006). Journal of Climate, 21 (10), 2283–2296. [Official version]
