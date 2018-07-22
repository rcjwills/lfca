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
