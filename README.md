# forced-patterns

Forced-pattern filtering is a method to identify spatial patterns (linear combinations of empirical orthogonal functions (EOFs)) with the maximum ratio of signal to noise (with signal defined as variance that is agreed upon across an ensemble). The resulting forced patterns (FPs) isolate the patterns of forced change within climate model ensembles (where each simulation is subject to the same external forcing). This method is presented in Wills et al. (2020, J. Climate, in review).

Here, we provide an example FP filtering script in Matlab and Python (COMING SOON), with application to the CESM Large Ensemble (Kay et al. 2015). The script run_FP_filtering_example.m runs an example that creates Figs. 1 and 5 of the associated paper (Wills et al. 2020). The method is contained within forced_pattern_analysis.m.

The data file for the example includes seasonal-mean surface temperature data from the CESM Large Ensemble (Kay et al. 2015) and can be downloaded here: atmos.uw.edu/~rcwills/data/CESM-LE_hist-rcp85_1920-2019_3monthly_ts.mat. Note that this data has been reduced to a resolution of ~5°, and there are small deviations from the published figures for this reason. All use of this data should reference the appropriate publication:

Kay, J., and Coauthors, 2015: The Community Earth System Model (CESM) large ensemble project: A community resource for studying climate change in the presence of internal climate variability. Bull. Am. Meteorol. Soc., 96 (8), 1333–1349.

Reference for Method: Wills, R.C.J., D.S. Battisti, K.C. Armour, T. Schneider, and C. Deser, 2020: Pattern recognition methods to separate forced responses from internal variability in climate model ensembles and observations, in review at Journal of Climate.
