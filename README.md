# CrossSectionFitting
Fitting codes for double Jpsi cross section measurement
Applied in CMSSW_12_4_0, cmsenv is required

## Order of the fitting
1. 1D fit to acquire the distributions on the ctau dimensions
2. Modify the ctau shape in the 4D fit code according to the last line
3. 4D fit to acquire the total number of different components, and the shape of the Jpsi mass peak
4. Modify the mass peak shape in the binned 4D fit code according to the last line
5. Binned 4D fit to acquire the number of different components in each bin (automatic scripts are provided)
6. Modify the signal numbers array in the template fit code according to the last line
7. Template fit to acquire the SPS/DPS fraction
8. Iteration of weight mixing shall be applied
