# OASIS3_LMEanalyses
Code used to analyze longitudinal amyloid deposition and brain network function using OASIS-3 dataset.


**amyloid_analyses.m**

A script that imports data, subsets it based on availability of longitudinal data and presence of current neurological diagnoses other than probably AD/very mild cognitive impairment, calculates progressor vs non progressor labels, calculates mixed effect models (random intercepts and slopes for each subject) and uses FDR correction for the main effect of group and group x time since baseline interaction

**FC_analyses.m**

imports data, subsets it based on availability of longitudinal data, presence of current neurological diagnoses other than probably AD/very mild cognitive impairment, registration QC and mean FD calculates mixed effect models (random intercepts and slopes for each subject) for the main effect of group and group x time since baseline interaction

**permutation_tests.R.zip**

PLEASE UNZIP AND UPLOAD ACTUAL SCRIPTS
- permutation testing for functional connectivity mixed effect models

**amyloid_fittedmeanbrainplots.m**

scripts for calculating colormaps for mapping out amyloid levels in the Desikan-Killiany Freesurfer parcellation

**correlation_analyses.m**

scripts calculating the first principal component for amyloid deposition slopes and correlating these slopes with functional connectivity (FC) slopes
plotting scatterplots of the amyloid slopes against baseline age and against FC slopes


code from https://github.com/peterzhukovsky/fMRI_GLM_validation was modified for some of the above scripts used in these analyses. 
