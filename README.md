# Segmentation-Lasso-Survival

This repository contains all code necessary to reproduce the simulation study presented in the manuscript "Unsupervized identification of prognostic copy-number
alterations using segmentation and lasso regularization"

# Data

Data generated for the simulation study is available under the Data directory. Scripts used to produce this data are available in the scripts/ directory, as 1_Data.R

# Processing

Different methods compared are implemented under the Methods.R and RunMethods.R scripts. They are applied on the simulated data in files 2a and 2b. Their results are then interpreted with the Criteria.R functions in the 3_synthesis scripts, and plotted in the 4_plots.R file.
Intermediary results are available in the Data directory.
The user should be warned that running the methods on the "long" simulation scenarios is computationnaly intensive, and requires large memory space. We do not recommend using those scripts on a personal computer.
