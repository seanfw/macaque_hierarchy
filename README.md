# Calculation of Macaque Cortical Hierarchy

This repository contains the R and MATLAB code and associated data used specifically for the calculation of the macaque cortical hierarchy as part of the analysis for the published paper:

Froudist-Walsh, Sean, Daniel P. Bliss, Xingyu Ding, Lucija Rapan, Meiqi Niu, Kenneth Knoblauch, Karl Zilles, Henry Kennedy, Nicola Palomero-Gallagher, and Xiao-Jing Wang. "A dopamine gradient controls access to distributed working memory in the large-scale monkey cortex." Neuron 109, no. 21 (2021): 3500-3520.

## Contributors
The code and data were prepared by Seán Froudist-Walsh and Kenneth Knoblauch.

## Original Connectivity Data
The original connectivity data used in this analysis are publicly available and can be downloaded from [corenets.org](http://www.corenets.org).

## Repository Overview
This repository consists of code to calculate hierarchical positions of brain regions based on their connectivity patterns, primarily using the beta-binomial model.

There are two key scripts in the repository: `betabinomial_hierarchy.m` and `hierarchy_regressions.R`. The former is a MATLAB script for pre-processing the connectivity data, creating the design matrix, and saving the data in an appropriate format for the hierarchical modeling. The latter is an R script which performs beta-binomial regressions to calculate the hierarchical positions of each brain region.

## betabinomial_hierarchy.m

This MATLAB script is used to pre-process the connectivity data, create the design matrix, and save the data for further hierarchical modeling.


## Files
The repository is structured as follows:

**`/raw_data`**:
- `infragranular.txt`: raw data file - number of cell bodies projecting from infragranular layers for each connection.
- `source.txt`: raw data file - ordered list of source areas for studied connections
- `inj_site.txt`: raw data file - ordered list of injected areas for studied connections
- `supragranular.txt`: raw data file - number of cell bodies projecting from infragranular layers for each connection.

**`/src`**:
- `betabinomial_hierarchy.m`: Main MATLAB script for processing the connectivity data.
- `hierarchy_regressions.R`: This R script reads in data fits a beta-binomial model to the data.

**`/processed_data`**:
- `FLN_SLN_89_areas.mat`: Processed connectivity data (FLN, SLN) for the 89-area cortical graph, with missing values for areas without injections.
- `beta_binomial_data.mat`: MATLAB data file used as input for the R script.
- `betaBinHierVals.mat`: MATLAB data file is the output of the R script, containing the coefficients of the beta-binomial model.
- `hierarchy_89_areas.mat`: estimated hierarchy values for the full 89-area cortical graph
- `beta_bin_hierarchy_subgraph.mat`: hierarchy and connectivity data (FLN/SLN) for the full 40-area subgraph


### How to use:

1. Make sure you have the necessary data files in the correct directories. The script expects text files for the target and source data, as well as for the supragranular and infragranular data, to be located in a directory named `raw_data`.

2. Run the script in MATLAB up to the comment "%% Calculate beta-binomial fit in R (hierarchy_regressions.R)". This will generate several figures and also save a `.mat` file containing the relevant data in the directory named `processed_data`.

## hierarchy_regressions.R

This R script performs a beta-binomial regression to estimate the hierarchical position of each brain region, based on the connectivity data pre-processed by `betabinomial_hierarchy.m`.

### How to use:

1. First, run `betabinomial_hierarchy.m` up to the comment "%% Calculate beta-binomial fit in R (hierarchy_regressions.R)" to generate the design matrix and save it in a `.mat` file.

2. Run `hierarchy_regressions.R` in R. Make sure you have the correct path to the `.mat` file generated in step 1. This script will output a set of hierarchical values for each region, which can then be used to order the regions in a hierarchical manner. The output of the script is saved in a `.mat` file in a directory named `processed_data`.

3. Run the remaining portion of `betabinomial_hierarchy.m`, using the output from `hierarchy_regressions.R`.

Note: Before running both scripts, please ensure that your MATLAB and R environments have necessary packages installed and paths correctly set. Also, the scripts assume a certain directory structure, so please ensure your data files are organized accordingly.

If you encounter any issues or need further clarification, please raise an issue on the repository.

## License
This code and data are available under a CC-BY-NC-4.0 license. If you use the code or data in your research, please cite our paper:

Froudist-Walsh, Sean, et al. "A dopamine gradient controls access to distributed working memory in the large-scale monkey cortex." Neuron 109, no. 21 (2021): 3500-3520.

## Contact
For any issues or concerns, please contact the authors.
