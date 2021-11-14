# FF-roost-structure
Data and code repository for Lunn et al (2021) Counterintuitive scaling between population abundance and local density: implications for modelling transmission of infectious diseases in bat populations

To empirically evaluate the relationship between total roost abundance and tree-level measures of abundance and density – the scale most likely to be relevant for virus transmission – we collected a 13-month dataset of tree-roosting Pteropus spp. from 2,522 spatially referenced trees across eight roosts, and evaluated whether roost features at different scales (roost-level, subplot-level, tree-level) were predictive of these local density dynamics. 

The code is organised into an RProject containing four files:
- 00_functions_VSubmission.R - behind-the-scenes custom functions to help with data processing and visualisation
- 01_create-kernel-density-and-Voroni-polygons_VSubmission.R - spatial code to 1) create kernel density estimates from tree spatial point patterns weighted by bat index values in trees, and 2) calculate and visualise crown area of trees in subplots
- 02_format-raw-data_VSubmission.R - code to format the raw data into structures suitible for analyes, and to save processed data for later reference
- 03_fit-models_VSubmission.R - code to run general addative mixed models (GAMMs) addressing the main aim
- 04_generate-figures_VSubmission.R - code to generate figures presented in the published manuscript
** It is important to load the code as the RProject so that root folders in setwd are correct **

Data is organised into an RProject containing two folders
- /Raw : raw data input for code '01_create-kernel-density-and-Voroni-polygons_VSubmission.R' and '02_format-raw-data_VSubmission.R'
- /Processed : processed data input for code '03_fit-models_VSubmission.R' and '04_generate-figures_VSubmission.R'
More information on each file and data within files is given in the data_README file

Output is organised into an RProject containing:
- /Figures : figures presented in the main text or supplementary information in the published paper
- /Model output : GAM model output for both response variables (Tree-level 3-D density and Tree-level abundance) saved as:
     - Raw GAM output (files beginning with MODEL_)
     - Comparison tables of candidate models, ranked by AIC (files ending in _comparison)
     - Coefficient values from best ranked models (containing model_output)
Model output is organised by response variable, then bat species (main directory is species combined, sub-directories are BFF=black flying-fox, GHFF=grey-headed flying-fox, LRFF=little red flying-fox.
Models are fitted generalized additive mixed models (GAMMs) with restricted maximum likelihood (REML) estimation, Poisson (tree-level abundance) or Gamma (tree-level 3-D density) error distributions (both with a log link) and random effects of roost site, subplot and survey session with the mgcv package in R 
More details on model structure are given in the published manuscript
