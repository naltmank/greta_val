# greta_val

Repository for analyses and figure production conducted by NAK and LH for Greta Aeby and Val Paul. Data were collected by Greta Aeby in Belize across three timepoints: the years 2019, 2022, and 2024. Data were cleaned by LH in a separate script not replicated here. Data themselves were unmanipulated; the only cleaning was to unify site naming conventions that changed over the course of the experiment. 

## Repository Structure

The repository itself is split into three main sections: data, output, and scripts. Data and outputs are organized in subfolders. Scripts are kept in the main branch to simplify pathing. 

### "data" folder

The data folder contains the two main datasets we were provided for analyses: 

1. benthic_cover.csv contains cover data for broad benthic categories (hard coral, soft coral, sponges, macroalgae, and other) that were collected in-situ using the point intercept method. Point intercept counts are converted to proportional cover (i.e. ranging from 0 to 1 instead of 0 to 100) in the 01_benthic cover script.

2. coral_demographics.csv contains coral density data from targeted coral surveys. It contains counts of coral colonies identified to the lowest taxonomic level.

### "output" folder

The output folder contains all figures that were produced in the scripts. We use a fully reproducible workflow so that figure production should be exactly replicated across devices and operating systems.

### Scripts

Each script is clearly labeled and can be run independently. Scripts are meant to be run in an RProject from a cloned version of this repository. Pathing reflects this assumption and may need to be updated if the project is not cloned. 

01_benthic_cover.R runs analyses and figure production for the broad benthic cover data ("benthic_cover.csv").

02_demographic_data.R runs all univariate analyses of coral community composition based on "coral_demographics.csv". This script includes the GLMM of beta dispersion.
