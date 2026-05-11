# greta_val

Repository for analyses and figure production conducted by NAK and LH for Greta Aeby and Val Paul. Data were collected by Greta Aeby in Belize across three timepoints: the years 2019, 2022, and 2024. Data were cleaned by LH. Data themselves were unmanipulated; the only cleaning was to unify site naming conventions that changed over the course of the experiment. 

## Repository Structure

The repository itself is split into three main sections: data, output, and scripts. Data and outputs are organized in subfolders. Scripts are kept in the main branch to simplify pathing. 

### "data" folder

The data folder contains the two main datasets we were provided for analyses: 

1. benthic_cover.csv contains cover data for broad benthic categories (hard coral, soft coral, sponges, macroalgae, and other) that were collected in-situ using the point intercept method. Point intercept counts are converted to proportional cover (i.e. ranging from 0 to 1 instead of 0 to 100) in the 01_benthic cover script.

2. coral_demographics.csv contains coral density data from targeted coral surveys. It contains counts of coral colonies identified to the lowest taxonomic level.

However, there is a subfolder called "raw" which contains the raw coral species counts. This is only called in 00_clean_demo_data.R, which cleans the site names, dates, etc. and produces the coral_demographics.csv dataset that is used downstream.

### "output" folder

The output folder contains all figures and tables that were produced in the scripts. We use a fully reproducible workflow so that figure production should be exactly replicated across devices and operating systems.

### Scripts

Each script is clearly labeled and can be run independently. Scripts are meant to be run in an RProject from a cloned version of this repository. Pathing reflects this assumption and may need to be updated if the project is not cloned. 

00_clean_demo_data.R (author: LH) cleans the raw coral density data, producing "coral_demographics.csv"

01_benthic_cover.R (author: NAK) runs analyses and figure production for the broad benthic cover data ("benthic_cover.csv").

02_demographic_data.R (author: NAK) runs all univariate analyses of coral community composition based on "coral_demographics.csv". This script includes the GLMM of beta dispersion.

03_coral_multivariate.R (author: LH) runs the nMDS, PERMANOVA, and SIMPER analyses on the coral community, and produces the relevant figures and tables in the "output" folder.
