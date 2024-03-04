# Caenorhabditis elegans Survival Comparison Analysis

## Description
This repository contains R source code and data for replicating the analysis described in the article "A New Survival Comparison Procedure Designed to Study the Caenorhabditis elegans Model." The analysis focuses on comparing survival rates using a novel procedure specifically designed for the Caenorhabditis elegans model.

## Features
- Implementation of survival comparison procedures in R tailored for the Caenorhabditis elegans model.
- Data preprocessing scripts to clean and prepare the survival data.
- Statistical modeling scripts to perform survival analysis and comparison.
- Customization flexibility: Modify the data file paths in the scripts to analyze different datasets.
- Multiple analysis scripts:
  1. **Survival Curve Classification**: Classifies survival curves based on a novel procedure.
  2. **Variability Analysis of Nematode Counting Errors**: Analyzes the variability of counting errors in nematodes.
  3. **Data Augmentation**: Augments the database using simulated data that considers counting error variability.


  
## Installation
1. Clone this repository to your local machine.
2. Ensure you have R (version 4.1.2 or higher) installed on your system.
3. The following R packages are required to run the scripts:
   - FactoMineR
   - factoextra
   - nFactors
   - readxl
   - tidyverse
   - missMDA
   - ggforce
   - rstatix
   - ggpubr
   - GGally

   You can install them by executing `install.packages("package_name")` in the R console.



## Usage
### Running the Scripts via Terminal/Command Prompt
- Open your terminal/command prompt.
- Navigate to the directory containing the R scripts.
- Run the appropriate R scripts for the analysis you wish to perform using the `Rscript` command, for example:
  - For PCA analysis: `Rscript PCA.R`
  - For ANOVA analysis: `Rscript ANOVA.R var_name` (e.g., `Rscript ANOVA.R GEN_AB`)
  - For correlation analysis: `Rscript correlation_analysis.R`

### Using RStudio
- Open RStudio.
- Open the desired R script (`PCA.R`, `ANOVA.R`, or `correlation_analysis.R`).
- Execute the script line by line or by selecting the entire script and running it.

## Data
Please note that the dataset included in this repository is a reduced and anonymized version of the actual data used for the article. 

## Authors
- [Grollemund Paul-Marie](https://github.com/pmgrollemund/)
- Comte Elise

## License
This project is licensed under a [Creative Commons License](https://creativecommons.org/) allowing free reuse of the research work. 
