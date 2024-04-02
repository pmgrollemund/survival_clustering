# Caenorhabditis elegans Survival Comparison Analysis (under construction)

## Description
This repository contains R source code and data for replicating the analysis described in the article "A clustering-based survival comparison procedure designed to study the Caenorhabditis elegans Model". The analysis focuses on comparing survival curves using a novel procedure specifically designed for the Caenorhabditis elegans model.

## Features
- Implementation of survival comparison procedures in R tailored for the Caenorhabditis elegans model.
- Data preprocessing scripts to clean and prepare the survival data.
- Statistical modeling scripts to perform survival analysis and comparison.
- Customization flexibility: Modify the data file paths in the script call to analyze different datasets, without changing the R scripts. 
- Multiple analysis scripts:
  1. **Survival Curve Clustering**: Classifies survival curves based on the work "The discriminative functional mixture model for a comparative analysis of bike sharing systems" of C. Bouveyron, E. Côme et J. Jacques.
  2. **Variability Analysis of Nematode Counting Errors**: Analyzes the variability of counting errors in nematodes.
  3. **Data Augmentation**: Augments the database using simulated data that considers counting error variability.


  
## Installation
1. Clone this repository to your local machine.
2. Ensure you have R (version 4.1.2 or higher) installed on your system.
3. The following R packages are required to run the scripts:
   - rjson
   - tictoc
   - cobs
   - fda, funFEM
   - ggplot2, ggpubr, ggthemes
   - mclust
   - MLmetrics
   - purrr
   - survival
   - truncnorm
   - utils

   You can install them by executing `install.packages("package_name")` in the R console.



## Usage
Before running any scripts, make sure to modify the `wd` option in the file `options.json`. The path specified should correspond to the root directory of this Git repository.
When using the scripts with a new dataset, ensure that you correctly modify the entries in the file `specific_options_test.json` to match the specifics of your dataset.

### Running the Scripts via Terminal/Command Prompt
- Open your terminal/command prompt.
- Navigate to the directory containing the R scripts.
- Run the appropriate R scripts for the analysis you wish to perform using the `Rscript` command, for example:
  - For Counting Errors Variability analysis: `Rscript Variability.R ../data/calibration_data.ods`
  - For Survival Clustering analysis: `Rscript pipeline_survival_clustering.R options/options.json options/specific_options_test.json`

### Using RStudio
- Open RStudio.
- Open the desired R script (`Variability.R` or `pipeline_survival_clustering.R`).
- Execute the script line by line or by selecting the entire script and running it.
- Make sure that the 'args' object, at the start of the script, contains the required values, i.e. the paths to the option files or database files.

### Understanding Options
#### For `options.json`:
- `"wd"`: Working directory for R (the root directory of this Git repository).<sup>1</sup>
- `"verbose"`: Indicates whether the scripts should display information.
- `"script_functions"`: Path to an R script containing necessary functions for the procedure.<sub>*Do not modify unless necessary*</sub>
- `"choice_functions"`: Path to an R script for assigning appropriate functions based on user choices in the options file.
- `"error_distribution_file"`: Contains the result of the `Variability.R` function that analyzes counting error variability.
- `"type_data"`: Specifies the file extension of the data file.
- `"type_transformation"` and `"type_pretreatment"`: Indicate data preprocessing methods.
- `"n_cut"`: Number of periods to segment the time domain of survival curves.
- `"n_functional_base"`: Number of basis functions to use for the B-spline basis.
- `"simulation"`: Indicates whether to use simulated data.
- `"simulated_data.time_number"`: Number of evaluation time points for simulated data.
- `"alpha_simulated_data"`: Transparency level of simulated data in graphical results.
- `"type_clust"`: Method of clustering for functional data.
- `"criteria_funFEM"`: Criterion for selecting the best modeling of clustered survival data.
- `"n_rep"`: Number of repetitions for clustering method.
- `"min_K"` and `"max_K"`: Minimum and maximum number of groups to consider for clustering.


#### For `specific_options_test.json`:
- `"data_file"`: Path to the file containing the dataset.<sup>1</sup>
- `"sheet"`: Name of the sheets to import.<sup>1</sup>
- `"control_file"`: Path to the file containing control group data. If the same as `"data_file"`, use "same".<sup>1</sup>
- `"control_sheets"`: Sheets corresponding to the control group.<sup>1</sup>
- `"treat_group"`: Sheets related to different experimental conditions.<sup>1</sup>
- `"treat_names"`: Names assigned to each experimental condition.<sup>1</sup>
- `"folder_path"`: Path where analysis results will be saved.

<sup>0</sup> *Do not modify unless necessary.* 

<sup>1</sup> *Make sure to set this option to the correct directory to avoid errors.* 

## Data
Please note that the dataset included in this repository is a reduced and anonymized version of the actual data used for the article. 

## Authors
- [Grollemund Paul-Marie](https://github.com/pmgrollemund/)
- Comte Elise

## License
This project is licensed under a [Creative Commons License](https://creativecommons.org/) allowing free reuse of the research work. 
