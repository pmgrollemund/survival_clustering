################################################################################
################ Survival curve classification 
####### Elise Comte and Paul-Marie Grollemund
####### 2021-04-29 : first implementation
####### 2023-05-10 : major update
####### 2024-04-02 : make the script reusable
################################################################################
### Usage from terminal :
# move to the folder of the script and run : 
#     Rscript pipeline_survival_clustering.R options specific_options
# where 'options' is the path of your general options for the run of the script
# and 'specific_options' is for the specific options related to your dataset.
# Let try : 
#     Rscript pipeline_survival_clustering.R options/options.json options/specific_options_test.json
################################################################################

#### Clean up ----
rm(list=ls())

# Required initial packages ----
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(tictoc))

## Get shell options (manual) ----
args <- commandArgs(trailingOnly=TRUE)
# setwd("~/Documents/MCF/Recherche/UMRF/Elegans/Survival_classification/git_survival_clustering/survival_classification/")
# args <- c("./options/options.json","./options/specific_options_test.json")
# setwd("~/Documents/MCF/Recherche/UMRF/Elegans/Survival_classification/git_survival_clustering/survival_classification/")
# args <- c("./options/options_article.json","./options/article/specific_options_Cardin_OP_VS_OPAF_3.json")
if(length(args)>1){
  # general options
  json_file <- args[1]
  json_data <- fromJSON(file=json_file)
  for(i in 1:length(json_data)){
    assign(names(json_data)[i], json_data[[i]])
  }
  # specific options
  json_file_bis <- args[2]
  json_data_bis <- fromJSON(file=json_file_bis)
  for(i in 1:length(json_data_bis)){
    assign(names(json_data_bis)[i], json_data_bis[[i]])
  }
}else{
  stop("You need to specify the path to the options files. \nCheck the folder 'options'.")
}
setwd(wd)

#### Initialization ----
source("scripts/initialization.R")

# Launch ----
if(verbose){
  cat(part,"\n",
      "Survival curve clustering \n",
      " version ",version," \n",
      part,"\n",sep="")
}
tic("Running time")


#### Initialization ----
if(verbose) cat(section," Initialization \n",sep="")

# Required packages ----
if(verbose) cat(task," Import packages \n",sep="")
suppressPackageStartupMessages(library(cobs))
suppressPackageStartupMessages(library(fda,quietly = TRUE))
suppressPackageStartupMessages(library(funFEM))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(mclust,quietly=TRUE))
suppressPackageStartupMessages(library(MLmetrics))
suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(truncnorm))
suppressPackageStartupMessages(library(utils))

# Create output dir ----
if(verbose) cat(task," Create output directory \n",sep="")
if(!dir.exists(folder_path))
  dir.create(path = folder_path)

results_path <- paste0(folder_path,"/res")
if(!dir.exists(results_path))
  dir.create(path = results_path)

img_path <- paste0(folder_path,"/img")
if(!dir.exists(img_path))
  dir.create(path = img_path)

# Required functions ----
if(verbose) cat(task," Import functions \n",sep="")
source(script_functions)
source(choice_functions)

#### Data ----
if(verbose) cat(section," Database \n",sep="")

# Import data ----
if(verbose) cat(task," Import data \n",sep="")
if(verbose) cat(subtask,"Read from database \n")
data <- import_data(data_file, sheet)
control <- import_data(control_file, control_sheets)

# Pretreat ----
if(verbose) cat(task,"Pretreat \n")
if(verbose) cat(subtask,"Time \n")
time_max <- compute_time_max(data)
Time <- seq(0,time_max,length.out = simulated_data.time_number)

if(verbose) cat(subtask,"Survival curves \n")
data_pretreated <- pretreatment(data)
if(verbose) cat(subtask,"Survival curves (control group) \n")
control_pretreated <- pretreatment(control,data_pretreated,control_sheets) 

# Transformation ----
if(verbose) cat(subtask,"Data transforming : survival differences \n")
data_transformed <- transformation(data_pretreated,control_pretreated)

#### Clustering ----
if(verbose) cat(section,"Clustering \n")
if(verbose) cat(task,"Compute clusters \n")
data_cluster <- do_clustering(data_transformed)

if(verbose) cat(task,"Compute clusters probabilities \n")
final_data <- Proba_cluster(data_cluster,threshold_deviance)

#### Post-process ----
if(verbose) cat(section,"Post-process \n")

# Graphical results ----
if(verbose) cat(task,"Export graphical results \n")
compute_graphical_output()

#### Clean up ----
if(verbose) cat(section,"Ending \n")
if(verbose) cat(task,"Save \n")
Path <- paste(results_path,"/Param_curves.RData",sep="") 
Path_image <- paste(results_path,"/Renv.RData",sep="") 
save(data_pretreated, final_data, file = Path)
save.image(file = Path_image)

if(verbose) cat(task,"Clean up \n")
rm(list=ls())
