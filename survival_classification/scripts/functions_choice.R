################################################################################
################ Determine which functions to use for 
################   "Survival curve classification"
######## Elise Comte and Paul-Marie Grollemund
####### 2021-05-06 : first implementation
####### 2024-04-02 : make the script reusable
################################################################################

# Some objects
if(control_file == "same"){
  control_file <- data_file
}

# Some pathes 
Clustering_path <- paste0(results_path,"/Clustering.RData") 
Proba_path <- paste0(results_path,"/Proba.RData")

# import_data 
if (type_data == "xlsx") import_data <- import_xlsx 

# pretreatments ----
if (type_pretreatment == "spline") {
  if (simulation) {
    pretreatment <- pretreat_splines_simuldata
  } else {
    pretreatment <- pretreat_splines # do not use
  }
}

# transformation ----
if (type_transformation == "deviance") transformation <- deviance_transformation 

# clustering ----
if (type_clust == "funFEM") {
  if (criteria_funFEM == "ICL") do_clustering <- funFEM_ICL
}

# proba_cluster ----
if (type_clust == "funFEM") Proba_cluster <- Pb_funFEM
