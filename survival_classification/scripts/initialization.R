################################################################################
################ Initialization for "Survival curve classification"
####### Elise Comte and Paul-Marie Grollemund
####### 2023-05-04 : first implementation
################################################################################
# Version -----
version <- "1.2.1"

#### Verbose function ----
replaceMessage <- function(x, width = 80) # 80
{
 message("\r", rep(" ", times = width - length(x)), "\r", appendLF = F)
 message(x, appendLF = F)
}

# Verbose -----
part <- "######################################"
section <- "####"
subsection <- "##"
task <- "-" 
subtask <- "\t"

# Get current timestamp ----
current_timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
# current_timestamp <- "test"

folder_path <- paste(folder_path,"_",current_timestamp,sep="")

# Working directory ----
# ---- Should only be used for testing/debugging 
# setwd("/home/pmgrolle/Documents/MCF/Recherche/UMRF/Elegans/Survival_classification/code/classification/")
