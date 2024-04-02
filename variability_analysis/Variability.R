################################################################################
################### Analyse the error count variability
######## Elise Comte - Paul-Marie Grollemund
####### 2021-07-19 : first implementation
####### 2023-05-23 : major update
####### 2024-04-02 : make the script reusable
################################################################################
### Usage from terminal :
# move to the folder of the script and run : 
#     Rscript Variability.R data_path
# where 'data_path' is the path of your calibration data file. 
# Let try : 
#     Rscript Variability.R ../data/calibration_data.ods
################################################################################
#### Clean up ----
rm(list=ls())

## Get shell options (manual) ----
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  # data file path
  data_file_path <- args[1]
}else{
  stop("You need to specify the path to the data file.")
}

#### Required packages ----
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(readxl))

#### Required functions ----
completion_NA <- function(df,df_weights){
  number_NA <- sum(is.na(df$freq))
  while(number_NA>0){
    number_NA_debut <- sum(is.na(df$freq))
    
    index_na <- which(is.na(df$freq))
    for(i in sample(index_na)){
      current <- df[i,]
      
      tmp <- df[
        df$mean %in% (df$mean[i]+c(-1,0,1)) &
          df$error %in% (df$error[i]+c(-1,0,1)),
      ]
      tmp <- tmp[!is.na(tmp$freq),]
      
      if(nrow(tmp) >0){
        distances <- abs(current$mean - tmp$mean) + abs(current$error - tmp$error) 
        weights <- rep(0, length(distances))
        for(j in 1:length(weights)){
          weights[j] <- df_weights$weight[df_weights$distance == distances[j]]  
        }
        
        df$freq[i] <- crossprod(tmp$freq,weights)
      }
      
    }
    number_NA <- sum(is.na(df$freq))
    if(number_NA_debut == number_NA) break
  }
  return(df)
}
smooth <- function(df,df_weights){
  df0 <- df
  
  for(i in sample(1:nrow(df))){
    current <- df[i,]
    
    tmp <- df0[
      df$mean %in% (df$mean[i]+c(-1,0,1)) &
        df$error %in% (df$error[i]+c(-1,0,1)),
    ]
    
    distances <- abs(current$mean - tmp$mean) + abs(current$error - tmp$error) 
    weights <- rep(0, length(distances))
    for(j in 1:length(weights)){
      weights[j] <- df_weights$weight[df_weights$distance == distances[j]]  
    }
    
    df$freq[i] <- crossprod(tmp$freq,weights)
    
  }
  
  means <- unique(df$mean)
  for(j in 1:length(means)){
    df$freq[df$mean == means[j]] <- 
      df$freq[df$mean == means[j]] / 
      sum(df$freq[df$mean == means[j]])
  }
  
  return(df)
}

#### Required objects ----
df_weights <- data.frame(
  distance = 0:5
)
df_weights$weight <- exp(-df_weights$distance)


#### Create the output dir ----
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path("..","results",paste("variability_analysis_run",timestamp,sep="_"))

if(!dir.exists(file.path("..","results"))) 
  dir.create(file.path("..","results"))

dir.create(output_dir)



#### Data importation ----
# Import ----
# Looking for sheets in calibration data file
data_file_sheets <- readODS::ods_sheets(data_file_path)
n_experimenters <- length(data_file_sheets)

# Initialize object
calibration_raw_data <- list(); length(calibration_raw_data) <- n_experimenters
names(calibration_raw_data) <- data_file_sheets

# Import raw data from each sheet 
for(i in 1:n_experimenters){
  calibration_raw_data[[i]] <- readODS::read_ods(data_file_path,sheet = data_file_sheets[i])
}


# Preteat ----
# Design a suitable data.frame 
calibration_data <- data.frame(
  count = NULL,
  mean = NULL,
  day = NULL,
  rep = NULL
)

for(i in 1:n_experimenters){
  # Extract count matrix
  counts_matrix <- calibration_raw_data[[i]][,-(1:2)]
  n_counts <- ncol(counts_matrix)
  
  # Compute mean count to estimate the true number of nematodes into the well
  mean_tmp <- round(rowMeans(counts_matrix))
  
  tmp <- data.frame(
    count = as.vector(t(as.matrix(counts_matrix))),
    mean = rep(mean_tmp,each=n_counts),
    day = rep(calibration_raw_data[[i]]$Day,each=n_counts),
    rep = rep(calibration_raw_data[[i]]$Rep,each=n_counts)
  )  
  
  calibration_data <- rbind(calibration_data,tmp)
}

# Estimate error counts ----
calibration_data$count_error <- calibration_data$count - calibration_data$mean


#### Compute error distribution ---- 
# Determine the range for error and mean ----
errors <- min(calibration_data$count_error):max(calibration_data$count_error)
means <- min(calibration_data$mean):max(calibration_data$mean)


# Initialize the distribution object ----
error_distribution <- expand.grid(
  mean=means,error=errors
)
error_distribution$freq <- rep(0,nrow(error_distribution))

# Increment the error distribution ----
for(i in 1:nrow(error_distribution)){
  index <- which(calibration_data$mean == error_distribution$mean[i] & 
                   calibration_data$count_error == error_distribution$error[i]
  )
  error_distribution$freq[i] <- length(index)
}

# Normalize each marginal distribution ----
for(mean in means){
  error_distribution$freq[error_distribution$mean == mean] <- 
    error_distribution$freq[error_distribution$mean == mean] / 
    sum(error_distribution$freq[error_distribution$mean == mean])
}

# Treat with missing values ----
error_distribution$freq[error_distribution$freq %in% c(0,NaN)] <- NA


df <- error_distribution
complete_error_distribution <- completion_NA(error_distribution,df_weights)
smooth_error_distribution <- smooth(complete_error_distribution,df_weights)
smooth_error_distribution <- smooth(smooth_error_distribution,df_weights)


#### Graphical outputs ---- 
# Compute the plot
p <- ggplot(smooth_error_distribution, aes(mean, error, fill= freq)) + 
  geom_tile() + scale_fill_viridis()+ 
  theme_pubclean(base_size = 15)

# Save the plot
ggsave(p,filename = file.path(output_dir,"error_distribution.pdf"),device = "pdf",
       width = 9,height = 6)


# Save the distribution 
save(smooth_error_distribution,file = file.path(output_dir,"error_distribution.RData"))








