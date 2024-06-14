################################################################################
################ Replica validation
####### Paul-Marie Grollemund
####### Muriel Bonnet, Guy Febvre
####### Franck Peyre, Vincent Romera 
####### 20??-??-?? : first implementation 
####### 2022-06-22 : v1.4, replica validation
#######                shows sample sizes, and provides "*"-code according to the p-value
####### 2022-06-30 : v1.5, works with an additional sheet "initialisation"
#######                Control group: "Statut1_Cond" and other groups: "Statut2_Cond"
#######                Also provides color (Couleur_Condition) and a numerical code (Num_Condition)
####### 2024-05-28 : v2.0, program redesign
#######                remove only problematic pairs "Puits"x"Replica"
####### 2024-06-05 : v2.0.1, change legends for some plots
####### 2024-06-10 : v2.0.2, Add function for removing special characters in file 
#######                name and add of Axis X name and step in "TO BE CHANGED"
####### 2024-06-14 : v2.1.0, include option in the data file
################################################################################
### Usage :
# Be sure the folder of the script contains the file "Fiche_Comptage.ods"
#    This file must consists of the sheets: "Resultats", "Initialisation" and 
#    "Options".
# Can also be run throught terminal with the Rscript command 
################################################################################

################################################################################
#### Clean up ----
################################################################################
rm(list=ls())

################################################################################
#### Initialization ----
################################################################################

#/-----------------------------------------------------------------------------\
#------------ START : TO BE CHANGED
#\-----------------------------------------------------------------------------/
# Working directory ----
# write here the working directory where you wish results would be saved
wd <- NULL 
if(!is.null(wd)) setwd(wd)

# Data file name ----
data_file_name <- "../data/random_survival.ods"
data_sheet <- "Resultats"
initialization_sheet <- "Initialisation"
option_sheet <- "Options"
#/-----------------------------------------------------------------------------\
#------------ END : TO BE CHANGED
#\-----------------------------------------------------------------------------/

# Verbose options 
version <- "2.1.0"
part <- paste(rep("-",50),collapse = "")
section <- paste(rep("-",10),collapse = "")
task <- paste(rep("-",2),collapse = "") 
subtask <-  paste(rep("-",1),collapse = "")

# Major options ----
default_option <- data.frame(
  line_size = 0.5,
  threshold_pvalue = 0.05,
  verbose = TRUE,
  height = 10,
  width = 20,
  breakTime = 2,
  max_x_ticks = 60,
  graphics_title_comparison = "Comparison of all conditions",
  xlab = "Time (Day)",
  ylab = "Survival probability"
)

# Import data
options <- readODS::read_ods(path = data_file_name, sheet = option_sheet,skip = 1)

# Assign values
for(i in 1:nrow(options)){
  assign(options$Nom[i],options$Valeur[i])
}

# Assign default values if values are not provided 
for(i in 1:length(default_option)){
  if( !exists(names(default_option)[i]) ){
    assign(names(default_option)[i],default_option[[i]])
  }else{
    if( is.na(get(names(default_option)[i])) ){
      assign(names(default_option)[i],default_option[[i]])
    }
  }
}
breakTime_save <- breakTime

# Launch ----
if(verbose)
  cat(part,"\n",
      "Replica validation \n",
      " version ",version," \n",
      'Run on "',data_file_name,'"\n',
      part,"\n\n",sep="")

# Initialization ----
if(verbose) cat(section,"Initialization \n")

# Required initial packages ----
if(verbose) cat(task,"Import packages \n")
suppressPackageStartupMessages(library("readODS"))
suppressPackageStartupMessages(library("survival"))
suppressPackageStartupMessages(library("survminer"))
suppressPackageStartupMessages(library("dplyr"))

# Create output dir ----
if(verbose) cat(task,"Create an output directory \n")
current_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path(getwd(),
                        paste("Result_Replica_Validation_",current_time,sep="")
)
if(!dir.exists(output_dir))
  dir.create(path = output_dir)

validation_output_dir <- file.path(output_dir,"Validation")
if(!dir.exists(validation_output_dir))
  dir.create(path = validation_output_dir)
survival_curves_output_dir <- file.path(output_dir,"Survival_curves")
if(!dir.exists(survival_curves_output_dir))
  dir.create(path = survival_curves_output_dir)
comparison_output_dir <- file.path(output_dir,"Survival_comparison")
if(!dir.exists(comparison_output_dir))
  dir.create(path = comparison_output_dir)

# Required functions ----
if(verbose) cat(task,"Import functions \n")

# Removing special characters
clean_filename <- function(filename) {
  cleaned <- gsub("[^[:alnum:]_]", "_", filename)
  return(cleaned)
}

# "replica_structure" is an object informing about the strata "Replica x Puits"
# the following function enables to change this object according to which 
# replicas are removed 
update_replica_structure <- function(replica_structure,puits_to_remove,replica_to_remove){
  index_puit <- which(names(replica_structure) == puits_to_remove)
  
  replica_structure[[index_puit]][replica_to_remove] <- NA
  replica_structure[[index_puit]][-(1:replica_to_remove)] <- 
    replica_structure[[index_puit]][-(1:replica_to_remove)] - 1
  # replica_structure[[index_puit]] <- 
  #   replica_structure[[index_puit]][!is.na(replica_structure[[index_puit]])]
  
  if(length(replica_structure)>index_puit){
    for(k in (index_puit+1):length(replica_structure))
      replica_structure[[k]] <- replica_structure[[k]] - 1
  }
  
  return(replica_structure)
}

# Compute a "*"-code for a given pvalue
get_pvalue_code <- function(pvalue){
  if(pvalue <= 0.001) return("***")
  if(pvalue <= 0.01)  return(" **")
  if(pvalue <= 0.05)  return("  *")
  if(pvalue <= 0.1)   return("  .")
  if(pvalue >  0.1)   return("   ")
}

# Define a transparent palette from a color
make_transparent_palette <- function(color,n){
  palette <- rep(color,n)
  palette_rgb <- col2rgb(palette)
  for(k in 1:length(palette)){
    palette[k] <- 
      rgb(palette_rgb[1,1],palette_rgb[2,1],palette_rgb[3,1],
          max=255,alpha=((k+1)*100/(length(palette)+1))*255/100  
      )
  }
  return(palette)
}

# Plot survival curves
SurvPlot<-function(pv, fit, title, legendsTitle, legendsLabs, palette, xlim, 
                   breakTime,line_size,xlab,ylab)
{
  p<-ggsurvplot(fit,
                pval = pv,
                pval.size = 4,
                pval.coord = c(xlim/2, 0.9),
                conf.int = FALSE,
                title = title,
                xlab = xlab,
                ylab = ylab,
                surv.median.line = 'hv',
                ggtheme = theme_bw(),
                legend.title = legendsTitle,
                legend.labs = legendsLabs,
                palette = palette,
                legend=c(0.85,0.65),
                xlim=c(0,xlim),
                break.time.by= breakTime,
                font.x=c("bold"),
                font.y=c("bold"),
                font.tickslab=c("bold"),
                size=line_size
  )
  return(p)
}

if(verbose) cat('\n')
################################################################################
#### Data and sheets ----
################################################################################
if(verbose) cat(section,"Data and sheets \n")

# Import data ----
if(verbose) cat(task,"Import data \n")
if(verbose) cat(subtask,"Sheet: Resultats \n")
data <- read_ods(path = data_file_name, sheet = data_sheet)
data <- as.data.frame(data)
nrow_data <- nrow(data) 

# About initialisation sheet ----
if(verbose) cat(subtask,"Sheet: Initialisation \n")

# Import  ----
init <- read_ods(path = data_file_name, sheet = initialization_sheet)
init <- as.data.frame(init)
nrow_init <- nrow(init)

# Determine role of the conditions ----
if(verbose) cat(subtask,"Determine roles of the conditions specified in the sheet 'Initialisation.'\n")
message_role_condition <- NULL

# Determine the condition to remove ----
conditions_to_remove <- init$Nom_Condition[which(init$Statut1_Condition == 0)]

if(length(conditions_to_remove)>0){ # message if there is at least a condition to remove
  data <- data[!(data$Souche %in% conditions_to_remove),]  
  nrow_data <- nrow(data) 
  message_tmp <- paste("\t The following conditions are removed: ",
                       paste0("\n\t\t",conditions_to_remove,collapse=""),"\n",
                       sep="")
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}else{ # message if there is no condition to remove
  message_tmp <- "\t Condition to remove: none.\n "
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}

# Get the control group ----
control_condition <- init$Nom_Condition[init$Statut1_Condition == "t"]

if(length(control_condition) == 0){ # message if there is no control group found
  message_tmp <- "\t None control group found.\n"
  message_role_condition <- c(message_role_condition,message_tmp)
  write(message_role_condition,file.path(output_dir,"Role_condition.txt"))
  stop(message_tmp)
}else{ # message indicating the control group found
  message_tmp <- paste("\t Control condition: \n\t\t",control_condition,"\n",sep="")
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}

# Determine the condition to work with ----
conditions_to_treat <- init$Nom_Condition[which(init$Statut1_Condition != 0)]

if(length(conditions_to_treat) == 0){ # message if all condition are removed
  message_tmp <- "All available conditions are removed.\n\n"
  message_role_condition <- c(message_role_condition,message_tmp)
  write(message_role_condition,file.path(output_dir,"Role_condition.txt"))
  stop(message_tmp)
}else{ # message to indicate the conditions to compare to the control group
  condition <- init[init$Nom_Condition %in% conditions_to_treat,1:2]
  nrow_init <- nrow(condition)
  message_tmp <- paste("\t Conditions to compare to the control group: \n\t\t",
                       paste(condition[condition[,2] != control_condition,2],"\n\t\t",collapse=""),
                       "\n",
                       sep="")
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}

# Determine the condition to validate ----
if("to_validate" %in% colnames(init)){
  condition_to_not_validate <- init$Nom_Condition[init$to_validate == 0]
  condition_to_not_validate <- 
    condition_to_not_validate[condition_to_not_validate %in% conditions_to_treat]
  
  if(length(condition_to_not_validate)>0){
    message_tmp <- "\t The following condition will not be validated:"
    message_tmp <- paste(message_tmp,
                         paste("\n\t\t",condition_to_not_validate,collapse=""),
                         "\n",
                         sep="")
    
    if(verbose) cat(message_tmp)
  }
  condition_to_validate <- init$Nom_Condition[init$to_validate == 1 & init$Statut1_Condition != 0]
}else{
  condition_to_validate <- conditions_to_treat
  condition_to_not_validate <- NULL
}


# Determine the condition to work with ----
conditions_to_pairwise_comparison <- init$Nom_Condition[which(init$Statut2_Condition == 1)]
conditions_to_pairwise_comparison <- 
  conditions_to_pairwise_comparison[conditions_to_pairwise_comparison %in% conditions_to_treat]

if(length(conditions_to_pairwise_comparison) <2){ # message if there is no pairwise comparison to do
  message_tmp <- paste("\t No pairwise comparison, since less than 2 condition have the value '1' in ",
                       "the column 'Statut2_Condition' of the sheet 'Initialisation'.\n\n",
                       sep="")
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}else{ # message to indicate the condition to compair in one-to-one comparison
  condition <- init[init$Nom_Condition %in% conditions_to_treat,1:2]
  nrow_init <- nrow(condition)
  
  message_tmp <- paste("\t Conditions for pairwise comparison: \n\t\t",
                       paste(conditions_to_pairwise_comparison,"\n\t\t",collapse=""),
                       sep="")
  message_role_condition <- c(message_role_condition,message_tmp)
  if(verbose) cat(message_tmp)
}

# Save messages about condition role ----
write(message_role_condition,file.path(output_dir,"Role_condition.txt"))
if(verbose) cat('\n')

################################################################################
#### Replica Validation ----
################################################################################
if(verbose) cat(section,"Replica validation \n")
message_replica_validation <- NULL

# Determine replicas to remove ----
# For each experimental conditions
data$to_remove <- rep(NA,nrow_data)
if(!is.null(condition_to_not_validate))
  data$to_remove[data$Souche %in% condition_to_not_validate] <- FALSE

for(condition in condition_to_validate){
  # Get data about this condition 
  data_condition <- data[data$Souche == condition,]
  n_puits <- length(unique(data_condition$Puits))
  n_replica_puits <- sapply(unique(data_condition$Puits),
                            function(puit) 
                              length(unique(data_condition[data_condition$Puits == puit,"Replica"])
                              )
  )
  replica_structure <- lapply(n_replica_puits,
                              function(n) 1:n)
  for(i in 2:length(replica_structure)) 
    replica_structure[[i]] <- replica_structure[[i]] + max(replica_structure[[i-1]])
  
  # Determine if there is an overall difference between replicas (Puits x Replica)
  pval <- survdiff(formula=Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition)
  replica_removed <- NULL
  pvalues <- pval$pvalue
  contributions_replica <- NULL
  
  # If there is a difference, it is required to remove some replica that 
  # are very different from the majority. To do that, we consider remove step 
  # by step the replica related to the high contribution to the pvalue.
  # Recall that the contribution is related to (O-E)^2/V
  while(pval$pvalue < threshold_pvalue){
    # Compute the contributions the statistic value
    contributions <- (pval$obs-pval$exp)^2 / diag(pval$var)
    
    # Determine the one which mostly contributes
    index_replica <- which.max(contributions)
    
    # Determine the contribution of this replica to the p-value
    contribution_replica <- 
      round(contributions[index_replica] / sum(contributions) * 100,2)
    
    # Determine Puit and Replica for the aforementioned replica
    index_replica <- lapply(replica_structure,
                            function(vec) which(vec == index_replica)
    )
    puits_to_remove <- names(unlist(index_replica))
    replica_to_remove <- unlist(index_replica)
    
    replica_structure <- 
      update_replica_structure(replica_structure,puits_to_remove,replica_to_remove)
    
    
    # Remove the raw data related to the replica to remove
    data_condition <- data_condition[ 
      !(data_condition$Replica == replica_to_remove & data_condition$Puits == puits_to_remove)
      ,]
    
    # Determine if there is an overall difference between replicas (Puits x Replica)
    pval <- survdiff(formula=Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition)
    pvalues <- c(pvalues,pval$pvalue)
    contributions_replica <- c(contributions_replica,contribution_replica)
    
    replica_removed <- rbind(replica_removed,
                             c(replica_to_remove,puits_to_remove))
  }
  
  # Mark replicas to remove
  data_condition <- data[data$Souche == condition,]
  
  replica_puits <- paste( data_condition$Replica, data_condition$Puits)
  replica_puits_removed <- paste(replica_removed[,1],replica_removed[,2])
  data_condition$removed <- as.factor(replica_puits %in% replica_puits_removed)
  
  data$to_remove[data$Souche == condition] <- as.character(data_condition$removed)
  
  # If some replicas were removed 
  if(!is.null(replica_removed)){
    # Display the number of the replica for the user
    message_tmp <- paste(subtask," '",condition,"': The following replicas are removed.\n",
                         sep="")
    for(k in 1:nrow(replica_removed))
      message_tmp <- paste(message_tmp,
                           '\t Replica ',replica_removed[k,1],' - ',replica_removed[k,2],
                           "\n\t    - contribution to anomalous heterogeneity: ",contributions_replica[k],"%",
                           "\n\t    - p-value before removing: ",get_pvalue_code(pvalues[k]),
                           " ",pvalues[k],"\n",
                           sep="")
    message_tmp <- paste(message_tmp,
                         "\t-> Final p-value: ",get_pvalue_code(pvalues[k+1]),
                         " ",pvalues[k+1],"\n\n",
                         sep="")
    message_replica_validation <- c(message_replica_validation,message_tmp)
    if(verbose) cat(message_tmp)
    
    # Plot the survival curves to check 
    n_replica_structure <- length(unlist(replica_structure))
    tmp <- which( !is.na(unlist(replica_structure)))
    palette <- rep("red",n_replica_structure)
    palette[tmp] <- "gray"
    labs <- rep("removed",n_replica_structure)
    labs[tmp] <- "validated"
    
    fit <- survfit(Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition)
    time_max <- max(data_condition$Duree_de_Vie,na.rm=T)
    
    if(time_max / breakTime > max_x_ticks){
      breakTime <- ceiling(time_max/max_x_ticks)
    }else{
      breakTime <- breakTime_save
    }
    
    plot_removed_replicas <- 
      ggsurvplot(fit,
                 legend="right",
                 palette = palette,
                 title = paste(condition,": Replicas validation"),
                 conf.int = FALSE,
                 xlab = xlab,
                 legend.title = "Removed : red",
                 # legend.labs = c("Removed","Keep"),
                 ylab = ylab,
                 ggtheme = theme_bw(),
                 xlim=c(0,time_max),
                 break.time.by= breakTime,
                 font.x=c("bold"),
                 font.y=c("bold"),
                 font.tickslab=c("bold"),
                 size=line_size
      )
    
    fit <- survfit(Surv(Duree_de_Vie, Statut)~removed, data=data_condition)
    time_max <- max(data_condition$Duree_de_Vie,na.rm=T)
    
    if(time_max / breakTime > max_x_ticks){
      breakTime <- ceiling(time_max/max_x_ticks)
    }else{
      breakTime <- breakTime_save
    }
    
    title <- paste(condition,": Replication validation (average)")
    plot_removed_replicas_average <- 
      SurvPlot(TRUE, fit, title,"", c("validated","removed"), 
               c("gray","red"), time_max, breakTime,line_size = line_size,
               xlab,ylab
      )
    
    # Save the plots
    file_name <- paste("Validation_replicas_",clean_filename(condition),".pdf",sep="") #Add of use of removing special characters function
    suppressWarnings(suppressMessages(
      ggsave(plot=plot_removed_replicas$plot,device = "pdf",
             filename = file.path(validation_output_dir,file_name),
             height = height,width=width,units = "cm"
      )
    ))
    
    file_name <- paste("Validation_replicas_average_",clean_filename(condition),".pdf",sep="") #Add of use of removing special characters function
    suppressWarnings(suppressMessages(
      ggsave(plot=plot_removed_replicas_average$plot,device = "pdf",
             filename = file.path(validation_output_dir,file_name),
             height = height,width=width,units = "cm")
    ))
  }else{
    # Indicate that none replicas were removed
    message_tmp <- paste(subtask," '",condition,"': All Replicas are validated.",
                         "\n\t-> p-value: ",get_pvalue_code(pval$pvalue),
                         " ",round(pval$pvalue,8),"\n\n",
                         sep="")
    message_replica_validation <- c(message_replica_validation,message_tmp)
    if(verbose) cat(message_tmp)
  }
}

# Definitely remove replicas ----
data <- data[data$to_remove != TRUE,]
data$to_remove <- NULL

# Save messages about replica validation ----
write(message_replica_validation,file.path(validation_output_dir,"Replica_validation.txt"))

################################################################################
#### Survival curve for each condition ----
################################################################################
if(verbose) cat(section,"Survival curve for each condition \n")
message_survival_curve <- NULL

# Get survival curves ----
# For each experimental conditions
time_max <- max(data$Duree_de_Vie,na.rm = T)+1

if(time_max / breakTime > max_x_ticks){
  breakTime <- ceiling(time_max/max_x_ticks)
}else{
  breakTime <- breakTime_save
}

for(condition in conditions_to_treat){
  
  # Get data about this condition 
  data_condition <- data[data$Souche == condition,]
  fit <- survfit(Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition)
  
  color <- init$Couleur_Condition[init$Nom_Condition == condition]
  palette <- make_transparent_palette(color,length(fit$strata))
  
  # Plot the survival curves
  survival_curves_condition <- 
    ggsurvplot(fit,legend="none",
               title = condition,
               palette = palette,
               pval.coord = c(time_max/2, 0.9),
               conf.int = FALSE,
               xlab = xlab,
               ylab = ylab,
               ggtheme = theme_bw(),
               xlim=c(0,time_max),
               break.time.by= breakTime,
               font.x=c("bold"),
               font.y=c("bold"),
               font.tickslab=c("bold"),
               size=line_size
    )
  
  # Save the plot
  file_name <- paste("Survival_curves_",clean_filename(condition),".pdf",sep="") #Add of use of removing special characters function
  suppressWarnings(suppressMessages(
    ggsave(plot=survival_curves_condition$plot,device = "pdf",
           filename = file.path(survival_curves_output_dir,file_name),
           height = height,width=width,units = "cm")
  ))
  
  # Summary of the survival curve : to save
  n_caractere <- 
    length(unlist(strsplit(survival_curves_condition[["plot"]][["labels"]][["title"]],split="")))
  
  message_survival_curve <- 
    c(message_survival_curve,
      capture.output({
        cat("\n/",rep("-",35-n_caractere)," ", survival_curves_condition[["plot"]][["labels"]][["title"]],
            " ", rep("-",35),"\\\n",sep="")
        print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition))
        cat("\n","------- Summary -------\n")
        print(summary(fit)$table)
        cat("\\",rep("-",70),"/\n\n",sep="")
      })
    )
  
  # Summary of the survival curve : to print
  if(verbose){
    cat("\n/",rep("-",35-n_caractere)," ", survival_curves_condition[["plot"]][["labels"]][["title"]],
        " ", rep("-",35),"\\\n",sep="")
    print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Puits+Replica, data=data_condition))
    cat("\n","------- Summary -------\n")
    print(summary(fit)$table)
    cat("\\",rep("-",70),"/\n\n",sep="")
  }
}

# Save messages about survival curves ----
write(message_survival_curve,file.path(survival_curves_output_dir,"Summary.txt"))

################################################################################
#### Condition comparison (control) ----
################################################################################
if(verbose) cat(section,"Condition comparison to the control group \n")
message_control_comparison <- NULL

# Determine the other groups ----
other_conditions <- init$Nom_Condition[init$Statut1_Condition == 1]

# Determine length of condition names
n_caractere_control <- 
  length(unlist(strsplit(control_condition,split="")))

# Define an object to save the pvale codes
pval_codes <- data.frame(condition = init$Nom_Condition[init$Nom_Condition %in% other_conditions])
pval_codes$pval_code <- rep(NA,nrow(pval_codes))
pval_codes$sample_size_condition <- rep(NA,nrow(pval_codes))

# For each comparison "group x control-group" ----
for(condition in other_conditions){
  data_comparison <- data[data$Souche %in% c(control_condition,condition),]
  
  # Define the message
  n_caractere <- 
    length(unlist(strsplit(condition,split="")))
  
  message_control_comparison <- 
    c(message_control_comparison,
      capture.output({
        cat("\n/",rep("-",35-n_caractere-n_caractere_control)," ", 
            condition," <<--->> ",control_condition,
            " ", rep("-",35-11),"\\\n",sep="")
        print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison))
        cat("\\",rep("-",70),"/\n\n",sep="")
      })
    )
  if(verbose){
    cat("\n/",rep("-",35-n_caractere-n_caractere_control)," ", 
        condition," <<--->> ",control_condition,
        " ", rep("-",35-11),"\\\n",sep="")
    print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison))
    cat("\\",rep("-",70),"/\n\n",sep="")
  }
  
  
  # Paste sample size and pvalue code to labels
  sample_size_control_condition <- 
    nrow(data_comparison[data_comparison$Souche == control_condition,])
  sample_size_condition <- nrow(data_comparison[data_comparison$Souche == condition,])
  pval <- survdiff(formula=Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison)
  pval_code <- get_pvalue_code(pval$pvalue)
  
  pval_codes$pval_code[pval_codes$condition == condition] <- pval_code
  pval_codes$sample_size_condition[pval_codes$condition == condition] <- sample_size_condition
  
  # Get label and color with the same order as in init sheet
  color_labels_order <- 
    init[init$Nom_Condition %in% c(control_condition,condition),c("Nom_Condition","Couleur_Condition")]
  labels_order <- color_labels_order$Nom_Condition
  palette <- color_labels_order$Couleur_Condition
  
  # Paste sample size and pvalue code to labels
  labels_order[labels_order != control_condition] <- 
    paste(labels_order[labels_order != control_condition], 
          " (n=",sample_size_condition,") ",pval_code,sep="")
  
  new_control_condition <- paste(labels_order[labels_order == control_condition], 
                                 " (n=",sample_size_control_condition,")",sep="")
  labels_order[labels_order == control_condition] <- new_control_condition  
  
  
  # Fit the condition comparison
  data_comparison$Souche[data_comparison$Souche %in% condition] <- labels_order[labels_order != new_control_condition]
  data_comparison$Souche[data_comparison$Souche %in% control_condition] <- labels_order[labels_order == new_control_condition]
  data_comparison$Souche <- factor(data_comparison$Souche,levels=labels_order)
  fit <- survfit(Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison)
  time_max <- max(data_comparison$Duree_de_Vie,na.rm=T)
  title <- paste(condition, " x ",control_condition,sep="")
  
  if(time_max / breakTime > max_x_ticks){
    breakTime <- ceiling(time_max/max_x_ticks)
  }else{
    breakTime <- breakTime_save
  }
  
  # Plot the survival curves
  survival_comparison <- suppressWarnings(
    SurvPlot(TRUE, fit, title,"Condition", labels_order, 
             palette, time_max, breakTime,line_size,
             xlab,ylab)
  )
  
  # Save the plot
  file_name <- paste("Control_comparison_",clean_filename(condition),".pdf",sep="") #Add of use of removing special characters function
  suppressWarnings(suppressMessages(ggsave(plot=survival_comparison$plot,device = "pdf",
                                           height = height,width=width,units = "cm",
                                           filename = file.path(comparison_output_dir,
                                                                file_name))))
}

# Save messages about control comparison ----
write(message_control_comparison,file.path(comparison_output_dir,"Control_comparison.txt"))


# Get label and color with the same order as in init sheet
color_labels_order <- init[init$Nom_Condition %in% data$Souche,c("Nom_Condition","Couleur_Condition")]
labels_order <- color_labels_order$Nom_Condition
palette <- color_labels_order$Couleur_Condition

# Paste sample size and pvalue code to labels
labels_order[labels_order != control_condition] <- 
  paste(labels_order[labels_order != control_condition], 
        " (n=",pval_codes$sample_size_condition,") ",pval_codes$pval_code,sep="")
color_labels_order$new_label <- labels_order

# Compare all conditions to control group ----
data_tmp <- data
for(condition in other_conditions)
  data_tmp$Souche[data$Souche %in% condition] <- color_labels_order$new_label[color_labels_order$Nom_Condition==condition]
data_tmp$Souche[ data_tmp$Souche %in% control_condition ] <- 
  new_control_condition
labels_order[labels_order == control_condition] <- new_control_condition  

data_tmp$Souche <- factor(data_tmp$Souche,levels=labels_order)
fit <- survfit(Surv(Duree_de_Vie, Statut)~Souche, data=data_tmp)
time_max <- max(data_tmp$Duree_de_Vie,na.rm=T)

if(time_max / breakTime > max_x_ticks){
  breakTime <- ceiling(time_max/max_x_ticks)
}else{
  breakTime <- breakTime_save
}

# Plot the survival curves ----
survival_comparison <- 
  ggsurvplot(fit,
             legend=c(0.85,0.65),
             title = graphics_title_comparison,
             palette=palette,
             pval.coord = c(time_max/2, 0.9),
             conf.int = FALSE,
             xlab = xlab,
             ylab = ylab,
             legend.title = "Condition",
             legend.labs = labels_order,
             ggtheme = theme_bw(),
             xlim=c(0,time_max),
             surv.median.line = 'hv',
             break.time.by= breakTime,
             font.x=c("bold"),
             font.y=c("bold"),
             font.tickslab=c("bold"),
             data=data,
             size=line_size
  )

# Save the plot ----
file_name <- "All_condition_comparison.pdf"
suppressWarnings(suppressMessages(ggsave(plot=survival_comparison$plot,device = "pdf",
                                         height = height,width=width,units = "cm",
                                         filename = file.path(comparison_output_dir,
                                                              file_name))))

################################################################################
#### Pairwise condition comparison (no control) ----
################################################################################
if(verbose) cat(section,"Pairwise Condition comparison \n")
message_control_comparison <- NULL

other_conditions <- init$Nom_Condition[
  init$Statut2_Condition == 1 & init$Statut1_Condition == 1
]
new_other_conditions <- color_labels_order$new_label[color_labels_order$Nom_Condition %in% other_conditions]

if(length(other_conditions) > 1){
  # For each comparison "group x other group"
  for(i in 1:(length(other_conditions)-1)){
    condition1 <- other_conditions[i]
    new_condition1 <- new_other_conditions[i]
    for(j in (i+1):length(other_conditions)){
      condition2 <- other_conditions[j]
      new_condition2 <- new_other_conditions[j]
      
      data_comparison <- data[data$Souche %in% c(condition1,condition2),]
      
      n_caractere1 <- 
        length(unlist(strsplit(condition1,split="")))
      n_caractere2 <- 
        length(unlist(strsplit(condition2,split="")))
      
      message_control_comparison <- 
        c(message_control_comparison,
          capture.output({
            cat("\n/",rep("-",35-n_caractere1-n_caractere2)," ", 
                condition1," <<--->> ",condition2,
                " ", rep("-",35-11),"\\\n",sep="")
            print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison))
            cat("\\",rep("-",70),"/\n\n",sep="")
          })
        )
      if(verbose){
        cat("\n/",rep("-",35-n_caractere1-n_caractere2)," ", 
            condition1," <<--->> ",condition2,
            " ", rep("-",35-11),"\\\n",sep="")
        print(survdiff(formula=Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison))
        cat("\\",rep("-",70),"/\n\n",sep="")
      }
      
      fit <- survfit(Surv(Duree_de_Vie, Statut)~Souche, data=data_comparison)
      time_max <- max(data_comparison$Duree_de_Vie,na.rm=T)
      title <- paste(condition1, " x ",condition2,sep="")
      
      color1 <- init$Couleur_Condition[init$Nom_Condition == condition1]
      color2 <- init$Couleur_Condition[init$Nom_Condition == condition2]
      palette <- c(color1,color2)
      labels <- c(condition1,condition2)
      palette <- palette[order(labels)]
      labels <- sort(labels)
      
      if(time_max / breakTime > max_x_ticks){
        breakTime <- ceiling(time_max/max_x_ticks)
      }else{
        breakTime <- breakTime_save
      }
      
      # Plot the survival curves
      survival_comparison <- 
        SurvPlot(TRUE, fit, title,"Condition", labels, 
                 palette, time_max, breakTime,line_size,
                 xlab,ylab)
      
      # Save the plot
      file_name <- paste("Pairwise_comparison_",clean_filename(condition1),"_",condition2,".pdf",sep="")
      suppressWarnings(suppressMessages(ggsave(plot=survival_comparison$plot,device = "pdf",
                                               height = height,width=width,units = "cm",
                                               filename = file.path(comparison_output_dir,
                                                                    file_name))))
    }
  }
  
  # Save messages about control comparison ----
  write(message_control_comparison,file.path(comparison_output_dir,"Pairwise_comparison.txt"))
}else{
  if(verbose) cat("None pairwise comparison required.\n")
}
if(verbose) cat("\n")

################################################################################
#### End messages ----
################################################################################
if(verbose) cat(section,"End messages \n")

# Save validated ----
if(verbose) cat(task,"Save Validated data \n")

write_ods(data,sheet = data_sheet,
          path = file.path(output_dir,"Validated_data.ods"))
write_ods(init,sheet = initialization_sheet,append = TRUE,
          path = file.path(output_dir,"Validated_data.ods"))

# Save directory----
if(verbose) cat(task,"Save directory \n")
cat("Results and messages are saved in the directory :\n",output_dir,"\n",
    part,"\n"
)

