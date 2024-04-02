################################################################################
################ Functions for "Survival curve classification"
######## Elise Comte and Paul-Marie Grollemund
####### 2021-05-06 : first implementation
####### 2024-04-02 : make the script reusable
################################################################################

#### Data importation ----
import_xlsx <- function(path,sheet) {
  data <- list() ; length(data) <- length(sheet)
  if (length(sheet)>0){
    for (i in 1:length(sheet)) {
      data[[i]] <- readxl::read_excel(path, sheet = sheet[i])
    }
  }
  names(data) <- sheet
return(data)
}

#### Useful functions ----
# Return the maximal time accros data
compute_time_max <- function(data) {
  res_max <- 0
  for (i in 1:length(data)) {
    res_max <- max(res_max,max(data[[i]]$time))
  }
  return(res_max)
}

# Computhe Kaplan Meier estimate
KM <- function (data){
  return(survminer::surv_fit(Surv(time, status)~rep.intra, data=data))
}

# Return a ggsurvplot from KM model
fit_plot <- function(data){
  p <- 
    survminer::ggsurvplot(KM(data),
                          pval = FALSE,
                          pval.size = 2,
                          conf.int = FALSE,
                          title = "",
                          xlab = "time (d)",
                          ylab = "survival probability",
                          surv.median.line = 'hv',
                          legend.title = "Replicate",
                          break.time.by=2,
                          font.x=c("bold"),
                          font.y=c("bold"),
                          font.tickslab=c("bold"))
  return(p)
}

#### Test simulation ----
# Change format for the observed  data
real_dataframe <- function(data) {
  # Wells number 
  nWells <- 0
  for (i in unique(data$rep.inter)) {
    nWells <- nWells + length(unique(data$rep.intra[which(data$rep.inter == i)]))
  }
  
  # Days number
  time_max <- max(data$time)
  
  # rep.intra and rep.inter
  rep.intra <- list() ; length(rep.intra) <- (time_max + 2) * nWells
  rep.inter <- list() ; length(rep.inter) <- (time_max + 2) * nWells
  ind <- 1
  for (i in unique(data$rep.inter)) {
    dat_tmp <- data[which(data$rep.inter == i),]
    ind_positive <- ind+(time_max+2)*length(unique(dat_tmp$rep.intra))-1
    rep.intra[ind:ind_positive] <- 
      rep(unique(dat_tmp$rep.intra), each = time_max+2)
    rep.inter[ind:ind_positive] <- i
    ind <- ind_positive + 1
  }
  
  # Define curves
  curves <-
    data.frame(
      time = rep(0:(time_max + 1), nWells),
      n_nematode = rep(NA, (time_max + 2) * nWells),
      status = rep(data$status[1], (time_max + 2) * nWells),
      strain = factor(rep(data$strain[1], (time_max + 2) * nWells)),
      Wells = factor(rep(1:nWells, each = (time_max + 2))),
      rep.intra = factor(unlist(rep.intra)),
      rep.inter = factor(unlist(rep.inter)),
      num_duplic = "observed_data"
    )
  
  # Fill the object
  rep_intra <- unique(data$rep.intra)
  rep_inter <- unique(data$rep.inter)
  Wells <- 1
  for (i in 1:length(unique(data$rep.inter))) {
    data_rep.inter <- data[which(data$rep.inter == unique(data$rep.inter)[i]),]
    for (j in 1:length(unique(data_rep.inter$rep.intra))) {
      data_tmp <-
        data_rep.inter[data_rep.inter$rep.intra == unique(data_rep.inter$rep.intra)[j], ]
      time <- unique(data_tmp$time)
      
      # First days
      index_negative <- which(curves$time == 0 & curves$Wells == Wells)
      index_positive <- which(curves$time == data_tmp$time[1] & curves$Wells == Wells)
      curves$n_nematode[index_negative:index_positive] <- max(data_tmp$nematode_nb)
      n_deaths <- length(which(data_tmp$time == data_tmp$time[1])) # PMG : 2023-05-22 me semble etre une erreur et mal compter le nombre de morts au 1er instant
      
      # Other days
      for (k in 1:(max(1,length(time) - 1))) {
        index_negative <- index_positive + 1
        index_positive <-
          which(curves$time == time[k + 1] & curves$Wells == Wells)
        if(length(index_positive)>0){
          curves$n_nematode[index_negative:index_positive] <-
            curves$n_nematode[index_negative - 1] - n_deaths
        }
        n_deaths <- length(which(data_tmp$time == time[k + 1]))
      }
      
      # Last day
      index_negative <- index_positive + 1
      index_positive <- which(curves$time == (time_max + 1) & curves$Wells == Wells)
      if(length(index_negative)>0){
        curves$n_nematode[index_negative:index_positive] <- 0
      }
      Wells <- Wells + 1
    }
  }
  
  # Deal with NA values
  index_na <- which(is.na(curves$n_nematode))
  if(length(index_na)>0){
    for(i in index_na) curves$n_nematode[i] <- curves$n_nematode[i-1]
  }
  
  # Return result
  curves$Wells <- factor(curves$Wells)
  return(curves)
}

# Nematode count simulation
sim_n_nematode <- function(n_nematode,smooth_error_distribution){
  error <- 0
  if(n_nematode >= 35){
    error <- 
      sample(smooth_error_distribution$error[smooth_error_distribution$mean == 35],
             1,
             prob = smooth_error_distribution$freq[smooth_error_distribution$mean == 35])
  }
  if(n_nematode < 35 &  n_nematode >= 10){
    error <- 
      sample(smooth_error_distribution$error[smooth_error_distribution$mean == n_nematode],
             1,
             prob = smooth_error_distribution$freq[smooth_error_distribution$mean == n_nematode])
  }
  if(n_nematode < 10 & n_nematode > 0){
    error <- sample(c(-1,0,1),1,prob=c(0.2,0.6,0.2))
  }
  if(n_nematode == 0){
    error <- sample(c(0,1),1,prob=c(0.95,0.5))
  }
  
  new_n_nematode <- n_nematode + error
  if(new_n_nematode<0) new_n_nematode <- 0
  
  return(new_n_nematode)
}

# Generate a dataset
simul_data <- function(data, n) {
  if (n == 0) {
    return(data)
  }
  # load(error_sd_file)
  load(error_distribution_file)
  eps <- 0.2
  time_size <- nrow(data)
  Simul <-
    data.frame(
      time = rep(data$time, n),
      n_nematode = rep(NA, n*time_size),
      status = rep(data$status[1], n*time_size),
      strain = rep(data$strain[1], n*time_size),
      Wells = rep("duplic", each = n*time_size),
      rep.intra = rep(data$rep.intra[1], n*time_size),
      rep.inter = rep(data$rep.inter[1], n*time_size),
      num_duplic = rep(1:n, each = time_size)
    )
  
  Simul <- rbind(Simul, data)
  Simul$num_duplic <- factor(Simul$num_duplic)
  for (i in 1:n) {
    
    # First day
    n_nematode <- data$n_nematode[1]
  
    new_n_nematode <- sim_n_nematode(n_nematode,smooth_error_distribution)
    Simul$n_nematode[(i - 1)*time_size + 1] <- new_n_nematode
    
    # Other days
    for (j in 2:time_size) {
      n_nematode_before <- Simul$n_nematode[(i - 1)*time_size + j - 1]
      n_nematode <- data$n_nematode[j]
      
      new_n_nematode <- sim_n_nematode(n_nematode,smooth_error_distribution)
      
      if(n_nematode_before < new_n_nematode){
        Simul$n_nematode[(i - 1)*time_size + 1:(j-1)][
          Simul$n_nematode[(i - 1)*time_size + 1:(j-1)] < new_n_nematode
        ] <- new_n_nematode
      }
      Simul$n_nematode[(i - 1)*time_size + j] <- new_n_nematode
    }
  }
  return(Simul)
}

# Change the format of the dataframe
reformat_dataframe <- function(Data) {
  
  # Initialize 
  newdata <- data.frame(
    nematode_nb = 1:Data$n_nematode[1],
    time = rep(NA, Data$n_nematode[1]),
    status = rep(Data$status[1], Data$n_nematode[1]),
    strain = rep(Data$strain[1], Data$n_nematode[1]),
    rep.intra = rep(Data$rep.intra[1], Data$n_nematode[1]),
    rep.inter = rep(Data$rep.inter[1], Data$n_nematode[1]),
    num_duplic = rep(Data$num_duplic[1], Data$n_nematode[1])
  )
  nb_nematode <- unique(Data$n_nematode)
  
  # For each nematode
  ind_time_negative <- 1
  for (i in 2:length(nb_nematode)) {
    time <- Data$time[which(Data$n_nematode == nb_nematode[i])][1] - 1
    n_deaths <- nb_nematode[i-1] - nb_nematode[i]
    ind_time_positive <- ind_time_negative + n_deaths - 1
    newdata$time[ind_time_negative:ind_time_positive] <- time
    ind_time_negative <- ind_time_positive + 1
  }
  
  # Result
  return(newdata)
}

# Change the format of the dataframe
change_data <- function(Data) {
  
  # Initialize
  n <- n_simulated_data
  simulated_data <- list() ; length(simulated_data) <- length(Data)
  
  # For each experimental condition
  for (i in 1:length(Data)) {
    real_data <- real_dataframe(Data[[i]])
    for (j in unique(real_data$Wells)) {
      data_by_Wells <- real_data[which(real_data$Wells == j),]
      data_by_Wells$n_nematode[data_by_Wells$n_nematode<0] <- 0
      if (names(Data)[[i]] %in% control_sheets) { # PMG 2024-04-02 
        simulation <- simul_data(data_by_Wells,n) # PMG 2023-05-19 : j'ai mis n à la place de 0
      } else {
        simulation <- simul_data(data_by_Wells,n)
      }
      
      for (k in unique(simulation$num_duplic)) {
        data_tmp <- simulation[which(simulation$num_duplic == k),]
        
        # PMG (2023-05-19) : rajoute quelque chose au cas ou le nbre de nematode soit 
        #        constant ce qui pose un pb dans le calcul
        if( all(data_tmp$n_nematode == data_tmp$n_nematode[1]) ){
          data_tmp$n_nematode[length(data_tmp$n_nematode)] <- 
            max(data_tmp$n_nematode[length(data_tmp$n_nematode)] - 2,0)
          data_tmp$n_nematode[1] <- data_tmp$n_nematode[1]+1
        }
        
        reformat <- reformat_dataframe(data_tmp)
        if (k == unique(simulation$num_duplic)[1]) {
          reformat_simulation <- reformat
        } else {
          reformat_simulation <- rbind(reformat_simulation, reformat)
        }
      }
      if (j == unique(real_data$Wells)[1]) {
        newdata <- reformat_simulation
      } else {
        newdata <- rbind(newdata, reformat_simulation)
      }
      simulated_data[[i]] <- newdata
    }
  }
  
  # Result
  names(simulated_data) <- names(data)
  return(simulated_data)
}

#### Pretreatment ----
# Estimate a non-increasing spline version of the survival curve 
MSpline <- function(Data) {
  len <- length(unique(Data$rep.intra))
  pred <- matrix(NA, nrow = (simulated_data.time_number+1), ncol = len)
  MSE_curves <- list() ; length(MSE_curves) <- len
  for (i in 1:len){
    d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$time,
                    y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$surv)
    
    if(nrow(d) == 1){ # PMG 2023-05-23 : je rajoute en cas de pb 
      j <- sample((1:len)[-i],1) 
      d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[j]),]$time,
                      y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[j]),]$surv)
    }
    if(nrow(d) == 2){ # PMG 2023-05-23 : je rajoute en cas de pb
      d <- rbind(d[1,],
                 c(0.3*d$x[1] + 0.7 *d$x[2],0.3*d$y[1] + 0.7 *d$y[2]),
                 c(0.7*d$x[1] + 0.3 *d$x[2],0.7*d$y[1] + 0.3 *d$y[2]),
                 d[2,])
    }
    if(nrow(d) == 3){ # PMG 2023-05-23 : je rajoute en cas de pb
      d <- rbind(d[1,],
                 c(0.5*d$x[1] + 0.5 *d$x[2],0.5*d$y[1] + 0.5 *d$y[2]),
                 d[2,],d[3,])
    }
    
    # Constraints matrix 
    Constraints <- matrix(NA, nrow = 2*simulated_data.time_number+2, ncol = 3)
    for (j in 1:simulated_data.time_number) {
      # <= 1 
      Constraints[j,] <- c(-1,Time[j],1)
      # >= 0 
      Constraints[(simulated_data.time_number + j),] <- c(1,Time[j],0)
    }
    # = 1 at time 0
    Constraints[2*simulated_data.time_number+1,] <- c(0,0,1)
    # = 0 at last time
    Constraints[2*simulated_data.time_number+2,] <- c(0,time_max,0)
    ## Interpolation
    suppressWarnings(Cspline <-
                       cobs(d$x, d$y, constraint = "decrease", lambda = 0, 
                            degree = 2, pointwise = Constraints,
                            print.warn = FALSE, print.mesg = FALSE))
    preds <- predict(Cspline, interval="none", z = Time)
    pred[-(simulated_data.time_number+1),i] <- preds[,2]
    pred[simulated_data.time_number+1,i] <- MSE(Cspline$fitted,d$y)
  }
  return(pred)
}




#### Pretreat Splines avec simul_data
pretreat_splines_simuldata <- function(Data,data_pretreated=NULL,control_sheets=NULL) {
  if (deparse(substitute(Data)) == "control") { # control data
    P_surv_curves <-data_pretreated[data_pretreated$cond %in% control_sheets,]
  } else { # data related to treatments
    
    # Data simulation
    simulated_data <- change_data(Data)
    n_Wells <- 0
    for (i in 1:length(simulated_data)) {
      for (j in 1:length(unique(simulated_data[[i]]$num_duplic))) {
        data_tmp <-
          simulated_data[[i]][which(simulated_data[[i]]$num_duplic == unique(simulated_data[[i]]$num_duplic)[j]), ]
        n_Wells <- n_Wells + length(unique(data_tmp$rep.intra))
      }
    }
    
    # Initialisation
    Approx <- list() ; length(Approx) <- simulated_data.time_number*n_Wells
    name <- list() ; length(name) <- n_Wells
    Wells <- list() ; length(Wells) <- n_Wells*simulated_data.time_number
    Cond <- list() ; length(Cond) <- n_Wells*simulated_data.time_number
    MSE_cobs <- data.frame(MSE = rep(NA, n_Wells), Curves = rep(NA, n_Wells))
    Treat <- list() ; length(Treat) <- n_Wells*simulated_data.time_number
    
    # Survival expansion
    iter <- 1
    pb <- txtProgressBar(min = 0, max = n_Wells * simulated_data.time_number, style = 3)
    counter <- 0
    for (j in 1:length(simulated_data)) {
      for (l in 1:length(unique(simulated_data[[j]]$num_duplic))) {
        data_tmp <-
          simulated_data[[j]][which(simulated_data[[j]]$num_duplic == 
                                  unique(simulated_data[[j]]$num_duplic)[l]),]
        N <- MSpline(data_tmp)
        for (i in 1:length(unique(data_tmp$rep.intra))) {
          name[iter] <-
            paste(names(simulated_data)[[j]], "rep", i, "duplic",
                  unique(simulated_data[[j]]$num_duplic)[l], sep = "_")
          MSE_cobs$MSE[iter] <- N[simulated_data.time_number + 1, i]
          MSE_cobs$Curves[iter] <- name[[iter]]
          for (k in 1:simulated_data.time_number) {
            Approx[((iter - 1) * simulated_data.time_number) + k] <- N[k, i]
            counter <- counter + 1
            setTxtProgressBar(pb, counter)
          }
          iter <- iter + 1
        }
      }
    }
    close(pb)
    
    # Compute each well
    for (i in 1:n_Wells) {
      for (j in 1:simulated_data.time_number) {
        Wells[((i - 1) * simulated_data.time_number) + j] <- name[[i]]
      }
    }
    
    # Compute the experimental condition name
    iter <- 1
    for (i in 1:length(simulated_data)) {
      i_Treat <- NA
      for (k in 1:length(treat_group)) {
        for (l in 1:length(treat_group[[k]])) {
          if (treat_group[[k]][l] == names(simulated_data)[[i]]) {
            i_Treat <- k
          }
        }
      }
        
      for (l in 1:length(unique(simulated_data[[i]]$num_duplic))) {
        data_tmp <-
          simulated_data[[i]][which(simulated_data[[i]]$num_duplic == unique(simulated_data[[i]]$num_duplic)[l]),]
        for (k in 1:length(unique(data_tmp$rep.intra))) {
          for (j in 1:simulated_data.time_number) {
            Cond[iter] <- names(simulated_data)[[i]]
            
            if(!is.na(i_Treat)){
              Treat[iter] <- treat_names[i_Treat]  
            }else{
              Treat[iter] <- NA
            }
            iter <- iter + 1
          }
        }
      }
    }
    
    # Compute the survival curves : P_surv_curves
    P_surv_curves <- data.frame(
      time = rep(Time, n_Wells),
      surv = as.numeric(as.matrix(Approx)),
      cond = unlist(Cond),
      Wells = unlist(Wells),
      Treat = unlist(Treat)
    )
    P_surv_curves <- P_surv_curves[!is.na(P_surv_curves$Treat),]
  }
  
  # Result 
  return(P_surv_curves)
}



#### Transformation ----
deviance_transformation <- function(data_pretreated,control_pretreated){
  # Initialize
  data_mean <- data.frame(time = unique(control_pretreated$time), 
                          surv = rep(NA,length(unique(control_pretreated$time))))
  
  # Compute the average survival curve for the control group
  for (i in 1:length(data_mean$surv)){
    data_mean$surv[i] <- 
      mean(control_pretreated$surv[which(control_pretreated$time==control_pretreated$time[i])])
  }
  
  # Define the deviance curves
  data_mean2 <- data.frame(time = rep(data_mean$time,length(unique(data_pretreated$Wells))), 
                           surv = rep(data_mean$surv,length(unique(data_pretreated$Wells))))
  data_deviance <- data_pretreated
  data_deviance$surv <- data_pretreated$surv - data_mean2$surv
  
  # Return the result
  return(data_deviance)
}

# Dataframe surv par Wells
transform <- function(data){
  data2 <- matrix(nrow = length(unique(data$time)),ncol = length(unique(data$Wells)))
  for (i in 1:ncol(data2)){
    data2[,i] <- data$surv[which(data$Wells == unique(data$Wells)[i])]
  }
  return(data2)
}

#### funFEM modifications ----

funFEM_modified <- function (fd, K = 2:6, model = "AkjBk", crit = "bic", init = "kmeans", 
                             Tinit = c(), maxit = 50, eps = 1e-06, disp = FALSE, lambda = 0, 
                             graph = FALSE) 
{
  call = match.call()
  MOD = c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", 
          "AkB", "AjBk", "AjB", "ABk", "AB", "all")
  CRI = c("bic", "aic", "icl")
  INIT = c("user", "random", "kmeans", "hclust")
  if (any(!model %in% MOD)) 
    stop("Invalid model name\n", call. = FALSE)
  if (any(!crit %in% CRI)) 
    stop("Invalid criterion name\n", call. = FALSE)
  if (any(!init %in% INIT)) 
    stop("Invalid initialization name\n", call. = FALSE)
  if (any(1 %in% K)) 
    stop("K = 1 is not allowed at the moment!\n", call. = FALSE)
  if (length(model) == 1) 
    if (model == "all") 
      model = c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", 
                "AkBk", "AkB", "AjBk", "AjB", "ABk", "AB")
  resultat = vector("list", length = length(K) * length(model))
  bic = aic = icl = ll = nbprm = rep(NA, length(K) * length(model))
  it = 1
  for (k in 1:length(K)) {
    if (disp) 
      cat(">> K =", K[k], "\n")
    for (i in 1:length(model)) {
      try(resultat[[it]] <- 
            funFEM.main_modified(fd, K[k], init = init, maxit = maxit, eps = eps, 
                                 Tinit = Tinit, model = model[i],lambda = lambda, 
                                 graph = graph))
      if (length(resultat[[it]]) > 0) {
        try(bic[it] <- resultat[[it]]$bic)
        try(aic[it] <- resultat[[it]]$aic)
        try(icl[it] <- resultat[[it]]$icl)
        try(nbprm[it] <- resultat[[it]]$nbprm)
        try(ll[it] <- resultat[[it]]$ll)
        if (disp) {
          if (crit == "bic") 
            cat(model[i], "\t:\t bic =", resultat[[it]]$bic, 
                "\n")
          if (crit == "aic") 
            cat(model[i], "\t:\t aic =", resultat[[it]]$aic, 
                "\n")
          if (crit == "icl") 
            cat(model[i], "\t:\t icl =", resultat[[it]]$icl, 
                "\n")
        }
      }
      it = it + 1
    }
  }
  if (prod(is.na(bic)) == 1) 
    stop("No reliable results to return (empty clusters in all partitions)!")
  if (crit == "bic") {
    id_max = which.max(bic)
    crit_max = resultat[[id_max]]$bic
  }
  if (crit == "aic") {
    id_max = which.max(aic)
    crit_max = resultat[[id_max]]$aic
  }
  if (crit == "icl") {
    id_max = which.max(icl)
    crit_max = resultat[[id_max]]$icl
  }
  res = resultat[[id_max]]
  if (disp) 
    cat("The best model is", res$model, "with K =", res$K, 
        "(", crit, "=", crit_max, ")\n")
  res$crit = crit
  nm = length(model)
  res$allCriterions = data.frame(K, model, bic = matrix(bic, 
                                                        ncol = nm, byrow = T), aic = matrix(aic, ncol = nm, byrow = T), 
                                 icl = matrix(icl, ncol = nm, byrow = T), nbprm = matrix(nbprm, 
                                                                                         ncol = nm, byrow = T), ll = matrix(ll, ncol = nm, 
                                                                                                                            byrow = T))
  res$call = call
  class(res) = "fem"
  res
}
.estep_modified <- function (prms, fd, U) {
  threshold_error <- 1e-6
  Y = t(fd$coefs)
  n = nrow(Y)
  p = ncol(Y)
  K = prms$K
  mu = prms$mean
  prop = prms$prop
  D = prms$D
  d = K - 1
  QQ = matrix(NA, n, K)
  QQ2 = matrix(NA, n, K)
  T = matrix(NA, n, K)
  for (k in 1:K) {
    bk = D[k, p, p]
    mY = prms$my[k, ]
    YY = Y - t(matrix(rep(mY, n), p, n))
    projYY = YY %*% U %*% t(U)
    if (d == 1) {
      for (i in 1:n) {
        QQ[i, k] = 1/D[k, 1, 1] * sum(projYY[i, ]^2) + 
          1/D[k, p, p] * sum((YY[i, ] - projYY[i, ])^2) + 
          (p - d) * log(bk) + log(D[k, 1, 1]) - 2 * log(prop[k]) + 
          p * log(2 * pi)
      }
    } else {
      sY = U %*% ginv(D[k, (1:d), (1:d)]) %*% t(U)
      for (i in 1:n) {
        det_tmp <- det(D[k, (1:d),(1:d)]) # Avoid production of NaN
        if(det_tmp < threshold_error){
          det_tmp <- threshold_error
        }
        QQ[i, k] = projYY[i, ] %*% sY %*% as.matrix(projYY[i, 
        ], p, 1) + 1/bk * sum((YY[i, ] - projYY[i, 
        ])^2) + (p - d) * log(bk) + log(det_tmp) - 2 * log(prop[k]) + p * log(2 * pi)
      }
    }
  }
  A = -1/2 * QQ
  loglik = sum(log(rowSums(exp(A - apply(A, 1, max)))) + apply(A, 
                                                               1, max))
  for (k in 1:K) {
    T[, k] = 1/rowSums(exp(0.5 * (QQ[, k] * matrix(1, n, 
                                                   K) - QQ)))
  }
  list(T = T, loglik = loglik)
}
funFEM.main_modified <- 
  function (fd, K, model = "AkjBk", init = "kmeans", lambda = 0, Tinit = c(), 
            maxit = 50, eps = 1e-08, graph = F) 
  {
    Y = t(fd$coefs)
    n = nrow(Y)
    p = ncol(Y)
    Lobs = rep(c(-Inf), 1, (maxit + 1))
    if (init == "user") {
      T = Tinit
    } else if (init == "kmeans") {
      
      cond <- TRUE
      test_iter <- 1
      max_test_iter <- 100
      while(cond & test_iter <= max_test_iter){ # rajout d'une condition qui empeche de ne pas démarrer
      T = matrix(0, n, K)
      ind = kmeans(Y, K)$cluster
      for (k in 1:K) {
        T[which(ind == k), k] = 1
      }
       cond <- min(colSums(T)) <= 1
       test_iter <- test_iter + 1
      }
    } else if (init == "random") {
      T = t(rmultinom(n, 1, c(rep(1/K, K))))
    } else if (init == "hclust") {
      T = matrix(0, n, K)
      ind = cutree(hclust(dist(Y), method = "ward.D2"), K)
      for (k in 1:K) {
        T[which(ind == k), k] = 1
      }
    }
    V = funFEM:::.fstep(fd, T, lambda)
    prms = funFEM:::.mstep(fd, V, T, model = model)
    res.estep = .estep_modified(prms, fd, V)
    T = res.estep$T
    Lobs[1] = res.estep$loglik
    Linf_new = Lobs[1]
    for (i in 1:maxit) {
      V = funFEM:::.fstep(fd, T, lambda)
      prms = funFEM:::.mstep(fd, V, T, model = model)
      res.estep = .estep_modified(prms, fd, V)
      T = res.estep$T
      Lobs[i + 1] = res.estep$loglik
      if (i >= 2) {
        acc = (Lobs[i + 1] - Lobs[i])/(Lobs[i] - Lobs[i - 
                                                        1])
        Linf_old = Linf_new
        Linf_new = Lobs[i] + 1/(1 - acc) * (Lobs[i + 1] - 
                                              Lobs[i])
        if (abs(Linf_new - Linf_old) < eps | is.na(Linf_new)) {
          break
        }
      }
    }
    if (graph) {
      par(mfrow = c(1, 2))
      plot(as.data.frame(as.matrix(Y) %*% V[, 1:2]), col = max.col(T), 
           xlab = "axis 1", ylab = "axis 2", pch = 20)
      plot(Lobs[1:i], xlab = "iterations", ylab = "vraisemblance Espace observe", 
           col = 2, pch = 20)
    }
    cls = max.col(T)
    crit = funFEM:::.criteria(Lobs[(i + 1)], T, prms, n)
    W = inprod(fd$basis, fd$basis)
    U = t(W) %*% V
    res = list(model = model, K = K, cls = cls, P = T, prms = prms, 
               U = U, aic = crit$aic, bic = crit$bic, icl = crit$icl, 
               loglik = Lobs[2:(i + 1)], ll = Lobs[i + 1], nbprm = crit$nbprm,V=V)
  }

#### Do_clustering ----

funFEM_ICL <- function(data){
  #  Change the data structure
  data_without_simulation <- 
    data[grepl(x=data$Wells,pattern = "observed"),]
  data2 <- transform(data_without_simulation)
  
  # Bsplines expansion for the survival curbes
  Base <- create.bspline.basis(rangeval = c(0,time_max),
                               nbasis = n_functional_base)
  fd <- smooth.basis(Time,data2,Base)$fd
  Approx <- eval.fd(Time,fd)
  
  # Initialiaze object for the funFEM fitting
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk",
              "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  ICLs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  
  # Fit all the funFEM models 
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      res_try <- try(try_tmp <- funFEM_modified(fd,K=min_K:max_K,model=models[j],
                                       init='kmeans',lambda=0,disp=FALSE,crit = "icl"))
      if(class(res_try) == "fem"){
        model_choice_tmp <- try_tmp
        Ks[i,j] <- model_choice_tmp$K
        ICLs[i,j] <- model_choice_tmp$icl
      }
      
      # model_choice_tmp <- 
      #   suppressWarnings(funFEM_modified(fd,K=min_K:max_K,model=models[j],
      #                                    init='kmeans',lambda=0,disp=FALSE,crit = "icl"))
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  
  # Model choice
  choice_K <- apply(Ks,2,table)
  ICL <- apply(ICLs,2,mean)
  m <- models[which.max(ICL)]
  K <- as.numeric(names(which.max(choice_K[[which.max(ICL)]])))
  
  # Fit the best model
  fitted_model <- funFEM_modified(fd,K=K,model=m,init='kmeans',lambda=0,
                                  disp=FALSE,crit = "icl")
  
  # Cluster allocation 
  data_without_simulation$cls <- rep(fitted_model$cls,each=simulated_data.time_number)
  simulated_data <- 
    data[!(grepl(x=data$Wells,pattern = "observed")),]
  data2 <- transform(simulated_data)
  simulated_fd <- smooth.basis(Time,data2,Base)$fd
  res.estep <- .estep_modified(fitted_model$prms, simulated_fd, fitted_model$V)
  P <- res.estep$T
  simulated_data$cls = rep(max.col(P),each=simulated_data.time_number)
  
  # Put cluster allocation into the original data object
  data$Cls <- rep(NA,nrow(data))
  data$status <- rep(NA,nrow(data))
  
  Wells <- unique(data$Wells)
  for(well in Wells){
    if(grepl(x=well,pattern = "observed")){
      data$Cls[data$Wells == well] <- 
        data_without_simulation$cls[data_without_simulation$Wells == well]
      data$status[data$Wells == well] <- "observed"
    }else{
      data$Cls[data$Wells == well] <- 
        simulated_data$cls[simulated_data$Wells == well]
      data$status[data$Wells == well] <- "simulated"
    }
  }
  data$Cls <- factor(data$Cls)
  
  # Information to save
  P_final <- matrix(NA,ncol=ncol(P),nrow=nrow(P)+nrow(fitted_model$P))
  P_final[grepl(x=unique(data$Wells),pattern = "observed"),] <- fitted_model$P
  P_final[!(grepl(x=unique(data$Wells),pattern = "observed")),] <- P
  
  fitted_model$P_final <-  P_final
  
  save(K,Approx,fitted_model,file = Clustering_path)
  
  # Result to return
  return(data)
}

#### Compute_proba ----
Pb_funFEM <- function(Data){
  # Initialize ad load required object
  load(Clustering_path)
  threshold_deviance <- 0.01
  threshold_deviance2 <- 0.15
  
  # Define Cond
  my_Wells <- unique(Data$Wells)
  Cond <- list() ; length(Cond) <- length(my_Wells)
  for (i in 1:length(my_Wells)) {
    condition <- Data$cond[which(Data$Wells == my_Wells[i])][1]
    test <- strsplit(as.character(condition),split = "\n")[1]
    Cond[i] <- test[[1]]
  }
  
  # Compute the proba
  Proba_curves <- data.frame(fitted_model$P_final,cond = unlist(Cond))
  Proba_curves_bis <- data.frame(fitted_model$P_final)
  
  n_Treat <- length(treat_group)
  Proba_Treat <- data.frame(Pb = rep(NA,n_Treat*K), Cls = rep(NA,n_Treat*K), Treat = rep(NA,n_Treat*K))
  
  # For the first treatmeant
  Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[1]],K))
  for (j in 1:K){
    dat <- Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[1]),]
    if (length(treat_group[[1]])>1) {
      for (k in 2:length(unlist(treat_group[[1]]))) {
        dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[k]),])
      }
    }
    Proba$Pb[j] <- mean(dat[,j]) # 1/(dim(dat)[1])*sum(dat[,j])
    Proba$Cls[j] <- j
    Proba$Treat[j] <- treat_names[[1]]
  }
  Proba_Treat <- Proba
  
  # For the other ones
  for (i in 2:n_Treat){
    Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[i]],K))
    for (j in 1:K) {
      dat <- Proba_curves[which(Proba_curves$cond==treat_group[[i]][1]),]
      if (length(unlist(treat_group[[i]]))>1) {
        for (k in 2:length(unlist(treat_group[[i]]))) {
          dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[i]])[k]),])
        }
      }
      Proba$Pb[j] <- (1/(dim(dat)[1]))*sum(dat[,j])
      Proba$Treat[j] <- treat_names[[i]]
      Proba$Cls[j] <- j
    }
    Proba_Treat <- rbind(Proba_Treat,Proba)
  }
  
  # Cluster names
  name_cls <- list() ; length(name_cls) <- K
  sign <- matrix(NA, nrow=K, ncol = n_cut)
  
  means <- matrix(NA,ncol=n_cut,nrow=K)
  for (i in 1:K) {
    d <- Data[which(Data$Cls==i),]
    for (j in 1:n_cut) {
      t_m <- ((j-1)*time_max)/n_cut
      t_p <- ((j*time_max)/n_cut)
      means[i,j] <- mean(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
    }
  }
  
  number_sign <- floor(means / threshold_deviance2)
  
  symbol_sign <- sign(number_sign)
  symbol_sign[symbol_sign == -1] <- "-"
  symbol_sign[symbol_sign ==  0] <- "="
  symbol_sign[symbol_sign ==  1] <- "+"
  
  name_cls <- rep(NA,K)
  for(i in 1:K){
    name_i <- NULL
    for(j in 1:n_cut){
      if(number_sign[i,j] != 0){
        name_i <- paste(name_i,
                        paste0(rep(symbol_sign[i,j],each=abs(number_sign[i,j])),collapse = ""),
                        sep="/")
      } else{
        name_i <- paste(name_i,"=",sep="/")
      }
    }
    name_i <- gsub(x = name_i,pattern = "^/",replacement = "")
    name_cls[i] <-  name_i
  } 
  
  # #### PMG 2023-05-23 : rajout pour faire en sorte que le groupe témoin soit toujjours le premier
  control_group <- which.max(Proba_Treat$Pb[Proba_Treat$Treat %in% c("Témoin","Temoin","Control")])
  name_cls[control_group] <- "Control"
  name_cls[-control_group] <- paste("G",1:(K-length(control_group))," ",name_cls[-control_group],sep = "")
  # #### Fin rajout
  
  Proba_Treat$Cls <- unlist(name_cls)
  
  # Allocate the cluster names
  new_names <- list() ; length(new_names) <- nrow(data_cluster)
  for (i in 1:K) {
    new_names[which(data_cluster$Cls == i)] <- name_cls[[i]]
  }
  data4 <- data_cluster
  data4$Cls <- unlist(new_names)
  
  # Compute the average deviance for each cluster
  mean_deviance <- data.frame(time = rep(Time,K),
                           surv = rep(NA,K*simulated_data.time_number),
                           Cls = factor(rep(unique(data4$Cls),
                                            each = simulated_data.time_number)))
  my_mean <- list() ; length(my_mean) <- simulated_data.time_number*K
  for (i in 1:length(unique(data4$Cls))) {
    data_tmp <- data4[which(data4$Cls == unique(data4$Cls)[i]),]
    for (j in 1:simulated_data.time_number) {
      my_mean[(i-1)*simulated_data.time_number+j] <- mean(data_tmp$surv[which(data_tmp$time == Time[j])])
    }
  }
  mean_deviance$surv <- unlist(my_mean)
  
  # Save results
  path <- paste0(results_path,"/Proba_Treat.RData")
  save(Proba_Treat, file = path)
  save(Proba_curves, Proba_curves_bis, Proba_Treat, mean_deviance, file = Proba_path)
  
  # Return result
  return(data4)
}

#### Graphical results ----

compute_graphical_output <- function(){
  
  # Initialize
  color_function <- rainbow
  load(Clustering_path)
  load(Proba_path)
  
  # Pretraeat data
  data_pretreated$Treat[data_pretreated$Treat %in% c("Temoin","Témoin")] <- "Control"
  
  final_data$Treat[final_data$Treat %in% c("Temoin","Témoin")] <- "Control"
  final_data$Cls[final_data$Cls %in% c("Temoin","Témoin")] <- "Control"
  levels_cls <- unique(final_data$Cls)
  levels_cls <- c("Control",sort(levels_cls[levels_cls != "Control"]))
  final_data$Cls <- factor(final_data$Cls,levels = levels_cls)
                           
  data_transformed$Treat[data_transformed$Treat %in% c("Temoin","Témoin")] <- "Control"
  levels(mean_deviance$Cls)[match("Témoin",levels(mean_deviance$Cls))] <- "Control"
  
  periode_sep <- (1:(n_cut-1))*time_max/n_cut
  
  status_observed <- as.character(grepl(x = data_pretreated$Wells,pattern = "observed"))
  status_observed[status_observed == "TRUE"] <- "observed"
  status_observed[status_observed == "FALSE"] <- "simulated"
  
  Treat_levels <- unique(data_pretreated$Treat)
  Treat_levels <- c("Control",Treat_levels[Treat_levels != "Control"])
  data_pretreated$Treat <- factor(data_pretreated$Treat,
                                levels = Treat_levels)
  data_transformed$Treat <- factor(data_transformed$Treat,
                                levels = Treat_levels)
  
  
  # Graphical outputs : survival curves
  n <- length(unique(data_transformed$Treat))
  p <- ggplot(data_pretreated,aes(time,surv,col=Treat, group = Wells,
                                  alpha=status_observed,size=status_observed)) + 
    scale_size_manual(values=c(0.8,0.3)) + 
    scale_alpha_manual(values=c(1,alpha_simulated_data)) + 
    geom_line()+ 
    theme_pubclean(base_size = 15) +
    xlab("Days") + ylab("Survival probability") + 
    ggtitle("Survival curves") + 
    scale_color_manual(values = c("gray",color_function(n-1))) + 
    labs(col="",alpha="",size="") 
  Path <- paste0(img_path,"/Survival_curves.pdf")
  ggsave(Path,plot = p, width = 8,height = 8)
  
  # Graphical outputs : survival curves (only observed data)
  p <- ggplot(data_pretreated[status_observed == "observed",],aes(time,surv,col=Treat, group = Wells)) + 
    geom_line()+ 
    theme_pubclean(base_size = 15) +
    xlab("Days") + ylab("Survival probability") + 
    ggtitle("Survival curves") + 
    scale_color_manual(values = c("gray",color_function(n-1))) + 
    labs(col="",alpha="",size="") 
  Path <- paste0(img_path,"/Survival_curves_observed.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  # Graphical outputs : average survival curves according to treatment
  p <- ggplot(data_pretreated) +
    geom_line(aes(time,surv,col=Treat,group = Wells),
              size = 0.2, alpha = alpha_simulated_data) +
    xlab("Days") + ylab("Survival probability") +
    ggtitle("Survival probability and treatments") +
    scale_color_manual(values = c("gray",color_function(n-1))) + labs(col="") +
    geom_smooth(data = data_pretreated, aes(time, surv, group = Treat),
                span = 0.3,col="black",size=1.8*1.5)+
    geom_smooth(data = data_pretreated, aes(time, surv, col = Treat), 
                span = 0.3,size=1.8)+
    theme_pubclean(base_size = 15)
  Path <- paste0(img_path,"/Survival_mean_curves.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  
  # Graphical outputs : survival deviance curves
  p <- ggplot(data_transformed,aes(time,surv,col=Treat, group = Wells,
                                  alpha=status_observed,size=status_observed)) + 
    scale_size_manual(values=c(0.8,0.3)) + 
    scale_alpha_manual(values=c(1,alpha_simulated_data)) + 
    geom_line()+ 
    xlab("Days") + 
    ylab("Survival deviance") + 
    ggtitle("Survival deviance") + 
    scale_color_manual(values = c("gray",color_function(n-1))) +
    theme_pubclean(base_size = 15) +
    labs(col="",alpha="",size="") 
  Path <- paste0(img_path,"/Survival_deviance.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  # Graphical outputs : average survival deviance curves according to treatment
  p <- ggplot(data_transformed) +
    geom_line(aes(time,surv,col=Treat,group = Wells),
              size = 0.2, alpha = alpha_simulated_data) +
    xlab("Days") + ylab("Survival deviance") +
    ggtitle("Survival deviance and treatments") +
    scale_color_manual(values = c("gray",color_function(n-1))) + labs(col="") +
    geom_smooth(data = data_transformed, aes(time, surv, group = Treat),
                span = 0.3,col="black",size=1.8*1.5)+
    geom_smooth(data = data_transformed, aes(time, surv, col = Treat), 
                span = 0.3,size=1.8)+
    theme_pubclean(base_size = 15)
  Path <- paste0(img_path,"/Survival_deviance_mean.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  # Graphical outputs : average survival curves according to estimated clusters
  data_pretreated$Cls <- final_data$Cls
  data_pretreated$Cls <- factor(data_pretreated$Cls,levels=levels_cls)
  p <- ggplot(data_pretreated,aes(time,surv,col=Cls, group = Wells)) +
    geom_vline(xintercept = periode_sep,linetype="dashed",col="gray")+ 
    geom_line(size=0.3,alpha=alpha_simulated_data) +
    xlab("Days") + ylab("Survival probability") + 
    ggtitle("Survival clusters") + 
    scale_color_manual(values = c("gray",color_function(K-1)))+
    labs(col="")+ 
    geom_smooth(data = data_pretreated, aes(time, surv, group = Cls),
                span = 0.3,col="black",size=1.8*1.5)+
    geom_smooth(data = data_pretreated, aes(time, surv, group = Cls),
                span = 0.3,size=1.8)+
    theme_pubclean(base_size = 15) 
  Path <- paste0(img_path,"/Survival_clusters.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  # Graphical outputs : average survival deviance curves according to estimated clusters
  p <- ggplot(final_data, aes(time,surv,col=Cls,group = Wells)) + 
    geom_vline(xintercept = periode_sep,linetype="dashed",col="gray")+ 
    geom_line(size=0.3,alpha=alpha_simulated_data) +
    xlab("Days") + ylab("Survival deviance") + 
    ggtitle("Survival differences and clusters") + 
    scale_color_manual(values = c("gray",color_function(K-1))) +
    geom_smooth(data = final_data, aes(time, surv, group = Cls),
                span = 0.3,col="black",size=1.8*1.5)+
    geom_smooth(data = final_data, aes(time, surv, group = Cls),
                span = 0.3,size=1.8)+
    labs(col="")+ 
    theme_pubclean(base_size = 15)
  Path <- paste0(img_path,"/Survival_deviance_clusters.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",width = 8,height = 8))
  
  # Graphical outputs : Probability outputs 
  Proba_plot_path <- paste0(img_path,"/Cluster_distribution.pdf")
  
  Proba_Treat$Treat[Proba_Treat$Treat %in% c("Temoin","Témoin")] <- "Control"
  Proba_Treat$Cls[Proba_Treat$Cls %in% c("Temoin","Témoin")] <- "Control"
  Proba_Treat$Cls <- factor(Proba_Treat$Cls,levels=levels_cls)
  Proba_Treat$Treat <- factor(Proba_Treat$Treat,
                                 levels = Treat_levels)
  
  data.separe <- data.frame(x=seq(1.5, K - 0.5, length.out = (K-1)),
                            y=rep(0,(K-1)),xend=seq(1.5, K - 0.5, length.out = (K-1)),
                            yend=rep(1,(K-1)))
  
  p <- ggplot(Proba_Treat,aes(x=Cls,y=Pb,fill=Treat)) + geom_bar(stat="identity",position=position_dodge()) + 
    scale_y_continuous(limits = c(0,1)) + xlab("Clusters") + ylab("Frequency") +
    ggtitle("Experimental conditions and clusters")+
    scale_fill_manual(values=c("gray",color_function(length(table(Proba_Treat$Treat))-1))) + 
    labs(fill = "") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    geom_segment(data=data.separe, aes(x=x, y=y, xend=xend, yend=yend,fill=NULL),linetype="dashed", alpha=5)+ 
    theme_pubclean(base_size = 15)
  Path <- paste0(img_path,"/Association_Conditions_Clusters.pdf")
  suppressMessages(ggsave(Path,plot = p, device = "pdf",height = 5,width=10))
  
  Cls <- unique(Proba_Treat$Cls)
  for(cls in Cls){
    cls_name <- strsplit(cls,split=" ")[[1]][1]
    p <- ggplot(Proba_Treat[Proba_Treat$Cls == cls,]) + aes(x="",y=Pb,fill=Treat) + 
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      scale_fill_manual(values=c("gray",color_function(length(table(Proba_Treat$Treat))-1))) 
    Path <- paste0(img_path,"/Cluster_pie_",cls_name,".png")
    suppressMessages(ggsave(Path,plot = p, device = "png",height = 5,width=10))
  }
}