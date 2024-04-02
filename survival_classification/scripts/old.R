################################################################################
################ A draft containing functions that are no longer in use for 
################   "Survival curve classification"
######## Paul-Marie Grollemund
####### 2024-04-02 : first implementation
################################################################################

################################################################################
################################################################################
# functions_choice.R
################################################################################
################################################################################
# parametrisation ----
if (type_param == "Weibull") Nls <- NLS_wei 
if (type_param == "Lognormal") Nls <- NLS_LogN

################################################################################
################################################################################
# functions.R
################################################################################
################################################################################

#### Pretreatment ----
# Parametrisation avec Weibull, retourne p et lambda
Wei <- function(Data){
  W <- survival::survreg(Surv(Temps,Statut) ~ rep.intra, data = Data, dist = "weibull")
  lambda <- exp(W$icoef[1])
  p <- 1/W$scale
  Params <- c(p,lambda)
  return(Params)
}

# Parametrisation avec la cdf LogNormale, retourne mu et sd
LogN <- function(Data){
  Ln <- survival::survreg(Surv(Temps,Statut) ~ rep.intra, data = Data, dist = "lognormal")
  mu <- Ln$coef[1]
  sd <- Ln$scale
  rep2 <- Ln$coefficients[2]
  rep3 <- Ln$coefficients[3]
  Params <- c(mu,sd,rep2,rep3)
  return(Params)
}
# Renvoie 1 - fonction de repartition d'une distribution lognormale
FLogN <- function(X,mu,sd) {
  return(1-plnorm(X, meanlog = mu, sdlog = sd))
}
# Prediction sur Newdata (Weibull)
NLS_wei <- function(Data){
  pred <- matrix(NA, nrow = length_newdata, ncol = length(unique(Data$rep.intra)))
  for (i in 1:length(unique(Data$rep.intra))){
    d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$time,
                    y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$surv)
    nls <- nls(y~exp(-lambda^{-p} * x^p), data=d,
               start = list(lambda = Wei(Data)[2], p = Wei(Data)[1]))
    pred[,i] <- predict(nls,newdata=data.frame(x = Time))
  }
  return(pred)
}

# Prediction sur Newdata (Lognormal)
NLS_LogN <- function(Data){
  pred <- matrix(NA, nrow = length_newdata, ncol = length(unique(Data$rep.intra)))
  for (i in 1:length(unique(Data$rep.intra))){
    d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$time,
                    y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$surv)
    print(head(d))
    print(FLogN(d$x,LogN(Data)[1],LogN(Data)[2]))
    Coefs <- LogN(Data)
    if (i == 1) {
      nls <- nls(y~FLogN(x,Mu,Sd), data=d,
                 start = list(Mu = Coefs[1], Sd = Coefs[2]))
    }
    if (i == 2) {
      nls <- nls(y~FLogN(x,Mu,Sd), data=d,
                 start = list(Mu = Coefs[1]+Coefs[3], Sd = Coefs[2]))
      # while(class(nls) != "nls") {
      #   nls <- try(nls(y~FLogN(x,Mu,Sd), data=d,
      #                  start = list(Mu = Coefs[1]+Coefs[3], Sd = Coefs[2])),quietly=TRUE)
      # }
    }
    if (i == 3) {
      nls <- nls(y~FLogN(x,Mu,Sd), data=d,
                 start = list(Mu = Coefs[1], Sd = Coefs[2]))
    }
    pred[,i] <- predict(nls,newdata=data.frame(x = Time))
  }
  return(pred)
}


## Conspline
ConSpline <- function(Data) {
  suppressPackageStartupMessages(library(MLmetrics))
  pred <- NULL 
  for (i in 1:length(unique(Data$rep.intra))){
    d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$time,
                    y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$surv)
    ## Interpolation
    options(show.error.messages = FALSE)
    Cspline <- ConSpline::conspline(d$y, d$x, type = 2)
    options(show.error.messages = TRUE)
    # Ajustement
    Pred <- data.frame(time = d$x, surv = Cspline$fhat, Wells = rep(i,nrow(d)))
    Pred$surv[which(Pred$surv > 1)] <- 1
    Pred$surv[which(Pred$surv < 0)] <- 0
    Pred$MSE <- rep(MSE(Pred$surv,d$y),nrow(Pred))
    pred <- rbind(pred,Pred)
  }
  return(pred)
}

## Kaplan-Meier
K_m <- function(Data) {
  pred <- NULL
  for (i in 1:length(unique(Data$rep.intra))){
    d <- data.frame(x = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$time,
                    y = fit_plot(Data)$data.survplot[which(fit_plot(Data)$data.survplot$rep.intra==unique(Data$rep.intra)[i]),]$surv)
    # Ajustement
    Pred <- data.frame(time = d$x, surv = d$y, Wells = rep(i,nrow(d)))
    pred <- rbind(pred,Pred)
  }
  return(pred)
}

#### Pretreat conspline 
pretreat_conspline <- function(Data) {
  # Determination des parametres et de la prediction de la survie
  # cat("- Détermination des paramètres de chaque fonction de survie. \n")
  P_surv_curves <- NULL
  for (j in 1:length(Data)) {
    N <- ConSpline(Data[[j]])
    for (i in unique(N$Wells)) {
      N$Wells[which(N$Wells == i)] <- paste0(names(Data)[[j]],"_",i)
    }
    N$cond <- rep(names(Data)[[j]],nrow(N))
    for (k in 1:length(treat_group)) {
      for (l in 1:length(treat_group[[k]])){
        if (treat_group[[k]][l]==names(Data)[[j]]) i_Treat <- k
      }
    }
    for (k in 1:length(unique(N$Wells))) {
      N$Treat <- rep(treat_names[i_Treat],nrow(N))
    }
    P_surv_curves <- rbind(P_surv_curves,N)
  }
  lenlen <- length(unique(P_surv_curves$Wells))
  row_index <- list() ; length(row_index) <- lenlen
  MSE_cons <- list() ; length(MSE_cons) <- lenlen
  curves <- list() ; length(curves) <- lenlen
  for (i in 1:length(row_index)) {
    MSE_cons[i] <- unique(P_surv_curves$MSE[which(P_surv_curves$Wells == unique(P_surv_curves$Wells)[i])])
    curves[i] <- unique(P_surv_curves$Wells)[i]
  }
  MSE_conspline <- data.frame(Curves = unlist(curves), MSE = unlist(MSE_cons))
  # Path <- paste0(results_path,"/MSE_conspline.RData")
  # save(MSE_conspline, file = Path)
  P_surv_curves$MSE <- NULL
  return(P_surv_curves)
}

#### Pretreat Splines
pretreat_splines <- function(Data) {
  # cat("- Chargement des packages nécessaires. \n")
  suppressPackageStartupMessages(library(MLmetrics))
  suppressPackageStartupMessages(library(utils))
  #Determination de la taille des listes a creer
  # cat("- Calcul de la taille des listes à créer. \n")
  n_Wells <- 0
  for (i in 1:length(Data)) {
    n_Wells <- n_Wells + length(unique(Data[[i]]$rep.intra))
  }
  # Initialisation
  # cat("- Création des listes. \n")
  Approx <- list() ; length(Approx) <- simulated_data.time_number*n_Wells
  name <- list() ; length(name) <- n_Wells
  Wells <- list() ; length(Wells) <- n_Wells*simulated_data.time_number
  Cond <- list() ; length(Cond) <- n_Wells*simulated_data.time_number
  MSE_cobs <- data.frame(MSE = rep(NA, n_Wells), Curves = rep(NA, n_Wells))
  if (deparse(substitute(Data)) == "data") {
    Treat <- list() ; length(Treat) <- n_Wells*simulated_data.time_number
  }
  # Determination des parametres et de la prediction de la survie
  # cat("- Détermination des paramètres de chaque fonction de survie. \n")
  iter <- 1
  pb <- txtProgressBar(min = 0, max = n_Wells*simulated_data.time_number, style = 3)
  counter <- 0
  for (j in 1:length(Data)) {
    N <- MSpline(Data[[j]])
    for (i in 1:length(unique(Data[[j]]$rep.intra))) {
      name[iter] <- paste0(names(Data)[[j]],"_",i)
      MSE_cobs$MSE[iter] <- N[simulated_data.time_number+1,i]
      MSE_cobs$Curves[iter] <- name[[iter]]
      for (k in 1:simulated_data.time_number) {
        Approx[((iter-1)*simulated_data.time_number)+k] <- N[k,i]
        # update progress bar
        counter <- counter+1
        setTxtProgressBar(pb,counter)
      }
      iter <- iter + 1
    }
  }
  close(pb)
  ## P_surv_curves
  # cat("- Détermination de la variable Wells. \n")
  for (i in 1:n_Wells) {
    for (j in 1:simulated_data.time_number) {
      Wells[((i-1)*simulated_data.time_number)+j] <- name[[i]]
    }
  }
  # cat("- Détermination de la variable Cond. \n")
  if (deparse(substitute(Data)) == "data") {
    iter <- 1
    for (i in 1:length(Data)) {
      for (k in 1:length(treat_group)) {
        for (l in 1:length(treat_group[[k]])){
          if (treat_group[[k]][l]==names(Data)[[i]]) i_Treat <- k
        }
      }
      for (k in 1:length(unique(Data[[i]]$rep.intra))) {
        for (j in 1:simulated_data.time_number) {
          Cond[iter] <- names(Data)[[i]]
          Treat[iter] <- treat_names[i_Treat]
          iter <- iter + 1
        }
      }
    }
  } else {
    iter <- 1
    for (i in 1:length(Data)) {
      for (k in 1:length(unique(Data[[i]]$rep.intra))) {
        for (j in 1:simulated_data.time_number) {
          Cond[iter] <- names(Data)[[i]]
          iter <- iter + 1
        }
      }
    }
  }
  # cat("- Création de P_surv_curves. \n")
  if (deparse(substitute(Data)) != "data") {
    P_surv_curves <- data.frame(time = rep(Time,n_Wells),
                                surv = as.numeric(as.matrix(Approx)),
                                cond = unlist(Cond),
                                Wells = unlist(Wells))
  } else {
    P_surv_curves <- data.frame(time = rep(Time,n_Wells),
                                surv = as.numeric(as.matrix(Approx)),
                                cond = unlist(Cond),
                                Wells = unlist(Wells),
                                Treat = unlist(Treat))
    Path <- paste0(results_path,"/MSE_cobs.RData")
    save(MSE_cobs,file = Path)
  }
  return(P_surv_curves)
}

#### Pretreat Kaplan-Meier
pretreat_KM <- function(Data) {
  # Determination des parametres et de la prediction de la survie
  # cat("- Détermination des paramètres de chaque fonction de survie. \n")
  P_surv_curves <- NULL
  for (j in 1:length(Data)) {
    N <- K_m(Data[[j]])
    for (i in unique(N$Wells)) {
      N$Wells[which(N$Wells == i)] <- paste0(names(Data)[[j]],"_",i)
    }
    N$cond <- rep(names(Data)[[j]],nrow(N))
    for (k in 1:length(treat_group)) {
      for (l in 1:length(treat_group[[k]])){
        if (treat_group[[k]][l]==names(Data)[[j]]) i_Treat <- k
      }
    }
    for (k in 1:length(unique(N$Wells))) {
      N$Treat <- rep(treat_names[i_Treat],nrow(N))
    }
    P_surv_curves <- rbind(P_surv_curves,N)
  }
  return(P_surv_curves)
}

#### Pretraitement NLS
pretreat_nls <- function(Data) {
  #Determination de la taille des listes a creer
  # cat("- Calcul de la taille des listes à créer. \n")
  n_puits <- 0
  for (i in 1:length(Data)) {
    n_puits <- n_puits + length(unique(Data[[i]]$rep.intra))
  }
  # Initialisation
  # cat("- Création des listes. \n")
  NLS <- list() ; length(NLS) <- simulated_data.time_number*n_puits
  name <- list() ; length(name) <- n_puits
  Puits <- list() ; length(Puits) <- n_puits*simulated_data.time_number
  Cond <- list() ; length(Cond) <- n_puits*simulated_data.time_number
  if (deparse(substitute(Data)) != "control") {
    Mol <- list() ; length(Mol) <- n_puits*simulated_data.time_number}
  # Determination des parametres et de la prediction de la survie
  # cat("- Détermination des paramètres de chaque fonction de survie. \n")
  iter <- 1
  for (j in 1:length(Data)) {
    N <- Nls(Data[[j]])
    for (i in 1:length(unique(Data[[j]]$rep.intra))) {
      name[iter] <- paste0(names(Data)[[j]],"_",i)
      for (k in 1:simulated_data.time_number) {
        NLS[((iter-1)*simulated_data.time_number)+k] <- N[k,i]
      }
      iter <- iter + 1
    }
  }
  ## P_surv_curves
  # cat("- Détermination de la variable Puits. \n")
  for (i in 1:n_puits) {
    for (j in 1:simulated_data.time_number) {
      Puits[((i-1)*simulated_data.time_number)+j] <- name[[i]]
    }
  }
  # cat("- Détermination de la variable Cond. \n")
  if (deparse(substitute(Data)) != "control") {
    iter <- 1
    for (i in 1:length(Data)) {
      for (k in 1:length(treat_group)) {
        for (l in 1:length(treat_group[[k]])){
          if (treat_group[[k]][l]==names(Data)[[i]]) i_mol <- k
        }
      }
      for (k in 1:length(unique(Data[[i]]$rep.intra))) {
        for (j in 1:simulated_data.time_number) {
          Cond[iter] <- names(Data)[[i]]
          Mol[iter] <- treat_names[i_mol]
          iter <- iter + 1
        }
      }
    }
  } else {
    iter <- 1
    for (i in 1:length(Data)) {
      for (k in 1:length(unique(Data[[i]]$rep.intra))) {
        for (j in 1:simulated_data.time_number) {
          Cond[iter] <- names(Data)[[i]]
          iter <- iter + 1
        }
      }
    }
  }
  # cat("- Création de P_surv_curves. \n")
  if (deparse(substitute(Data)) == "control") {
    P_surv_curves <- data.frame(time = rep(Time,n_puits),
                                surv = as.numeric(as.matrix(NLS)),
                                cond = unlist(Cond),
                                puits = unlist(Puits))
  } else {
    P_surv_curves <- data.frame(time = rep(Time,n_puits),
                                surv = as.numeric(as.matrix(NLS)),
                                cond = unlist(Cond),
                                puits = unlist(Puits),
                                mol = unlist(Mol))
  }
  
  return(P_surv_curves)
}

#### Test simulation ----
simul_data_old <- function(data, n) {
  if (n == 0) {
    return(data)
  }
  load(sd_error_file)
  eps <- 0.2
  time_size <- nrow(data)
  Simul <-
    data.frame(
      time = rep(data$time, n),
      nvers = rep(NA, n*time_size),
      Statut = rep(data$Statut[1], n*time_size),
      Souche = rep(data$Souche[1], n*time_size),
      puits = rep("duplic", each = n*time_size),
      rep.intra = rep(data$rep.intra[1], n*time_size),
      rep.inter = rep(data$rep.inter[1], n*time_size),
      num_duplic = rep(1:n, each = time_size)
    )
  # cat("Simul sans data \n")
  # print(str(Simul))
  Simul <- rbind(Simul, data)
  # print("Simul avec data \n")
  # print(str(Simul))
  Simul$num_duplic <- factor(Simul$num_duplic)
  for (i in 1:n) {
    ## Indicateur de quand il faut changer p
    change_p_time <- FALSE
    change_p_time_bis <- FALSE
    ## Premier jour
    # cat("Premier jour \n")
    n_vers <- data$nvers[1]
    if (n_vers >= 35) {
      # cat("n_vers >= 35 \n")
      new_nvers <- n_vers +
        round(rnorm(1, mean = 0,
                    sd = ecart_type$sd[ecart_type$nvers == 35]))
      if(new_nvers > 0){
        Simul$nvers[(i - 1)*time_size + 1] <- new_nvers
      }else{
        Simul$nvers[(i - 1)*time_size + 1] <- new_nvers+1
      }
    }
    if (n_vers < 35 & n_vers >= 10) {
      # cat("n_vers < 35 & n_vers >= 10 \n")
      new_nvers <- n_vers +
        round(rnorm(1, mean = 0,
                    sd = ecart_type$sd[ecart_type$nvers == data$nvers[1]]))
      if(new_nvers > 0){
        Simul$nvers[(i - 1)*time_size + 1] <- new_nvers
      }else{
        Simul$nvers[(i - 1)*time_size + 1] <- new_nvers+1
      }
    }
    if (n_vers < 10) {
      # cat("n_vers < 10 \n")
      if (n_vers == 0) {
        # Simul$nvers[(i - 1)*time_size + j] <- n_vers
        Simul$nvers[(i - 1)*time_size + 1] <- n_vers # PMG 2023-05-23 j -> 1 non ?
      } else {
        new_nvers <- n_vers +
          round(rtruncnorm(1, a = - (n_vers + eps), b = Inf, 
                           mean = 0,
                           sd = ecart_type$sd[ecart_type$nvers == data$nvers[1]]))
        Simul$nvers[(i - 1)*time_size + 1] <- new_nvers
      }
      Simul$nvers[(i - 1)*time_size + 1] <- new_nvers
    }
    ## Autres jours
    for (j in 2:time_size) {
      n_vers <- data$nvers[j]
      nvers_before <- Simul$nvers[(i - 1)*time_size + j - 1]
      ## ?cart ? ne pas d?passer pour garder la d?croissance du nb de vers
      b_trunc <- nvers_before - n_vers
      ## nb de deces entre jour j-1 et jour j
      n_deces <- data$nvers[j-1] - data$nvers[j]
      ## si choix == TRUE alors on prend la loi normale
      ## si choix == FALSE on prend une distribution de dirac en 0
      if (change_p_time == FALSE){
        choix <- rbernoulli(1, p = 0.05)
      }
      if (n_deces > 1) {
        change_p_time <- TRUE
      }
      if (n_deces >= 3) {
        change_p_time_bis <- TRUE
      }
      if (change_p_time == TRUE & change_p_time_bis == FALSE) {
        # cat("p = 0.5 \n")
        choix <- rbernoulli(1, p = 0.5)
      } 
      if (change_p_time_bis == TRUE) {
        # cat("p = 0.9 \n")
        choix <- rbernoulli(1, p = 0.975)
      }
      if (choix == TRUE) {
        if (n_vers >= 35) {
          # cat("n_vers >= 35 \n")
          if (b_trunc <=3) {
            # cat("b_trunc <=5 \n")
            new_nvers <- n_vers +
              round(rtruncnorm(1, a = -n_vers, b = (b_trunc+eps), mean = 0,
                               sd = ecart_type$sd[ecart_type$nvers == 35]
              ))
          } else {
            new_nvers <- n_vers +
              round(rtruncnorm(1, a = -eps, b = (b_trunc+eps), mean = 0,
                               sd = ecart_type$sd[ecart_type$nvers == 35]
              ))
          }
          Simul$nvers[(i - 1)*time_size + j] <- new_nvers
        }
        if (n_vers < 35 & n_vers >= 10) {
          # cat("n_vers < 35 & n_vers >= 10 \n")
          new_nvers <- n_vers +
            round(rtruncnorm(1, a = - (n_vers + eps), b = (b_trunc+eps), 
                             mean = 0,
                             sd = ecart_type$sd[ecart_type$nvers == data$nvers[j]]))
          Simul$nvers[(i - 1)*time_size + j] <- new_nvers
        }
        if (n_vers < 10) {
          # cat("n_vers < 10 \n")
          if (n_vers == 0) {
            Simul$nvers[(i - 1)*time_size + j] <- n_vers
          } else {
            new_nvers <- n_vers +
              round(rtruncnorm(1, a = - (n_vers + eps), b = (b_trunc+eps), 
                               mean = 0,
                               sd = ecart_type$sd[ecart_type$nvers == data$nvers[j]]))
            Simul$nvers[(i - 1)*time_size + j] <- new_nvers
          }
        }
      } else {
        # cat("pas de deces \n")
        Simul$nvers[(i - 1)*time_size + j] <- nvers_before
      }
    }
  }
  return(Simul)
}

#### Clustering ----
funFEM_ICL_old <- function(data){
  # cat("- Importation du package utils.\n")
  library(utils)
  #Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data)
  # Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  # Choix du modele
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  ICLs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  # cat("- Classifications avec plusieurs paramètres différents.\n")
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)
      Ks[i,j] <- model_choice_tmp$K
      ICLs[i,j] <- model_choice_tmp$icl
      # update progress bar
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  # Choix du modele
  # cat("- Sélection des paramètres de la meilleure classification.\n")
  choice_K <- apply(Ks,2,table)
  ICL <- apply(ICLs,2,mean)
  m <- models[which.max(ICL)]
  K <- as.numeric(names(which.max(choice_K[[which.max(ICL)]])))
  # Ajustement du modele
  # cat("- Ajustement du modèle avec la meilleure classification.\n")
  fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=TRUE)
  # Cls
  # cat("- Création de la variable avec les classes. \n")
  CLS <- list() ; length(CLS) <- nrow(data)
  for (i in 1:length(fitted_model$cls)) {
    for(j in 1:simulated_data.time_number) {
      CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
    }
  }
  data$Cls <- unlist(CLS)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Approx,fitted_model,file = Clustering_path)
  return(data)
}
funFEM_ICL_old <- function(data){
  # cat("- Importation du package utils.\n")
  library(utils)
  #Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data)
  # Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  # Choix du modele
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  ICLs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  # cat("- Classifications avec plusieurs paramètres différents.\n")
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)
      Ks[i,j] <- model_choice_tmp$K
      ICLs[i,j] <- model_choice_tmp$icl
      # update progress bar
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  # Choix du modele
  # cat("- Sélection des paramètres de la meilleure classification.\n")
  choice_K <- apply(Ks,2,table)
  ICL <- apply(ICLs,2,mean)
  m <- models[which.max(ICL)]
  K <- as.numeric(names(which.max(choice_K[[which.max(ICL)]])))
  # Ajustement du modele
  # cat("- Ajustement du modèle avec la meilleure classification.\n")
  fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=TRUE)
  # Cls
  # cat("- Création de la variable avec les classes. \n")
  CLS <- list() ; length(CLS) <- nrow(data)
  for (i in 1:length(fitted_model$cls)) {
    for(j in 1:simulated_data.time_number) {
      CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
    }
  }
  data$Cls <- unlist(CLS)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Approx,fitted_model,file = Clustering_path)
  return(data)
}

## Fpca
Fpca_kmeans <- function(data) {
  # Nouvelle structure des donnees
  # cat("- Restructuration des données. \n")
  data2 <- list() ; length(data2) <- length(unique(data$Wells))
  x <- list() ; length(x) <- length(data2)
  for (i in 1:length(data2)){
    data2[[i]] <- data$surv[which(data$Wells == unique(data$Wells)[i])]
    x[[i]] <- Time
  }
  # ACP fonctionnelle
  # cat("- Réalisation de l'ACP fonctionnelle. \n")
  suppressWarnings(res <- fdapace::FPCA(data2,x))
  # Clustering avec kmeans 
  # cat("- Centrer et réduire les données. \n")
  PC <- scale(res$xiEst)
  # cat("- Calcul du nombre de clusters. \n")
  Nclust <- min(max_K,nrow(PC))
  test <- try(K <- NbClust::NbClust(PC,min.nc = min_K, max.nc = Nclust, method = "kmeans",index = "silhouette")$Best.nc[[1]], silent = TRUE)
  while(class(test) != "numeric" & Nclust > 2) {
    Nclust <- Nclust - 1
    test <- try(K <- NbClust::NbClust(PC, max.nc = Nclust, method = "kmeans",index = "silhouette")$Best.nc[[1]],silent=TRUE)
  }
  if (class(test) != "numeric") {K <- 2}
  # cat("- Application des kmeans. \n")
  Km <- kmeans(PC,centers = K)
  # Groupes
  # cat("- Création de la variable avec les classes. \n")
  Cls <- list() ; length(Cls) <- nrow(data)
  Cluster <- Km$cluster
  for (i in 1:length(Cluster)) {
    for (j in 1:simulated_data.time_number) {
      Cls[(i-1)*simulated_data.time_number+j] <- Cluster[i]
    }
  }
  data$Cls <- unlist(Cls)
  data$Cls <- factor(data$Cls)
  ## Sauvegarde
  # cat("- Sauvegarde des résultats. \n")
  save(K, Km, PC, res, file = Clustering_path)
  return(data)
}

Fpca_GMM <- function(data) {
  # cat("- Importation du package mclust. \n")
  suppressMessages(library(mclust,quietly=TRUE)) 
  # Nouvelle structure des donnees
  # cat("- Restructuration des données. \n")
  data2 <- list() ; length(data2) <- length(unique(data$Wells))
  x <- list() ; length(x) <- length(data2)
  for (i in 1:length(data2)){
    data2[[i]] <- data$surv[which(data$Wells == unique(data$Wells)[i])]
    x[[i]] <- Time
  }
  # ACP fonctionnelle
  # cat("- Réalisation de l'ACP fonctionnelle. \n")
  suppressWarnings(res <- fdapace::FPCA(data2,x))
  # Clustering avec GMM 
  # cat("- Centrer et réduire les données. \n")
  PC <- scale(res$xiEst)
  # cat("- Application de Gaussian Mixture Model. \n")
  Gmm <- Mclust(PC, G = min_K:min(max_K,nrow(PC)))
  K <- Gmm$G
  # Groupes
  # cat("- Création de la variable avec les classes. \n")
  Cls <- list() ; length(Cls) <- nrow(data)
  Cluster <- Gmm$classification
  for (i in 1:length(Cluster)) {
    for (j in 1:simulated_data.time_number) {
      Cls[(i-1)*simulated_data.time_number+j] <- Cluster[i]
    }
  }
  data$Cls <- unlist(Cls)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Gmm,PC,res,file = Clustering_path)
  return(data)
}

## Bsplines
Bsplines_kmeans <- function(data){
  #Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data)
  #Base de Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Centrer et réduire les données. \n")
  sCoefs <- scale(fd$coefs)
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  #Kmeans
  # cat("- Calcul du nombre de clusters. \n")
  Nclust <- min(max_K,nrow(sCoefs))
  test <- try(K <- NbClust::NbClust(sCoefs, min.nc = min_K, max.nc = Nclust, method = "kmeans", 
                                    index = "silhouette")$Best.nc[[1]], silent = TRUE)
  while(class(test) != "numeric" & Nclust > 2) {
    Nclust <- Nclust - 1
    test <- try(K <- NbClust::NbClust(sCoefs, max.nc = Nclust, method = "kmeans",index = "silhouette")$Best.nc[[1]],silent=TRUE)
  }
  if (class(test) != "numeric") {K <- 2}
  # cat("- Application des kmeans. \n")
  Km <- kmeans(t(sCoefs),centers = K)
  # Groupes
  # cat("- Création de la variable avec les classes. \n")
  Cls <- list() ; length(Cls) <- nrow(data)
  Cluster <- Km$cluster
  for (i in 1:length(Cluster)) {
    for (j in 1:simulated_data.time_number) {
      Cls[(i-1)*simulated_data.time_number+j] <- Cluster[i]
    }
  }
  data$Cls <- unlist(Cls)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Km,sCoefs,Approx, file = Clustering_path)
  return(data)
}

Bsplines_GMM <- function(data){
  suppressMessages(library(mclust,quietly = TRUE))
  
  
  data$Cls <- rep(1,nrow(data)) # PMG 2023-05-19 : faire un group control separé
  data$Cls[data$Treat != "Temoin"] <- NA # PMG 2023-05-19 : faire un group control separé
  data_without_control <- data[data$Treat != "Temoin",] # PMG 2023-05-19 : faire un group control separé
  n_Temoin <- length(unique(data$Wells[data$Treat =="Temoin"]))
  
  # Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data_without_control)
  # Base de Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Centrer et réduire les données. \n")
  sCoefs <- scale(fd$coefs)
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  # GMM
  # cat("- Application de Gaussian Mixture Model. \n")
  Gmm <- Mclust(t(sCoefs), G = max((min_K-1),2):(max_K-1))
  K <- Gmm$G
  # Groupes
  # cat("- Création de la variable avec les classes. \n")
  Cls <- list() ; length(Cls) <- nrow(data2)
  Cluster <- Gmm$classification
  for (i in 1:length(Cluster)) {
    for (j in 1:simulated_data.time_number) {
      Cls[(i-1)*simulated_data.time_number+j] <- Cluster[i]
    }
  }
  data$Cls[data$Treat != "Temoin"] <- unlist(Cls)+1 # PMG 2023-05-19 : faire un group control separé
  # data$Cls <- unlist(Cls) # PMG 2023-05-19 : faire un group control separé
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  
  #### PMG 2023-05-19 : faire un group control separé
  Gmm_tmp <- matrix(0,ncol=K+1,nrow=n_Temoin+nrow(  Gmm$z ))
  data_tmp_test <- data[,c("Treat","Wells")]
  test <- aggregate(data_tmp_test,by = list(data_tmp_test$Wells),FUN=function(v) v[1])
  test <- test$Treat
  index_not_control <- which(test != "Temoin")
  Gmm_tmp[ index_not_control , ] <- cbind(rep(0,nrow(Gmm$z)) ,Gmm$z)
  Gmm_tmp[ -index_not_control,1] <- 1
  Gmm$z <- Gmm_tmp
  #### PMG 2023-05-19 : faire un group control separé
  
  K <- K+1                                                                # PMG 2023-05-19 : faire un group control separé
  save(K,Gmm,sCoefs,Approx, file = Clustering_path)
  return(data)
}

## funFEM
funFEM_AIC_modified <- function(data){
  # cat("- Importation du package utils.\n")
  library(utils)
  
  
  data$Cls <- rep(1,nrow(data)) # PMG 2023-05-19 : faire un group control separé
  data$Cls[data$Treat != "Temoin"] <- NA # PMG 2023-05-19 : faire un group control separé
  data_without_control <- data[data$Treat != "Temoin",] # PMG 2023-05-19 : faire un group control separé
  n_Temoin <- length(unique(data$Wells[data$Treat =="Temoin"]))
  
  data2 <- transform(data_without_control)
  
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)  
  # Choix du modele
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  AICs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  # cat("- Classifications avec plusieurs paramètres différents.\n")
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)
      Ks[i,j] <- model_choice_tmp$K
      AICs[i,j] <- model_choice_tmp$aic
      # update progress bar
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  
  choice_K <- apply(Ks,2,table)
  AIC <- apply(AICs,2,mean)
  m <- models[which.min(AIC)]
  K <- as.numeric(names(which.max(choice_K[[which.min(AIC)]])))
  # Ajustement du modele
  # cat("- Ajustement du modèle avec la meilleure classification.\n")
  options(show.error.messages = FALSE)
  fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=FALSE)
  options(show.error.messages = TRUE)
  # Cls
  # cat("- Création de la variable avec les classes. \n")
  CLS <- list() ; length(CLS) <- nrow(data)
  for (i in 1:length(fitted_model$cls)) {
    for(j in 1:simulated_data.time_number) {
      CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
    }
  }
  # data$Cls <- unlist(CLS)
  data$Cls[data$Treat != "Temoin"] <- unlist(CLS)+1 # PMG 2023-05-19 : faire un group control separé
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  
  
  #### PMG 2023-05-19 : faire un group control separé
  P_tmp <- matrix(0,ncol=K+1,nrow=n_Temoin+nrow(  fitted_model$P ))
  data_tmp_test <- data[,c("Treat","Wells")]
  test <- aggregate(data_tmp_test,by = list(data_tmp_test$Wells),FUN=function(v) v[1])
  test <- test$Treat
  index_not_control <- which(test != "Temoin")
  P_tmp[ index_not_control , ] <- cbind(rep(0,nrow(fitted_model$P)) ,fitted_model$P)
  P_tmp[ -index_not_control,1] <- 1
  fitted_model$P <- P_tmp
  #### PMG 2023-05-19 : faire un group control separé
  
  
  K <- K+1                            
  save(K,Approx,fitted_model,file = Clustering_path)  
  
  
  # 
  # 
  # 
  # 
  # 
  # 
  # #Creation de la data.frame adaptee
  # # cat("- Restructuration des données. \n")
  # data2 <- transform(data)
  # # Bsplines
  # # cat("- Création d'une base de Bsplines. \n")
  # Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # # cat("- Évaluation des coefficients des courbes dans la base. \n")
  # fd <- smooth.basis(Time,data2,Base)$fd
  # # cat("- Estimation des courbes. \n")
  # Approx <- eval.fd(Time,fd)  
  # # Choix du modele
  # models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  # Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  # AICs <- matrix(NA,ncol=length(models),nrow=n_rep)
  # options(show.error.messages = FALSE)
  # # cat("- Classifications avec plusieurs paramètres différents.\n")
  # pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  # for(j in 1:length(models)){
  #   for(i in 1:n_rep){
  #     model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)
  #     Ks[i,j] <- model_choice_tmp$K
  #     AICs[i,j] <- model_choice_tmp$aic
  #     # update progress bar
  #     setTxtProgressBar(pb, (j-1)*n_rep+i)
  #   }
  # }
  # close(pb)
  # options(show.error.messages = TRUE)
  # # cat("- Sélection des paramètres de la meilleure classification.\n")
  # choice_K <- apply(Ks,2,table)
  # AIC <- apply(AICs,2,mean)
  # m <- models[which.min(AIC)]
  # K <- as.numeric(names(which.max(choice_K[[which.min(AIC)]])))
  # # Ajustement du modele
  # # cat("- Ajustement du modèle avec la meilleure classification.\n")
  # options(show.error.messages = FALSE)
  # fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=FALSE)
  # options(show.error.messages = TRUE)
  # # Cls
  # # cat("- Création de la variable avec les classes. \n")
  # CLS <- list() ; length(CLS) <- nrow(data)
  # for (i in 1:length(fitted_model$cls)) {
  #   for(j in 1:simulated_data.time_number) {
  #     CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
  #   }
  # }
  # data$Cls <- unlist(CLS)
  # data$Cls <- factor(data$Cls)
  # # cat("- Sauvegarde des résultats. \n")
  # save(K,Approx,fitted_model,file = Clustering_path)
  return(data)
}

funFEM_AIC <- function(data){
  # cat("- Importation du package utils.\n")
  library(utils)
  #Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data)
  # Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  # Choix du modele
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  AICs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  # cat("- Classifications avec plusieurs paramètres différents.\n")
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      try({model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)},silent=TRUE)
      Ks[i,j] <- model_choice_tmp$K
      AICs[i,j] <- model_choice_tmp$aic
      # update progress bar
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  ## Choix des paramètres
  # cat("- Sélection des paramètres de la meilleure classification.\n")
  choice_K <- apply(Ks,2,table)
  AIC <- apply(AICs,2,mean)
  m <- models[which.min(AIC)]
  K <- as.numeric(names(which.max(choice_K[[which.min(AIC)]])))
  # Ajustement du modele
  # cat("- Ajustement du modèle avec la meilleure classification.\n")
  fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=TRUE)
  # Cls
  # cat("- Création de la variable avec les classes. \n")
  CLS <- list() ; length(CLS) <- nrow(data)
  for (i in 1:length(fitted_model$cls)) {
    for(j in 1:simulated_data.time_number) {
      CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
    }
  }
  data$Cls <- unlist(CLS)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Approx,fitted_model,file = Clustering_path)
  return(data)
}

funFEM_BIC <- function(data){
  # cat("- Importation du package utils.\n")
  library(utils)
  #Creation de la data.frame adaptee
  # cat("- Restructuration des données. \n")
  data2 <- transform(data)
  # Bsplines
  # cat("- Création d'une base de Bsplines. \n")
  Base <- create.bspline.basis(rangeval = c(0,time_max), nbasis = n_functional_base)
  # cat("- Évaluation des coefficients des courbes dans la base. \n")
  fd <- smooth.basis(Time,data2,Base)$fd
  # cat("- Estimation des courbes. \n")
  Approx <- eval.fd(Time,fd)
  # Choix du modele
  models <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB")
  Ks <- matrix(NA,ncol=length(models),nrow=n_rep)
  BICs <- matrix(NA,ncol=length(models),nrow=n_rep)
  options(show.error.messages = FALSE)
  # cat("- Classifications avec plusieurs paramètres différents.\n")
  pb <- txtProgressBar(min = 0, max = length(models)*n_rep, style = 3)
  for(j in 1:length(models)){
    for(i in 1:n_rep){
      try({model_choice_tmp <- funFEM(fd,K=min_K:max_K,model=models[j],init='kmeans',lambda=0,disp=FALSE)},silent=TRUE)
      Ks[i,j] <- model_choice_tmp$K
      BICs[i,j] <- model_choice_tmp$bic
      # update progress bar
      setTxtProgressBar(pb, (j-1)*n_rep+i)
    }
  }
  close(pb)
  options(show.error.messages = TRUE)
  ## Choix des paramètres
  # cat("- Sélection des paramètres de la meilleure classification.\n")
  choice_K <- apply(Ks,2,table)
  BIC <- apply(BICs,2,mean)
  m <- models[which.min(BIC)]
  K <- as.numeric(names(which.max(choice_K[[which.min(BIC)]])))
  # Ajustement du modele
  # cat("- Ajustement du modèle avec la meilleure classification.\n")
  fitted_model <- funFEM(fd,K=K,model=m,init='kmeans',lambda=0,disp=TRUE)
  # Cls
  # cat("- Création de la variable avec les classes. \n")
  CLS <- list() ; length(CLS) <- nrow(data)
  for (i in 1:length(fitted_model$cls)) {
    for(j in 1:simulated_data.time_number) {
      CLS[(i-1)*simulated_data.time_number+j] <- fitted_model$cls[i]
    }
  }
  data$Cls <- unlist(CLS)
  data$Cls <- factor(data$Cls)
  # cat("- Sauvegarde des résultats. \n")
  save(K,Approx,fitted_model,file = Clustering_path)
  return(data)
}

#### Compute proba ----

## Prends un character et une liste de chacacter
## Renvoie TRUE si un élément de la liste est égal à name et FALSE sinon
Same_name <- function(name, names_list){
  if (name %in% names_list){
    return (TRUE)
  }
  for (i in 1:length(names_list)) {
    my_list <- names_list
    if (names_list[i] %in% my_list[-i]){
      return (TRUE)
    }
  }
  return (FALSE)
}

Pb_fpca_kmeans <- function(Data){
  ## Chargement des objets necessaires
  # cat("- Importation des données nécessaires. \n")
  load(Clustering_path)
  threshold_deviance <- 0.01
  ## Distance usuelle entre une courbe et un cluster
  # cat("- Calcul des distances entre une courbe et un cluster. \n")
  n_Wells <- length(unique(Data$Wells))
  dist <- matrix(NA,nrow = n_Wells, ncol = K)
  dist2 <- matrix(NA,nrow = nrow(dist), ncol = K)
  for (i in 1:n_Wells) {
    for (j in 1:K) {
      dist[i,j] <- sum((Km$centers[j, ] - PC[i,])^2)
      dist2[i,j] <- dist[i,j]^2
    }
  }
  Dist2 <- data.frame(dist2, cls = Km$cluster, row.names = unique(data_transformed$Wells))
  ## Calcul du nombre de courbes dans chaque cluster
  # cat("- Vérification qu'il n'y a pas de cluster avec une seule courbe. \n")
  TAB <- table(Dist2$cls)
  Tau <- TRUE
  i <- 1
  while ((Tau == TRUE) & (i <= length(TAB))){
    if (TAB[[i]]<2) {Tau <- FALSE} 
    i <- i+1
  }
  if (Tau == TRUE) {
    ## Calcul de tau
    # cat("- Calcul de tau. \n")
    tau <- list() ; length(tau) <- K
    for (i in 1:K){
      data <- Dist2[which(Dist2$cls == i),]
      tau[i] <-1/(dim(data)[1]-1)*sum(data[,i])
    }
    ## Calcul proba
    # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
    M <- matrix(NA,nrow=nrow(dist), ncol=K)
    Denom <- list() ; length(Denom) <- nrow(dist)
    Pbcurves <- matrix(NA,nrow=nrow(dist), ncol=K)
    for (i in 1:nrow(dist)) {
      for (j in 1:K) {
        M[i,j] <- (1/sqrt(tau[[Dist2$cls[i]]]))*exp((-1/(2*tau[[Dist2$cls[i]]]))*Dist2[i,j])
      }
      Denom[i] <- sum(M[i,])
      for (j in 1:K) {
        Pbcurves[i,j] <- M[i,j]/Denom[[i]]
      }
    }
  } else {
    M <- matrix(NA,nrow=nrow(dist), ncol=K)
    Denom <- list() ; length(Denom) <- nrow(dist)
    Pbcurves <- matrix(NA,nrow=nrow(dist), ncol=K)
    # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
    for (i in 1:nrow(dist)) {
      for (j in 1:K) {
        M[i,j] <- exp((-1/2)*Dist2[i,j])
      }
      Denom[i] <- sum(M[i,])
      for (j in 1:K) {
        Pbcurves[i,j] <- M[i,j]/Denom[[i]]
      }
    }
  }
  ## Creation cond
  my_Wells <- unique(Data$Wells)
  Cond <- list() ; length(Cond) <- length(my_Wells)
  for (i in 1:length(my_Wells)) {
    condition <- Data$cond[which(Data$Wells == my_Wells[i])][1]
    test <- strsplit(as.character(condition),split = "\n")[1]
    Cond[i] <- test[[1]]
  }
  ## Creation dataframe
  Proba_curves <- data.frame(Pb = Pbcurves, cond = unlist(Cond))
  Proba_curves_bis <- data.frame(Pbcurves)
  # print("Proba_curves")
  # print(head(Proba_curves))
  ## Proba Treatecules
  n_Treat <- length(treat_group)
  Proba_Treat <- data.frame(Pb = rep(NA,n_Treat*K), Cls = rep(NA,n_Treat*K), Treat = rep(NA,n_Treat*K))
  #Pour la 1ère Treatécule
  # cat("- Calcul des probabilités que la 1ère Treatécule soit dans chaque cluster. \n")
  Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[1]],K))
  for (j in 1:K){
    dat <- Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[1]),]
    if (length(treat_group[[1]])>1) {
      for (k in 2:length(unlist(treat_group[[1]]))) {
        dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[k]),])
      }
    }
    Proba$Pb[j] <- 1/(dim(dat)[1])*sum(dat[,j])
    Proba$Cls[j] <- j
    Proba$Treat[j] <- treat_names[[1]]
  }
  Proba_Treat <- Proba
  ## Les autres Treatécules
  # cat("- Calcul des probabilités que chaque Treatécule soit dans chaque cluster. \n")
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
      Proba$Cls[j] <- j
      Proba$Treat[j] <- treat_names[[i]]
    }
    Proba_Treat <- rbind(Proba_Treat,Proba)
  }
  # cat("Proba_Treat \n")
  # print(Proba_Treat)
  ## Noms clusters
  name_cls <- list() ; length(name_cls) <- K
  sign <- matrix(NA, nrow=K, ncol = n_cut)
  for (i in 1:K) {
    if (i==1) {
      # cat("- Caractérisation du 1er cluster. \n",sep="")
    } else {
      # cat(cat("- Caractérisation du",i),"ème cluster. \n",sep="")
    }
    d <- Data[which(Data$Cls==i),]
    for (j in 1:n_cut) {
      t_m <- ((j-1)*time_max)/n_cut
      t_p <- ((j*time_max)/n_cut)
      Mean <- median(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
      if (Mean <threshold_deviance & Mean>=-threshold_deviance) sign[i,j] <- "="
      if (Mean >= threshold_deviance) sign[i,j] <- "+"
      if (Mean < -threshold_deviance) sign[i,j] <- "-"
    }
    if (Proba_Treat$Pb[which(Proba_Treat$Cls == i & Proba_Treat$Treat == "Temoin")]>0.75) {
      name_cls[i] <- "Témoin"
    }
    else {
      name <- sign[i,1]
      for (k in 2:n_cut){
        name <- paste(name,sign[i,k],sep=",")
      }
      name_cls[i] <- name
    }
    if (i>1){
      # cat("- Vérification qu'il n'y ait pas deux fois le même nom. \n")
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
      # cat("Plusieurs clusters avec le même nom ? \n")
      # print(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)])))
      while(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)]))==TRUE){
        for (j in 2:i){
          for (l in 1:(j-1)){
            if (name_cls[[j]]==name_cls[[l]]){
              ## Valeur absolue de l'ecart entre les deux médianes
              Tab <- list(); length(Tab) <- n_cut
              for (m in 1:n_cut){
                dj <- Data[which(Data$Cls==j),]
                dl <- Data[which(Data$Cls==l),]
                t_m <- ((m-1)*time_max)/n_cut
                t_p <- (m*time_max)/n_cut
                Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
                Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
                Tab[m] <- abs(Meanj-Meanl)
              }
              ## Localisation du plus grand ecart entre les deux groupes
              k <- which.max(unlist(Tab))
              t_m <- ((k-1)*time_max)/n_cut
              t_p <- (k*time_max)/n_cut
              Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
              Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
              if (Meanj>=threshold_deviance) {
                if (Meanl>=Meanj){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if (Meanj<= -threshold_deviance) {
                if (Meanj<=Meanl){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } 
                else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in 2:n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if ((Meanj < threshold_deviance) & (Meanj > -threshold_deviance)){
                if (Meanj >= 0 & Meanl >= 0) {
                  if (Meanj >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl <= 0) {
                  if (Meanj <= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj >= 0 & Meanl <= 0) {
                  if (Meanj >= abs(Meanl)) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl >= 0) {
                  if (abs(Meanj) >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
              }
            }
          }
        }
      }
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
    }
  }
  Proba_Treat$Cls <- unlist(name_cls)
  # cat("- Modification des noms de chaque cluster dans le dataframe \n")
  new_names <- list() ; length(new_names) <- nrow(data_cluster)
  for (i in 1:K) {
    new_names[which(data_cluster$Cls == i)] <- name_cls[[i]]
  }
  data4 <- data_cluster
  data4$Cls <- unlist(new_names)
  ## Dataframe avec l'écart moyen de survie par classe
  # cat("- Calcul de l'écart moyen de survie par classe \n")
  # cat("K = \n")
  # print(K)
  # cat("unique(data4$Cls) _n")
  # print(unique(data4$Cls))
  mean_ecart <- data.frame(time = rep(Time,K),
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
  mean_ecart$surv <- unlist(my_mean)
  ## Sauvegarde
  # cat("- Sauvegarde des résultats. \n")
  path <- paste0(results_path,"/Proba_Treat.RData")
  save(Proba_Treat, file = path)
  save(Proba_curves, Proba_curves_bis, Proba_Treat, mean_ecart, file = Proba_path)
  return(data4)
}

Pb_fpca_GMM <- function(Data){
  ## Chargement des objets necessaires
  # cat("- Importation des données nécessaires. \n")
  load(Clustering_path)
  threshold_deviance <- 0.01
  # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
  ## Creation cond
  my_Wells <- unique(Data$Wells)
  Cond <- list() ; length(Cond) <- length(my_Wells)
  for (i in 1:length(my_Wells)) {
    condition <- Data$cond[which(Data$Wells == my_Wells[i])][1]
    test <- strsplit(as.character(condition),split = "\n")[1]
    Cond[i] <- test[[1]]
  }
  Proba_curves <- data.frame(Pb=Gmm$z, cond = unlist(Cond))
  Proba_curves_bis <- data.frame(Gmm$z)
  # cat("Proba_curves \n")
  # print(Proba_curves)
  ## Proba Treatecules
  n_Treat <- length(treat_group)
  Proba_Treat <- data.frame(Pb = rep(NA,n_Treat*K), Cls = rep(NA,n_Treat*K), Treat = rep(NA,n_Treat*K))
  #Pour la 1ère Treatécule
  # cat("- Calcul des probabilités que la 1ère Treatécule soit dans chaque cluster. \n")
  Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[1]],K))
  for (j in 1:K){
    dat <- Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[1]),]
    if (length(treat_group[[1]])>1) {
      for (k in 2:length(unlist(treat_group[[1]]))) {
        dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[k]),])
      }
    }
    Proba$Pb[j] <- 1/(dim(dat)[1])*sum(dat[,j])
    Proba$Cls[j] <- j
    Proba$Treat[j] <- treat_names[[1]]
  }
  Proba_Treat <- Proba
  # cat("Proba \n")
  # print(Proba_Treat)
  ## Les autres Treatécules
  # cat("- Calcul des probabilités que chaque Treatécule soit dans chaque cluster. \n")
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
  # cat("Proba_Treat \n")
  # print(Proba_Treat)
  ## Noms clusters
  name_cls <- list() ; length(name_cls) <- K
  sign <- matrix(NA, nrow=K, ncol = n_cut)
  for (i in 1:K) {
    # cat(cat("- Caractérisation du",i),"ème cluster. \n",sep="")
    d <- Data[which(Data$Cls==i),]
    for (j in 1:n_cut) {
      t_m <- ((j-1)*time_max)/n_cut
      t_p <- ((j*time_max)/n_cut)
      Mean <- median(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
      if (Mean <threshold_deviance & Mean>=-threshold_deviance) sign[i,j] <- "="
      if (Mean >= threshold_deviance) sign[i,j] <- "+"
      if (Mean < -threshold_deviance) sign[i,j] <- "-"
    }
    if (Proba_Treat$Pb[which(Proba_Treat$Cls == i & Proba_Treat$Treat == "Temoin")]>0.75) {
      name_cls[i] <- "Témoin"
    }
    else {
      name <- sign[i,1]
      for (k in 2:n_cut){
        name <- paste(name,sign[i,k],sep=",")
      }
      name_cls[i] <- name
    }
    if (i>1){
      # cat("- Vérification qu'il n'y ait pas deux fois le même nom. \n")
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
      # cat("Plusieurs clusters avec le même nom ? \n")
      # print(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)])))
      while(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)]))==TRUE){
        for (j in 2:i){
          for (l in 1:(j-1)){
            if (name_cls[[j]]==name_cls[[l]]){
              ## Valeur absolue de l'ecart entre les deux médianes
              Tab <- list(); length(Tab) <- n_cut
              for (m in 1:n_cut){
                dj <- Data[which(Data$Cls==j),]
                dl <- Data[which(Data$Cls==l),]
                t_m <- ((m-1)*time_max)/n_cut
                t_p <- (m*time_max)/n_cut
                Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
                Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
                Tab[m] <- abs(Meanj-Meanl)
              }
              ## Localisation du plus grand ecart entre les deux groupes
              k <- which.max(unlist(Tab))
              t_m <- ((k-1)*time_max)/n_cut
              t_p <- (k*time_max)/n_cut
              Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
              Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
              if (Meanj>=threshold_deviance) {
                if (Meanl>=Meanj){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if (Meanj<= -threshold_deviance) {
                if (Meanj<=Meanl){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } 
                else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in 2:n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if ((Meanj < threshold_deviance) & (Meanj > -threshold_deviance)){
                if (Meanj >= 0 & Meanl >= 0) {
                  if (Meanj >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl <= 0) {
                  if (Meanj <= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj >= 0 & Meanl <= 0) {
                  if (Meanj >= abs(Meanl)) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl >= 0) {
                  if (abs(Meanj) >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
              }
            }
          }
        }
      }
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
    }
  }
  Proba_Treat$Cls <- unlist(name_cls)
  # cat("- Modification des noms de chaque cluster dans le dataframe \n")
  new_names <- list() ; length(new_names) <- nrow(data_cluster)
  for (i in 1:K) {
    new_names[which(data_cluster$Cls == i)] <- name_cls[[i]]
  }
  data4 <- data_cluster
  data4$Cls <- unlist(new_names)
  ## Dataframe avec l'écart moyen de survie par classe
  # cat("- Calcul de l'écart moyen de survie par classe \n")
  mean_ecart <- data.frame(time = rep(Time,K),
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
  mean_ecart$surv <- unlist(my_mean)
  ## Sauvegarde
  # cat("- Sauvegarde des résultats. \n")
  path <- paste0(results_path,"/Proba_Treat.RData")
  save(Proba_Treat, file = path)
  save(Proba_curves, Proba_curves_bis, Proba_Treat, mean_ecart, file = Proba_path)
  return(data4)
}

Pb_Bsplines_kmeans <- function(Data){
  ## Chargement des objets necessaires
  # cat("- Importation des données nécessaires. \n")
  load(Clustering_path)
  threshold_deviance <-  0.01
  ## Creation cond
  my_Wells <- unique(Data$Wells)
  Cond <- list() ; length(Cond) <- length(my_Wells)
  for (i in 1:length(my_Wells)) {
    condition <- Data$cond[which(Data$Wells == my_Wells[i])][1]
    test <- strsplit(as.character(condition),split = "\n")[1]
    Cond[i] <- test[[1]]
  }
  ## Distance usuelle entre une courbe et un cluster
  # cat("- Calcul des distances entre une courbe et un cluster. \n")
  dist <- matrix(NA,nrow = ncol(sCoefs), ncol = K)
  dist2 <- matrix(NA,nrow = ncol(sCoefs), ncol = K)
  for (i in 1:ncol(sCoefs)) {
    for (j in 1:K) {
      dist[i,j] <- sum((Km$centers[j, ] - sCoefs[,i])^2)
      dist2[i,j] <- dist[i,j]^2
    }
  }
  Dist2 <- data.frame(dist2, cls = Km$cluster)
  ## Calcul du nombre de courbes dans chaque cluster
  # cat("- Vérification qu'il n'y a pas de cluster avec une seule courbe. \n")
  TAB <- table(Dist2$cls)
  Tau <- TRUE
  i = 1
  while ((Tau == TRUE) & (i <= length(TAB))){
    if (TAB[[i]]<2) {Tau <- FALSE} 
    i <- i+1
  }
  if (Tau == TRUE) {
    ## Calcul de tau
    # cat("- Calcul de tau. \n")
    tau <- list() ; length(tau) <- K
    for (i in 1:K){
      d <- Dist2[which(Dist2$cls == i),]
      tau[i] <-1/(dim(d)[1]-1)*sum(d[,i])
    }
    ## Calcul proba
    M <- matrix(NA,nrow=ncol(sCoefs), ncol=K)
    Denom <- list() ; length(Denom) <- ncol(sCoefs)
    Pbcurves <- matrix(NA,nrow=ncol(sCoefs), ncol=K)
    # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
    for (i in 1:ncol(sCoefs)) {
      for (j in 1:K) {
        M[i,j] <- (1/sqrt(tau[[Dist2$cls[i]]]))*exp((-1/(2*tau[[Dist2$cls[i]]]))*Dist2[i,j])
      }
      Denom[i] <- sum(M[i,])
      for (j in 1:K) {
        Pbcurves[i,j] <- M[i,j]/Denom[[i]]
      }
    }
  } else {
    ## Calcul proba
    M <- matrix(NA,nrow=ncol(sCoefs), ncol=K)
    Denom <- list() ; length(Denom) <- ncol(sCoefs)
    Pbcurves <- matrix(NA,nrow=ncol(sCoefs), ncol=K)
    # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
    for (i in 1:ncol(sCoefs)) {
      for (j in 1:K) {
        M[i,j] <- exp((-1/2)*Dist2[i,j])
      }
      Denom[i] <- sum(M[i,])
      for (j in 1:K) {
        Pbcurves[i,j] <- M[i,j]/Denom[[i]]
      }
    }
  }
  ## Creation dataframe
  Proba_curves <- data.frame(Pb = Pbcurves, cond = unlist(Cond))
  Proba_curves_bis <- data.frame(Pbcurves)
  ## Proba Treatecules
  n_Treat <- length(treat_group)
  Proba_Treat <- data.frame(Pb = rep(NA,n_Treat*K), Cls = rep(NA,n_Treat*K), Treat = rep(NA,n_Treat*K))
  #Pour la 1ère Treatécule
  # cat("- Calcul des probabilités que la 1ère Treatécule soit dans chaque cluster. \n")
  Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[1]],K))
  for (j in 1:K){
    dat <- Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[1]),]
    if (length(treat_group[[1]])>1) {
      for (k in 2:length(unlist(treat_group[[1]]))) {
        dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[k]),])
      }
    }
    Proba$Pb[j] <- 1/(dim(dat)[1])*sum(dat[,j])
    Proba$Cls[j] <- j
    Proba$Treat[j] <- treat_names[[1]]
  }
  Proba_Treat <- Proba
  ## Les autres Treatécules
  # cat("- Calcul des probabilités que chaque Treatécule soit dans chaque cluster. \n")
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
  ## Noms clusters
  name_cls <- list() ; length(name_cls) <- K
  sign <- matrix(NA, nrow=K, ncol = n_cut)
  for (i in 1:K) {
    # cat(cat("- Caractérisation du",i),"ème cluster. \n",sep="")
    d <- Data[which(Data$Cls==i),]
    for (j in 1:n_cut) {
      t_m <- ((j-1)*time_max)/n_cut
      t_p <- ((j*time_max)/n_cut)
      Mean <- median(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
      if (Mean <threshold_deviance & Mean>=-threshold_deviance) sign[i,j] <- "="
      if (Mean >= threshold_deviance) sign[i,j] <- "+"
      if (Mean < -threshold_deviance) sign[i,j] <- "-"
    }
    if (Proba_Treat$Pb[which(Proba_Treat$Cls == i & Proba_Treat$Treat == "Temoin")]>0.75) {
      name_cls[i] <- "Témoin"
    }
    else {
      name <- sign[i,1]
      for (k in 2:n_cut){
        name <- paste(name,sign[i,k],sep=",")
      }
      name_cls[i] <- name
    }
    if (i>1){
      # cat("- Vérification qu'il n'y ait pas deux fois le même nom. \n")
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
      # cat("Plusieurs clusters avec le même nom ? \n")
      # print(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)])))
      while(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)]))==TRUE){
        for (j in 2:i){
          for (l in 1:(j-1)){
            if (name_cls[[j]]==name_cls[[l]]){
              ## Valeur absolue de l'ecart entre les deux médianes
              Tab <- list(); length(Tab) <- n_cut
              for (m in 1:n_cut){
                dj <- Data[which(Data$Cls==j),]
                dl <- Data[which(Data$Cls==l),]
                t_m <- ((m-1)*time_max)/n_cut
                t_p <- (m*time_max)/n_cut
                Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
                Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
                Tab[m] <- abs(Meanj-Meanl)
              }
              ## Localisation du plus grand ecart entre les deux groupes
              k <- which.max(unlist(Tab))
              t_m <- ((k-1)*time_max)/n_cut
              t_p <- (k*time_max)/n_cut
              Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
              Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
              if (Meanj>=threshold_deviance) {
                if (Meanl>=Meanj){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"+")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if (Meanj<= -threshold_deviance) {
                if (Meanj<=Meanl){
                  name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[j,1])
                  for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[j]] <- name
                } 
                else {
                  name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                  if (k==1) {name <- paste(name,sign[l,1])
                  for (iter in 2:n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                  else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    name <- paste0(name,"-")
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                  }
                  name_cls[[l]] <- name
                }
              }
              if ((Meanj < threshold_deviance) & (Meanj > -threshold_deviance)){
                if (Meanj >= 0 & Meanl >= 0) {
                  if (Meanj >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl <= 0) {
                  if (Meanj <= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj >= 0 & Meanl <= 0) {
                  if (Meanj >= abs(Meanl)) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
                if (Meanj <= 0 & Meanl >= 0) {
                  if (abs(Meanj) >= Meanl) {
                    name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[j,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"-")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[j]] <- name
                  } else {
                    name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
                    if (k==1) {name <- paste(name,sign[l,1])
                    for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
                    else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                      name <- paste0(name,"+")
                      for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
                    }
                    name_cls[[l]] <- name
                  }
                }
              }
            }
          }
        }
      }
      # cat("unlist(name_cls) \n")
      # print(unlist(name_cls))
    }
  }
  Proba_Treat$Cls <- unlist(name_cls)
  # cat("- Modification des noms de chaque cluster dans le dataframe \n")
  new_names <- list() ; length(new_names) <- nrow(data_cluster)
  for (i in 1:K) {
    new_names[which(data_cluster$Cls == i)] <- name_cls[[i]]
  }
  data4 <- data_cluster
  data4$Cls <- unlist(new_names)
  cat("- Sauvegarde des résultats. \n")
  ## Dataframe avec l'écart moyen de survie par classe
  # cat("- Calcul de l'écart moyen de survie par classe \n")
  mean_ecart <- data.frame(time = rep(Time,K),
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
  mean_ecart$surv <- unlist(my_mean)
  ## Sauvegarde
  # cat("- Sauvegarde des résultats \n")
  path <- paste0(results_path,"/Proba_Treat.RData")
  save(Proba_Treat, file = path)
  save(Proba_curves, Proba_curves_bis, Proba_Treat, mean_ecart, file = Proba_path)
  return(data4)
}

Pb_Bsplines_GMM <- function(Data){
  ## Chargement des objets necessaires
  # cat("- Importation des données nécessaires. \n")
  load(Clustering_path)
  threshold_deviance <- 0.01
  threshold_deviance2 <- 0.15
  ## Creation cond
  my_Wells <- unique(Data$Wells)
  Cond <- list() ; length(Cond) <- length(my_Wells)
  for (i in 1:length(my_Wells)) {
    condition <- Data$cond[which(Data$Wells == my_Wells[i])][1]
    test <- strsplit(as.character(condition),split = "\n")[1]
    Cond[i] <- test[[1]]
  }
  ## Proba_curves
  # cat("- Calcul des probabilités qu'une courbe soit dans un cluster. \n")
  Proba_curves <- data.frame(Pb=Gmm$z, cond = unlist(Cond))
  Proba_curves_bis <- data.frame(Gmm$z)
  ## Proba Treatecules
  n_Treat <- length(treat_group)
  Proba_Treat <- data.frame(Pb = rep(NA,n_Treat*K), Cls = rep(NA,n_Treat*K), Treat = rep(NA,n_Treat*K))
  #Pour la 1ère Treatécule
  # cat("- Calcul des probabilités que la 1ère Treatécule soit dans chaque cluster. \n")
  Proba <- data.frame(Pb = rep(NA,K), Cls = rep(NA,K), Treat = rep(treat_names[[1]],K))
  for (j in 1:K){
    dat <- Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[1]),]
    if (length(treat_group[[1]])>1) {
      for (k in 2:length(unlist(treat_group[[1]]))) {
        dat <- rbind(dat,Proba_curves[which(Proba_curves$cond==unlist(treat_group[[1]])[k]),])
      }
    }
    Proba$Pb[j] <- 1/(dim(dat)[1])*sum(dat[,j])
    Proba$Cls[j] <- j
    Proba$Treat[j] <- treat_names[[1]]
  }
  Proba_Treat <- Proba
  ## Les autres Treatécules
  # cat("- Calcul des probabilités que chaque Treatécule soit dans chaque cluster. \n")
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
  ## Noms clusters
  name_cls <- list() ; length(name_cls) <- K
  sign <- matrix(NA, nrow=K, ncol = n_cut)
  
  
  
  moyennes <- matrix(NA,ncol=n_cut,nrow=K)
  for (i in 1:K) {
    d <- Data[which(Data$Cls==i),]
    for (j in 1:n_cut) {
      t_m <- ((j-1)*time_max)/n_cut
      t_p <- ((j*time_max)/n_cut)
      moyennes[i,j] <- mean(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
    }
  }
  
  number_sign <- floor(moyennes / threshold_deviance2)
  
  symbol_sign <- sign(number_sign)
  symbol_sign[symbol_sign == -1] <- "-"
  symbol_sign[symbol_sign ==  0] <- "="
  symbol_sign[symbol_sign ==  1] <- "+"
  
  name_cls <- rep(NA,K)
  name_cls[1] <- "Control"
  for(i in 2:K){
    name_i <- NULL
    for(j in 1:n_cut){
      if(number_sign[i,j] > 0){
        name_i <- paste(name_i,
                        paste0(rep(symbol_sign[i,j],each=number_sign[i,j]),collapse = ""),
                        sep="/")
      }else{
        
        name_i <- paste(name_i,"=",sep="/")
      }
    }
    name_i <- gsub(x = name_i,pattern = "^/",replacement = "")
    name_cls[i] <-  name_i
  } 
  
  
  name_cls[-1] <- paste("G",1:(K-1)," ",name_cls[-1],sep = "")
  
  
  
  # 
  # 
  # 
  # for (i in 1:K) {
  #   # cat(cat("- Caractérisation du",i),"ème cluster. \n",sep="")
  #   d <- Data[which(Data$Cls==i),]
  #   for (j in 1:n_cut) {
  #     t_m <- ((j-1)*time_max)/n_cut
  #     t_p <- ((j*time_max)/n_cut)
  #     Mean <- median(data.frame(d$surv[which(d$time>=t_m & d$time < t_p)])[,1])
  #     if (Mean <threshold_deviance & Mean>=-threshold_deviance) sign[i,j] <- "="
  #     if (Mean >= threshold_deviance) sign[i,j] <- "+"
  #     if (Mean < -threshold_deviance) sign[i,j] <- "-"
  #   }
  #   if (Proba_Treat$Pb[which(Proba_Treat$Cls == i & Proba_Treat$Treat == "Temoin")]>0.75) {
  #     name_cls[i] <- "Témoin"
  #   }
  #   else {
  #     name <- sign[i,1]
  #     for (k in 2:n_cut){
  #       name <- paste(name,sign[i,k],sep=",")
  #     }
  #     name_cls[i] <- name
  #   }
  #   
  #   
  #   
  #   if (i>1){
  #     # cat("- Vérification qu'il n'y ait pas deux fois le même nom. \n")
  #     # cat("unlist(name_cls) \n")
  #     # print(unlist(name_cls))
  #     # cat("Plusieurs clusters avec le même nom ? \n")
  #     # print(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)])))
  #     while(Same_name(name_cls[[i]],unlist(name_cls[1:(i-1)]))==TRUE){
  #       for (j in 2:i){
  #         for (l in 1:(j-1)){
  #           if (name_cls[[j]]==name_cls[[l]]){
  #             ## Valeur absolue de l'ecart entre les deux médianes
  #             Tab <- list(); length(Tab) <- n_cut
  #             for (m in 1:n_cut){
  #               dj <- Data[which(Data$Cls==j),]
  #               dl <- Data[which(Data$Cls==l),]
  #               t_m <- ((m-1)*time_max)/n_cut
  #               t_p <- (m*time_max)/n_cut
  #               Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
  #               Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
  #               Tab[m] <- abs(Meanj-Meanl)
  #             }
  #             ## Localisation du plus grand ecart entre les deux groupes
  #             k <- which.max(unlist(Tab))
  #             t_m <- ((k-1)*time_max)/n_cut
  #             t_p <- (k*time_max)/n_cut
  #             Meanj <- median(data.frame(dj$surv[which(dj$time>=t_m & dj$time < t_p)])[,1])
  #             Meanl <- median(data.frame(dl$surv[which(dl$time>=t_m & dl$time < t_p)])[,1])
  #             if (Meanj>=threshold_deviance) {
  #               if (Meanl>=Meanj){
  #                 name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                 if (k==1) {name <- paste(name,sign[j,1])
  #                 for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                 else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   name <- paste0(name,"+")
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                 }
  #                 name_cls[[j]] <- name
  #               } else {
  #                 name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                 if (k==1) {name <- paste(name,sign[l,1])
  #                 for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                 else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   name <- paste0(name,"+")
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                 }
  #                 name_cls[[l]] <- name
  #               }
  #             }
  #             if (Meanj<= -threshold_deviance) {
  #               if (Meanj<=Meanl){
  #                 name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                 if (k==1) {name <- paste(name,sign[j,1])
  #                 for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                 else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   name <- paste0(name,"-")
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                 }
  #                 name_cls[[j]] <- name
  #               } 
  #               else {
  #                 name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                 if (k==1) {name <- paste(name,sign[l,1])
  #                 for (iter in 2:n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                 else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   name <- paste0(name,"-")
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                 }
  #                 name_cls[[l]] <- name
  #               }
  #             }
  #             if ((Meanj < threshold_deviance) & (Meanj > -threshold_deviance)){
  #               if (Meanj >= 0 & Meanl >= 0) {
  #                 if (Meanj >= Meanl) {
  #                   name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[j,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"+")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[j]] <- name
  #                 } else {
  #                   name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[l,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"+")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[l]] <- name
  #                 }
  #               }
  #               if (Meanj <= 0 & Meanl <= 0) {
  #                 if (Meanj <= Meanl) {
  #                   name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[j,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"-")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[j]] <- name
  #                 } else {
  #                   name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[l,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"-")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[l]] <- name
  #                 }
  #               }
  #               if (Meanj >= 0 & Meanl <= 0) {
  #                 if (Meanj >= abs(Meanl)) {
  #                   name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[j,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"+")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[j]] <- name
  #                 } else {
  #                   name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[l,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"-")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[l]] <- name
  #                 }
  #               }
  #               if (Meanj <= 0 & Meanl >= 0) {
  #                 if (abs(Meanj) >= Meanl) {
  #                   name <- unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[j,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"-")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[j]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[j]] <- name
  #                 } else {
  #                   name <- unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[1]))
  #                   if (k==1) {name <- paste(name,sign[l,1])
  #                   for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}}
  #                   else {for (iter in 2:k){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                     name <- paste0(name,"+")
  #                     for (iter in (k+1):n_cut){name <- paste(name,unlist(lapply(strsplit(name_cls[[l]],','),function (v) v[iter])),sep=",")}
  #                   }
  #                   name_cls[[l]] <- name
  #                 }
  #               }
  #             }
  #           }
  #         }
  #       }
  #     }
  #     # cat("unlist(name_cls) \n")
  #     # print(unlist(name_cls))
  #   }
  #   
  #   
  #   
  #   
  #   
  # }
  # 
  # 
  # 
  
  
  Proba_Treat$Cls <- unlist(name_cls)
  # cat("- Modification des noms de chaque cluster dans le dataframe \n")
  new_names <- list() ; length(new_names) <- nrow(data_cluster)
  for (i in 1:K) {
    new_names[which(data_cluster$Cls == i)] <- name_cls[[i]]
  }
  data4 <- data_cluster
  data4$Cls <- unlist(new_names)
  ## Dataframe avec l'écart moyen de survie par classe
  # cat("- Calcul de l'écart moyen de survie par classe \n")
  mean_ecart <- data.frame(time = rep(Time,K),
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
  mean_ecart$surv <- unlist(my_mean)
  ## Sauvegarde
  # cat("- Sauvegarde des résultats. \n")
  path <- paste0(results_path,"/Proba_Treat.RData")
  save(Proba_Treat, file = path)
  save(Proba_curves, Proba_curves_bis, Proba_Treat, mean_ecart, file = Proba_path)
  return(data4)
}

#### Graphical outputs ----

plot_curves <- function(Data) {
  p <- ggplot(Data,aes(time,surv,col=Treat, group = Wells)) + 
    geom_line(size = 0.4, alpha = 0.4)+ 
    theme_pubclean(base_size = 15)
  return(p)
}
plot_vs <- function(){
  ## Graphique avce ConSpline
  n <- length(unique(data_transformed$Treat))
  
  ## Graphique avec Kaplan-Meier 
  p_km <- plot_curves(pretreat_KM(data)) + xlab("time en jours") + 
    ylab("") + scale_color_manual(values = rainbow(n)) + 
    theme(legend.position = "none") + ggtitle("Kaplan-Meier")
  
  ## Graphique avec conspline
  p <- plot_curves(pretreat_conspline(data)) + xlab("time en jours") + 
    ylab("Probabilité de survie") + scale_color_manual(values = rainbow(n)) + 
    labs(col="Treatécules") + theme(legend.position = "bottom") + ggtitle("Conspline")
  
  ## Graphique avec cobs
  p2 <- plot_curves(data_pretreated) + xlab("time en jours") + ylab("Probabilité de survie") + 
    scale_color_manual(values = rainbow(n)) + labs(col="Treatécules") + ggtitle("Cobs") + 
    theme(legend.position = "bottom")
  
  Path1 <- paste0(img_path,"/Conspline_KM.pdf")
  Plot1 <- ggpubr::ggarrange(p, p_km, ncol = 2, nrow = 1, common.legend = TRUE)
  suppressMessages(ggsave(Path1, plot = Plot1, device = "pdf"))
  
  Path2 <- paste0(img_path,"/Cobs_KM.pdf")
  Plot2 <- ggpubr::ggarrange(p2, p_km, ncol = 2, nrow = 1, common.legend = TRUE)
  suppressMessages(ggsave(Path2, plot = Plot2, device = "pdf"))
  
  Path_conspline <- paste0(results_path,"/MSE_conspline.RData")
  Path_cobs <- paste0(results_path,"/MSE_cobs.RData")
  load(Path_conspline)
  load(Path_cobs)
  
  ## Data_frame
  MSE_conspline$splines <- "Conspline"
  MSE_conspline$Mean <- mean(MSE_conspline$MSE)
  MSE_cobs$splines <- "Cobs"
  MSE_cobs$Mean <- mean(MSE_cobs$MSE)
  MSE_data <- rbind(MSE_conspline, MSE_cobs)
  MSE_data$splines <- factor(MSE_data$splines)
  
  ## Graphique
  plot_path <- paste0(img_path,"/MSE_conspline.pdf")
  plot_MSE <- ggplot(MSE_data, aes(Curves, MSE, col = splines)) + geom_point() + 
    geom_hline(aes(yintercept = Mean, col = splines), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle("MSE pour cobs et conspline")
  suppressMessages(ggsave(plot_path, plot = plot_MSE, device = "pdf"))
}



