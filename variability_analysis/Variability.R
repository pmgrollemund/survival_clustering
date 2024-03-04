################################################################################
################### Analyse the error count variability
######## Elise Comte - Paul-Marie Grollemund
######## 2021-07-19  - 2023-05-23
################################################################################

#### Clean up ----
rm(list=ls())

#### Define the appropriate working directory ----
setwd(".")

#### Required packages ----
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
library(ggpubr)
library(viridis)

#### Required functions ----
completion_NA <- function(df){
  number_NA <- sum(is.na(df$freq))
  while(number_NA>0){
    number_NA_debut <- sum(is.na(df$freq))
    
    index_na <- which(is.na(df$freq))
    for(i in sample(index_na)){
      current <- df[i,]
      
      tmp <- df[
        df$moyenne %in% (df$moyenne[i]+c(-1,0,1)) &
          df$ecart %in% (df$ecart[i]+c(-1,0,1)),
      ]
      tmp <- tmp[!is.na(tmp$freq),]
      
      if(nrow(tmp) >0){
        distances <- abs(current$moyenne - tmp$moyenne) + abs(current$ecart - tmp$ecart) 
        poids <- rep(0, length(distances))
        for(j in 1:length(poids)){
          poids[j] <- df_poids$poids[df_poids$distance == distances[j]]  
        }
        
        df$freq[i] <- crossprod(tmp$freq,poids)
      }
      
    }
    number_NA <- sum(is.na(df$freq))
    if(number_NA_debut == number_NA) break
  }
  return(df)
}
lisse <- function(df){
  df0 <- df
  
  for(i in sample(1:nrow(df))){
    current <- df[i,]
    
    tmp <- df0[
      df$moyenne %in% (df$moyenne[i]+c(-1,0,1)) &
        df$ecart %in% (df$ecart[i]+c(-1,0,1)),
    ]
    
    distances <- abs(current$moyenne - tmp$moyenne) + abs(current$ecart - tmp$ecart) 
    poids <- rep(0, length(distances))
    for(j in 1:length(poids)){
      poids[j] <- df_poids$poids[df_poids$distance == distances[j]]  
    }
    
    df$freq[i] <- crossprod(tmp$freq,poids)
    
  }
  
  moyennes <- unique(df$moyenne)
  for(j in 1:length(moyennes)){
    df$freq[df$moyenne == moyennes[j]] <- 
      df$freq[df$moyenne == moyennes[j]] / 
      sum(df$freq[df$moyenne == moyennes[j]])
  }
  
  return(df)
}


#### Create the output dir ----
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path("results",paste("ANOVA_run",timestamp,sep="_"))

if(!dir.exists("results")) 
  dir.create("results")

dir.create(output_dir)



#### Data importation ----
# Import ----

# Preteat ----


#### Compute means and deviations ----
# Estimate the number of nematode in the pit ----

# Estimate error counts ----


#### Graphical outputs ----



#### Chargement des données ----
setwd("/home/pmgrolle/Documents/MCF/Recherche/UMRF/Elegans/Survival_classification/code/analyse_variabilite/")

data <- readODS::read_ods(path = "../../data/Data_test_Elise/ComptagesX10.ods")
data <- data[,-(1:2)]

moyenne <- round(apply(data,1,mean),0)

df_data <- data.frame(
  mesure = as.vector(t(as.matrix(data))),
  moyenne = rep(moyenne,each=ncol(data)),
  jours = rep(c(0,1,4,5,6,7,8,12,13),each=3*ncol(data)),
  puit = rep(as.factor(1:3),each=ncol(data))
)

#### Chargement des données de la seconde acquisition de données ----
# PMG 2021/09/14
data_tmp <- readODS::read_ods(path = "../../data/Data_test_Elise/ComptagesX10bis.ods")
data_tmp <- data_tmp[,-(1:2)]

data_tmp_MB <- data_tmp[,1:10]
data_tmp_MB <- data_tmp_MB[-c(13:15,19:21),]
data_tmp_PV <- data_tmp[,14:23]

moyenne_MB <- round(apply(data_tmp_MB,1,mean),0)
moyenne_MB <- moyenne_MB[!is.na(moyenne_MB)]
moyenne_PV <- round(apply(data_tmp_PV,1,mean),0)

df_data <- rbind(
  df_data, 
  data.frame(
    mesure = as.vector(t(as.matrix(data_tmp_MB))),
    moyenne = rep(moyenne_MB,each=ncol(data_tmp_MB)),
    jours = rep(c(0,1,4,5,7,11),each=3*ncol(data_tmp_MB)),
    puit = rep(as.factor(1:3),each=ncol(data_tmp_MB)) 
  ),
  data.frame(
    mesure = as.vector(t(as.matrix(data_tmp_PV))),
    moyenne = rep(moyenne_PV,each=ncol(data_tmp_PV)),
    jours = rep(c(0,1,4,5,6,7,8,11),each=3*ncol(data_tmp_PV)),
    puit = rep(as.factor(1:3),each=ncol(data_tmp_PV))
  )
)


df_data$ecart <- df_data$mesure - df_data$moyenne

ecart_min <- min(df_data$ecart)
ecart_max <- max(df_data$ecart)

ecarts <- ecart_min:ecart_max
moyennes <- min(df_data$moyenne):max(df_data$moyenne)

distribution_ecart <- expand.grid(
  moyenne=moyennes,ecart=ecarts
)
distribution_ecart$freq <- rep(0,nrow(distribution_ecart))

for(i in 1:nrow(distribution_ecart)){
  index <- which(df_data$moyenne == distribution_ecart$moyenne[i] & 
          df_data$ecart == distribution_ecart$ecart[i]
        )
  distribution_ecart$freq[i] <- length(index)
}

for(moyenne in moyennes){
  distribution_ecart$freq[distribution_ecart$moyenne == moyenne] <- 
    distribution_ecart$freq[distribution_ecart$moyenne == moyenne] / 
    sum(distribution_ecart$freq[distribution_ecart$moyenne == moyenne])
}

distribution_ecart$freq[distribution_ecart$freq %in% c(0,NaN)] <- NA


ggplot(distribution_ecart, aes(moyenne, ecart, fill= freq)) + 
  geom_tile() + 
  theme_pubclean(base_size = 15)

df_poids <- data.frame(
  distance = 0:5
)
df_poids$poids <- exp(-df_poids$distance)

distribution_ecart_lisse <- distribution_ecart

df <- distribution_ecart
distribution_ecart_complete <- completion_NA(distribution_ecart)
distribution_ecart_lisse <- lisse(distribution_ecart_complete)
distribution_ecart_lisse <- lisse(distribution_ecart_lisse)



ggplot(distribution_ecart_lisse, aes(moyenne, ecart, fill= freq)) + 
geom_tile() + scale_fill_viridis()+ 
  theme_pubclean(base_size = 15)




save(distribution_ecart_lisse,file = "../classification/src/params_simul_PMG.RData")







sds <- rep(NA,length(unique(df_data$moyenne)))
for(i in 1:length(sds)){
  sds[i] <- sd(df_data$ecart[df_data$moyenne == sort(unique(df_data$moyenne)[i])])
}
plot(sort(unique(df_data$moyenne)),sds)


#### Quantiles des écarts ----
df_data$ecart <- df_data$mesure - df_data$moyenne
df_data$quan_min <- NA
df_data$quan_max <- NA

for (i in 1:length(unique(df_data$moyenne))) {
  Ind <- which(df_data$moyenne == unique(df_data$moyenne)[i])
  Kant <- quantile(df_data$ecart[Ind],probs = c(0.05,0.95))
  df_data$quan_min[Ind] <- Kant[[1]]
  df_data$quan_max[Ind] <- Kant[[2]]
}

ggplot(df_data) + geom_boxplot(aes(x=as.factor(moyenne),y=ecart))


fun_exp <- function(x) return(2.5+exp((x/3)-9.7))
line_exp <- data.frame(
  X = seq(0,35,1),
  Y_sup = fun_exp(seq(0,35,1)),
  Y_inf = -fun_exp(seq(0,35,1))
)

## Nombres de points hors de l'intervalle d?fini
count_point <- 0
for (i in 1:nrow(df_data)) {
  ind <- which(line_exp$X == df_data$moyenne[i])
  # cat("ind \n")
  # print(ind)
  if ((df_data$ecart[i] < line_exp$Y_inf[ind]) || 
      (df_data$ecart[i] > line_exp$Y_sup[ind])) {
    count_point <- count_point + 1
  }
}
## Pourcentage
Pourcent <- count_point/nrow(df_data)*100

ec <- ggplot(df_data,aes(x=moyenne,y=ecart)) + geom_point() + 
  geom_jitter(width = 0.2) + 
  geom_smooth(method = "loess", formula = y~x) +
  geom_ribbon(aes(ymin=quan_min, ymax=quan_max), fill="blue", alpha=0.2) +
  geom_line(data = line_exp, aes(X,Y_sup)) + 
  geom_line(data = line_exp, aes(X,Y_inf))
ec

ggsave(file="test_EC/Img/ecarts.pdf", plot = ec, device = "pdf")

## Test estimation param?tres lois normales
n_rep <- 20
data_test <- data.frame(nvers = rep(line_exp$X, n_rep), 
                        ecart = rep(NA, n_rep*nrow(line_exp)))

count_ext_point <- 0
for (i in 1:nrow(data_test)) {
  ind <- which(line_exp$X == data_test$nvers[i])
  # sd <- 1/(1.96^2)*abs(line_exp$Y_sup[ind]-line_exp$Y_inf[ind])
  sd <- 1/4*abs(line_exp$Y_sup[ind]-line_exp$Y_inf[ind])
  data_test$ecart[i] <- rnorm(1, mean = 0, sd = sd)
  if ((data_test$ecart[i]<line_exp$Y_inf[ind])|| (data_test$ecart[i]>line_exp$Y_sup[ind])) {
    count_ext_point <- count_ext_point + 1
  }
}
## Pourcentage de point qui ne sont pas dans l'intervalle
pourcent_ext <- count_ext_point/nrow(data_test)*100


test_plot <- ggplot(data = data_test, aes(x=nvers,y=ecart)) + geom_point() + 
  geom_jitter(width = 0.2) + 
  # geom_ribbon(data = line_exp, aes(ymin=Y_inf, ymax=Y_sup), fill="blue", alpha=0.2) +
  geom_line(data = line_exp, aes(X,Y_sup), col = 'red') + 
  geom_line(data = line_exp, aes(X,Y_inf), col = 'red')
test_plot 

ggsave(file="test_EC/Img/simul_ecarts.pdf", plot = test_plot, device = "pdf")
  
ecart_type <- data.frame(
  nvers = line_exp$X,
  sd = 1/4*abs(line_exp$Y_sup-line_exp$Y_inf)
)
save(ecart_type, file = "./test_EC/Res/params_simul.RData")

#### Estimation de la proba des écarts de comptages ---- 
ecart_moyenne <- matrix(NA,ncol=diff(range(df_data$ecart))+1,nrow=diff(range(df_data$moyenne))+1)
rownames(ecart_moyenne) <- max(df_data$moyenne):min(df_data$moyenne)
colnames(ecart_moyenne) <- min(df_data$ecart):max(df_data$ecart)

for(i in 1:nrow(df_data)){
  index_ligne <- which(as.numeric(rownames(ecart_moyenne)) == df_data$moyenne[i])
  index_colonne <- which(as.numeric(colnames(ecart_moyenne)) == df_data$ecart[i])
  
  if(is.na(ecart_moyenne[index_ligne,index_colonne])){
    ecart_moyenne[index_ligne,index_colonne] <- 1
  } else {
    ecart_moyenne[index_ligne,index_colonne] <- ecart_moyenne[index_ligne,index_colonne] + 1 
  }
}

tmp <- ecart_moyenne[,which(as.numeric(colnames(ecart_moyenne))<5 & as.numeric(colnames(ecart_moyenne))>-6)]

## Test calcul proba 
Probas <- tmp
for (i in 1:nrow(tmp)) {
  n <- sum(tmp[i,], na.rm = T)
  for (j in 1:ncol(tmp)) {
    Probas[i,j] <- tmp[i,j]/n
  }
}
Probas <- round(Probas, digits = 3)

## Moyennes des probas
PJC <- apply(Probas, 2, mean, na.rm = T)
scale(PJC, center = F, scale =T)
Probas_mean <- PJC/sum(PJC)
Proba_mean <- data.frame(X = c(-5,-4,-3,-2,-1,0,1,2,3,4), Probas = Probas_mean)
plot(c(-5,-4,-3,-2,-1,0,1,2,3,4),Probas_mean, xlab = "Écart", ylab = "Probabilité", type = 'l')

#### Graphiques densités ----

## Remplissage des NAs
tmp2 <- tmp
ncols <- ncol(tmp2)
for (i in 1:nrow(tmp2)) {
  if (!sum(is.na(tmp2[i,])) == ncols) {
    if (is.na(tmp2[i,1])) {
      for (j in 1:(which(!is.na(tmp2[i,]))[[1]] - 1)) {
        tmp2[i,j] <- 0
      }
    }
    if (is.na(tmp2[i,ncols])) {
      for (j in (tail(which(!is.na(tmp2[i,])),1) + 1):ncols) {
        tmp2[i,j] <- 0
      }
    }
  }
}

data <- tmp2[rowSums(is.na(tmp2)) == 0, ]


Proba <- data
for (i in 1:nrow(data)) {
  n <- sum(data[i,], na.rm = T)
  for (j in 1:ncol(data)) {
    Proba[i,j] <- data[i,j]/n
  }
}
Proba <- round(Proba, digits = 3)



## freq
freq <- NULL
for (i in 1:nrow(Proba)) {
  freq <- cbind(freq, t(Proba[i,]))
}
freq <- t(freq)

## Nvers
nvers <- list() ; length(nvers) <- nrow(Proba)
for (i in 1:nrow(Proba)) {
  nvers[[i]] <- rep(as.numeric(rownames(Proba))[i],ncol(tmp))
}


## Data
Data <- data.frame(Freq = freq, 
                   Ecart = rep(c(-5,-4,-3,-2,-1,0,1,2,3,4),nrow(Proba)),
                   Nvers = factor(unlist(nvers)))


ploplot <- ggplot(Data, aes(Ecart, Freq)) +
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(.~Nvers) + geom_density(stat="identity", aes(y = Freq))

ploplot

## Castex
Castex <- NULL
for (i in 1:nrow(Proba)) {
  Castex <- rbind(Castex, df_data[which(df_data$moyenne == 
                                          as.numeric(rownames(Proba))[i]),])
}
Castex$Mean_ecart <- rep(unlist(Probas_mean), nrow(Castex)/10)
Castex <- Castex[(which(Castex$ecart < 5 & Castex$ecart > -6)),]


Castex$moyenne <- factor(Castex$moyenne)
ploplot2 <- ggplot(data = Castex, aes(x = ecart, col = moyenne, group = moyenne)) + 
  geom_density(aes(x = ecart, col = moyenne, group = moyenne)) +
  scale_color_manual(values = viridis::plasma(nrow(Proba))) + xlim(-5,4) 
# + geom_point(data = Proba_mean, aes(x = X, y = Probas))

ploplot2


#### Tests de Shapiro ----

data <- tmp[rowSums(is.na(tmp)) < ncol(tmp), ]
Pval <- list() ; length(Pval) <- nrow(data)
for (i in 1:nrow(data)) {
  Pval[i] <- shapiro.test(data[i,])$p.value
}

PAC <- data[which(Pval>0.05),]

## Transformation en probas
Prob <- PAC
for (i in 1:nrow(PAC)) {
  n <- sum(PAC[i,], na.rm = T)
  for (j in 1:ncol(PAC)) {
    Prob[i,j] <- PAC[i,j]/n
  }
}
Prob <- round(Prob, digits = 3)

## Normalisation
my_sum <- apply(PAC, 1, sum, na.rm = T)
my_sum <- my_sum/10
PAC_bis <- PAC
for (i in 1:nrow(PAC)) {
  PAC_bis[i,] <- PAC[i,]/my_sum[i]
}

#### Paramètres des densités ----
Params <- data.frame(Mean = rep(NA, nrow(tmp)), Sd = rep(NA, nrow(tmp)), 
                     row.names = as.numeric(rownames((tmp))))
## Fonction qui calcule l'esperance
Hope <- function(Vec, Col, n) {
  Esp <- 0
  for (i in 1:length(Vec)) {
    if (!is.na(Vec[i])) {
      Esp <- Esp + (Vec[[i]] * Col[i])
    }
  }
  return(Esp/n)
  # return(c(floor(Esp/n), ceiling(Esp/n)))
}

COL <- as.numeric(colnames(PAC))
n <- ncol(PAC)
for (i in 1:nrow(PAC)) {
  Mu <- Hope(PAC_bis[i,],COL,n)
  # Mu <- sample(Hope(Prob[i,],COL,n),1,prob = c(0.5,0.5))
  Ind <- as.numeric(rownames(PAC)[i])
  Params$Mean[as.numeric(rownames(Params)) == Ind] <- Mu
  Params$Sd[as.numeric(rownames(Params)) == Ind] <- sd(PAC_bis[i,], na.rm = T)
}

ind_NA <- which(is.na(Params$Mean))
not_NA <- which(!is.na(Params$Mean))
for (i in 1:(length(ind_NA)-1)) {
  if (is.na(Params$Mean[ind_NA[i]])) {
    if ((ind_NA[i+1] - ind_NA[i]) != 1) {
      Params$Mean[ind_NA[i]] <- Params$Mean[ind_NA[i]-1] + 
        abs(Params$Mean[ind_NA[i]+1]-Params$Mean[ind_NA[i]-1])/2
      Params$Sd[ind_NA[i]] <- Params$Sd[ind_NA[i]-1] + 
        abs(Params$Sd[ind_NA[i]+1]-Params$Sd[ind_NA[i]-1])/2
    } else {
      Ind <- not_NA[which(not_NA > ind_NA[i])[1]]-1
      for (j in 1:(Ind-ind_NA[i]+1)) {
        Params$Mean[ind_NA[i] +j-1] <- Params$Mean[ind_NA[i]-1] + 
          (j*abs(Params$Mean[Ind+1]-Params$Mean[ind_NA[i]-1]))/(Ind-ind_NA[i]+2)
        Params$Sd[ind_NA[i] +j-1] <- Params$Sd[ind_NA[i]-1] + 
          (j*(Params$Sd[Ind+1]-Params$Sd[ind_NA[i]-1]))/(Ind-ind_NA[i]+2)
      }
    }
  }
}

## Densités obtenues
x = seq(-6,5,length.out = 100)
plot(x,dnorm(x,Params$Mean[1],Params$Sd[1]), type='l', ylim = c(0,0.6))
for (i in 2:nrow(Params)) {
  lines(x,dnorm(x,Params$Mean[i],Params$Sd[i]))
}

## Graphique des moyennes et écart-type
ggplot(Params) + geom_point(aes(as.numeric(rownames(Params)), Sd, col = 
                                  as.numeric(rownames(Params)))) + 
  geom_point(aes(as.numeric(rownames(Params)), Mean, col = 
                   as.numeric(rownames(Params)))) + 
  scale_x_reverse() + geom_smooth(aes(as.numeric(rownames(Params)), Sd, col = 
                                        as.numeric(rownames(Params)))) + 
  geom_smooth(aes(as.numeric(rownames(Params)), Mean, col = 
                    as.numeric(rownames(Params))))

## Graphiques des densités vraies données

## freq
freq <- NULL
for (i in 1:nrow(Prob)) {
  freq <- cbind(freq, t(Prob[i,]))
}
freq <- t(freq)

## Nvers
nvers <- list() ; length(nvers) <- nrow(Prob)
for (i in 1:nrow(Prob)) {
  nvers[[i]] <- rep(as.numeric(rownames(Prob))[i],ncol(tmp))
}

## Params
moy <- list() ; length(moy) <- nrow(Prob)
ecart_type <- list() ; length(ecart_type) <- nrow(Prob)
for (i in 1:length(moy)) {
  moy[i] <- Params$Mean[as.numeric(rownames(Params))==as.numeric(rownames(Prob))[i]]
  ecart_type[i] <- Params$Sd[as.numeric(rownames(Params))==as.numeric(rownames(Prob))[i]]
}


## Data
Densi <- data.frame(Freq = freq, 
                    Ecart = rep(c(-5,-4,-3,-2,-1,0,1,2,3,4),nrow(Prob)),
                    Nvers = factor(unlist(nvers)),
                    Mean = rep(unlist(moy), each = ncol(tmp)),
                    Sd = rep(unlist(ecart_type), each = ncol(tmp)))

ploplot <- ggplot(Densi, aes(Ecart, Freq)) +
  geom_bar(stat="identity",position=position_dodge()) +
  facet_grid(.~Nvers) + 
  geom_density(stat="identity", aes(y = Freq), col = "red") + 
  geom_line(aes(Ecart, dnorm(Ecart, Mean, Sd)), col = "blue")

ggsave(file = "test_EC/Img/Comparaison_densité.pdf", plot = ploplot, 
       device = "pdf", width = 14, height = 9)

## Graphique avec toutes les densités 
N <- 100
X = seq(-5, 4, length.out = N)
Density <- list() ; length(Density) <- nrow(Params)
for (i in 1:nrow(Params)) {
  Density[[i]] <- dnorm(X, Params$Mean[i], Params$Sd[i])
}
Dens <- data.frame(Nvers = rep(as.numeric(rownames(Params)),each = N), 
                   X = rep(X, nrow(Params)),
                   Dnorm = unlist(Density))
Dens$Nvers <- factor(Dens$Nvers)

all_densities <- ggplot(Dens, aes(x = X, y = Dnorm, col = Nvers)) + geom_line() + 
  scale_color_manual(values = rainbow(nrow(Params))) + xlab("Écart") + 
  ylab("Proba") + labs( col = "Nombre moyen de vers") + 
  ggtitle("Densités estimées par nombre de vers moyen")
ggsave(file = "test_EC/Img/densities_estimates.pdf", plot = all_densities, 
       device = "pdf")

save(df_data, Proba_mean, Params, file = "./test_EC/Res/Density.RData")
