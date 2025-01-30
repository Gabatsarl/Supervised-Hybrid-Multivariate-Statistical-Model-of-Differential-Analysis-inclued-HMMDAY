
## Load the packages

source("/Users/timakpo/Desktop/CC-Template/PostDoc/Redaction/library.R")

## Mis 
library(ape) ; library(vegan)
library(limma)
library("statmod")
library(edgeR)
library(ggplot2)

# Créer des données simulées
B<-20
LAMBDA<-matrix(0,nrow=B,ncol = 4)

Simulation_indepandant<-function(n,p, clayParam){
  
  ## Simulation des données 
  ## Creer une copule de Clayton of two variables
  claytonCop <- claytonCopula(param = clayParam, dim = p)
  
  ## Simulation of copula data
  sim_data <- rCopula(n, claytonCop)
  
  ## Transformation des données en données discrètes (par exemple, en Poisson)
  X <- apply(sim_data, 2, function(x) qpois(x, lambda = lambda0))
  gen_data<-
  X<-data.matrix(cbind(rep(1,n), X))
  
  Xbeta<-X%*% matrix(beta[1:(p+1)],ncol=1,nrow=p+1)
  groups <- factor(sample(c("A", "B"), n, replace = TRUE))
  
  ## Convertir la variable catégorielle en numérique (0 ou 1)
  groups <- as.numeric(cat_var) - 1  # "A" sera 0, "B" sera 1
  
  ## Simulation des erreurs et reconstitution de Y
  
  epsilon<-rnorm(n)
  Y<-Xbeta+beta[p+2]*condition+epsilon
  
  
  
  for (k in 1:B){
    sim_data_nb <- generate_nb_counts(nrow(gen_data), ncol(gen_data), 
                                  size_estimates, mu_estimates)
    m.D <- vegdist(t(sim_data_nb), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    design <- model.matrix(~groups+PCoA1+PCoA2,data)
    colnames(design)
    d <- DGEList(counts =sim_data_nb)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y,design=design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    # result<-topTable(out_ebayes,1)
    result_Limma <- topTable(out_ebayes, number = Inf)
    pval_pois<-result_Limma$P.Value # $adj.P.Val
    
    LAMBDA[k,1]<-median(qchisq(1-pval_pois,1))/qchisq(0.5,1)
    
    
    
    ## Avec une loi de Poisson
    
    
    # Créer des données simulées
    
    sim_data_poisson <- generate_poisson_counts(nrow(gen_data), ncol(gen_data),
                                                lambda_estimates,
                                                sqrt(lambda_estimates))
    m.D <- vegdist(t(sim_data_poisson), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    # Design 
    
    design <- model.matrix(~groups+PCoA1+PCoA2,data)
    colnames(design)
    d <- DGEList(counts =sim_data_poisson)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y,design=design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    # result<-topTable(out_ebayes,1)
    result_Limma <- topTable(out_ebayes, number = Inf)
    pval_pois<-result_Limma$P.Value # $adj.P.Val
    
    LAMBDA[k,2]<-median(qchisq(1-pval_pois,1))/qchisq(0.5,1)
    
    
    
    # ----- Distribution de Hurdle
    
    # Générer les données de présence/absence
    prob_zero <- 0.3  # Probabilité que le comptage soit zéro
    n_samples<-ncol(gen_data) ;   n_taxa<-nrow(gen_data)
    presence_absence <- rbinom(n_samples * n_taxa, 1, 1 - prob_zero)
    presence_absence_matrix <- matrix(presence_absence, nrow = n_samples, ncol = n_taxa)
    
    
    # Générer les comptages pour les valeurs non-nulles
    lambda_mean<-lambda_estimates
    sim_data_Hurdle <- matrix(0, nrow = n_samples, ncol = n_taxa)
    for (j in 1:n_taxa) {
      sim_data_Hurdle[, j] <- rpois(n_samples, lambda = lambda_mean) * 
        presence_absence_matrix[, i]
    }
    
    m.D <- vegdist(t(sim_data_Hurdle), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    # Design 
    
    design <- model.matrix(~groups+PCoA1+PCoA2,data)
    colnames(design)
    d <- DGEList(counts =sim_data_Hurdle)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y,design=design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    # result<-topTable(out_ebayes,1)
    result_Limma <- topTable(out_ebayes, number = Inf)
    pval_Hurdle<-result_Limma$P.Value # $adj.P.Val
    
    LAMBDA[i,3]<-median(qchisq(1-pval_Hurdle,1))/qchisq(0.5,1)
    
    
    
    ## Loi geometric GMH
    
    sim_data_GMH<-metaSPARSim(params2)
    
    m.D <- vegdist(t(sim_data_GMH$counts), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    ## ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design 
    
    design <- model.matrix(~groups+PCoA1+PCoA2,data)
    colnames(design)
    d <- DGEList(counts =sim_data_GMH$counts)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y,design=design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    ## result<-topTable(out_ebayes,1)
    result_Limma <- topTable(out_ebayes, number = Inf)
    pval_GMH<-result_Limma$P.Value # $adj.P.Val
    
    LAMBDA[k,4]<-median(qchisq(1-pval_GMH,1))/qchisq(0.5,1)
    
    
    
    
    
    ## Zero-Inflated Poisson (ZIP)
    
    ## Exemple de génération de données ZIP
    lambda_est <- mean(gen_data[gen_data > 0])  # Moyenne des comptages non-nuls
    pi_est <- sum(gen_data == 0) / length(gen_data)  # Proportion de zéros
    
    #pi_inflate <- 0.3
    #n <-n_samples*n_taxa
    # Créer la matrice de microbiome simulée
    sim_data_ZIP <- matrix(0, nrow = n_taxa, ncol = n_samples)
    
    ## Remplir la matrice avec les données ZIP
    #set.seed(123)  # Pour reproductibilité
    
    library(gamlss.dist)
    for (i in 1:n_taxa) {
      sim_data_ZIP[i, ] <- rZIP(n_samples, mu =mean(gen_data[i,])  , 
                                sigma = 0.01) # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
    }
    
    # Génération de zéros supplémentaires
    #zero_inflated <- rbinom(n, 1, pi_inflate)
    
    #sim_data_ZIP<-rpois(n, lambda_estimates) * (1 - zero_inflated)
    
    m.D <- vegdist(t(sim_data_ZIP), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    # Design 
    
    design <- model.matrix(~groups+PCoA1+PCoA2,data)
    colnames(design)
    d <- DGEList(counts =sim_data_pval_ZIP)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y,design=design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    # result<-topTable(out_ebayes,1)
    result_Limma <- topTable(out_ebayes, number = Inf)
    pval_ZIP<-result_Limma$P.Value # $adj.P.Val
    
    LAMBDA[k,5]<-median(qchisq(1-pval_ZIP,1))/qchisq(0.5,1)
    
    
    


  }
  
 return(LAMBDA) 
}


## Inclut y a modeliser







#########===================Simulation independante d'un echantillon =============##
##
##
##
##=========================================================================##




#rpois(n,lambda=)

#Y=g*epsilon<-mvtnorm::gg

##
beta_0<-0.4
beta_1<-0.3
beta<-c(beta_0,beta_1)
n<-100 # le nombre d'individus
p<-400 # le nombre de genes
sigma_eps <- diag(rep(1,p))
mean_eps<-rep(0,p)
epsilon <- rnorm(n=100, mean=0, sd=1)

## Conditions
Conditions <- c("A", "B")

## Nombre d'observations à simuler
n <- 100

## Simuler la variable catégorielle ordonnée
variable_ordonne <- factor(sample(Conditions, n, replace = TRUE),
                           ordered = FALSE, levels = c("A", "B"))
variable_ordonne<-relevel(variable_ordonne, ref = "A")

## Transformer en variable binaire (0 ou 1)
variable_ordonne <- as.numeric(variable_ordonne != "A")

#X=cbind(rep(1,n),variable_ordonne,epsilon)

y<-0





##===================== Simulation 1: Sans Condition =========================##
##
##
##
## ====== Methode 2 (hypothèse standard) =====================================##

library(MASS)  # pour mvrnorm

## Définir la matrice de corrélation
Sigma <- matrix(0.5, nrow = p, ncol = p)  # 0.5 de corrélation entre les colonnes
diag(Sigma) <- 1  # 1 sur la diagonale pour la corrélation parfaite avec soi-même

## Generer des variables latentes normales multivariées

latent_vars <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

## Transformer les variables latentes en variables de Poisson

X <- apply(latent_vars, 2, function(x) rpois(n, lambda = exp(x)))

## Simuler une variable catégorielle avec deux niveaux
cat_var <- factor(sample(c("A", "B"), n, replace = TRUE))

# Définir les coefficients pour chaque gène
beta <- runif(p, min = -2, max = 2)  # coefficients aléatoires entre -2 et 2

# Convertir la variable catégorielle en numérique (0 ou 1)
cat_var_numeric <- as.numeric(cat_var) - 1  # "A" sera 0, "B" sera 1

# Define the coef of "Condition"
gamma <- 1.5  # coef of conditions

# Générer y comme combinaison linéaire de X et de la variable catégorielle
y <- X %*% beta + gamma * cat_var_numeric +epsilon




## Ajustement du model

data=data.frame(groups=as.factor(cat_var_numeric))
m.D <- vegdist(X, "manhattan") ; 
result_GENERAL <- pcoa(m.D )

# ================================================##

data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]

library(doParallel)
library(limma)

Simul_pois<-function(y,X,data){
  
  for(i in 1: ncol(X)){
    data$Xg<-log10(X[,i]+1)
    data$Y<-y[,1]
    limma_design <- with(data, model.matrix(~ Xg+relevel(as.factor(groups),  ref = "0")+
                                              PCoA1+PCoA2))
    
    limma_fit <- lmFit(data$Y, limma_design)
    out_ebayes <- eBayes(limma_fit)
    
    if(inherits(limma_fit, "try-error")){ cat("NULL") }
    else {
      para <- limma_fit$coefficients ; 
      pval <- out_ebayes$p.value ;
      print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
    }
  }
  return(pval)
  
  
}


pval<-Simul_pois(y=y,X=X,data=data) #$PVALUE


## MB-GAN
# note: generate simulated data by Normal-to-Anything (NorTA) #

# libaray to implement NorTA #
library(SpiecEasi)









##===================== CAS 1 : Independance des données =========================##
##
##
##
## ====== -----------------------------------=====================================##



## Par la loi de Poisson
##
##
## Load the package copula
p<-3000
n<-100
lambda<-5

## Créer une copule de Clayton pour deux variables
claytonCop <- claytonCopula(2, dim = p)

## Simulation de données à partir de la copule
set.seed(123)
sim_data <- rCopula(n, claytonCop)

## Transformation des données en données discrètes (par exemple, en Poisson)
X <- apply(sim_data, 2, function(x) qpois(x, lambda = 5))

# Affichage de la matrice simulée
X


## simulation des donnes epsilon
epsilon<-rnorm(n)
beta<-runif(p,-1,1)
Y<-X%*%beta+epsilon
Y

## Mise en oeuvre du PCoA


# ================================================##

data<-data.frame(Y=Y,
                 PCoA1=result_GENERAL$vectors[,"Axis.1"],
                 PCoA2=result_GENERAL$vectors[,"Axis.2"])

## Ajustement du model

#library(doParallel)
#cl <- detectCores() %>% -1 %>% makeCluster  ###n'utiliser pas tous les coeurs (-1)
#registerDoParallel(cl)
library(limma) ; library(copula)
library(ape) ; library(vegan)
library(doParallel)

Ajustement_function<-function(y,X,condition){
  X<-X[,-1]
  m.D <- vegdist(X, "bray") ; 
  result_GENERAL <- pcoa(m.D )
  data_res<-data.frame(Y=y)
  data_res$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data_res$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  data_res$condition<-condition
  nG<-ncol(X)
  result<-matrix(0,ncol=2,nrow = nG)
  colnames(result)<-c("pvalue","logFC")
  
  for(i in 1:nG){
    data_res$Xg<-log10(X[,i]+1)
    
    limma_design <- with(data_res, model.matrix(~ Xg+relevel(as.factor(condition),
                                                        ref = "0")))
    limma_fit <- lmFit(data_res$Y, limma_design)
    out_ebayes <- eBayes(limma_fit)
    
    if(inherits(limma_fit, "try-error")){ cat("NULL") }
    else {
      para <- limma_fit$coefficients ; 
      pval <- out_ebayes$p.value ;
      #print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
      result[i,1]<-pval[2]
      result[i,2]<-topTable(out_ebayes,coef=2)$logFC
  }
  
  } 
  return(result)
  
}


a=Ajustement_function(y=Y,X=X,condition=)


####qqplot comparison---------------
library(cowplot)
library(plyr)
library(ggpubr)


myqqplot2 = function(pval){
  obs = -log10(sort(pval,decreasing=F))
  expect = -log10( 1:length(obs)/length(obs) )
  df = data.frame(Observed=obs, Expected=expect, pval=pval)
  
  # figure 1
  plot1 = ggscatter(df, x = "Expected", y = "Observed",
                    color = "#808080")+
    #palette ="lancet") +
    geom_abline(intercept = 0, slope = 1, color="red",lwd=1) +
    scale_x_continuous(name=expression(-Log[10](italic(p))~~Expected)) +
    scale_y_continuous(name=expression(-Log[10](italic(p))~~Observed),
                       breaks = c(1:9))+
    theme( text = element_text(family = "Times New Roman"),
           axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
           axis.text.y = element_text(size = 12, face = "bold"),
           axis.title.y = element_text(size = 16, face = "bold"),
           axis.title.x = element_text(size = 16, face = "bold"),
           plot.title = element_text(face = "bold", size = 16)
    )
  # figure 2
  return(plot1)
  
}


myqqplot2(pval=a[,1])



## Repetition des simulations B=200

## Loi de : POISSON

Simulation_Poisson<-function(B, # nombre de simulation,
                    n, # Taille d'échantillon
                    p, # nombre de metagenes
                    lambda0,
                    clayParam){ # parametre de la loi de Poisson
  beta<-rnorm(p+2)
  lambda1<-matrix(0,ncol=1,nrow = B)
  
  for( j in 1 : B) {
  
  ## Creer une copule de Clayton of two variables
  claytonCop <- claytonCopula(param = clayParam, dim = p)

  ## Simulation of copula data
  sim_data <- rCopula(n, claytonCop)
  
  ## Transformation des données en données discrètes (par exemple, en Poisson)
  X <- apply(sim_data, 2, function(x) qpois(x, lambda = lambda0))
  X<-data.matrix(cbind(rep(1,n), X))
  
  Xbeta<-X%*% matrix(beta[1:(p+1)],ncol=1,nrow=p+1)
  condition <- factor(sample(c("A", "B"), n, replace = TRUE))
  
  ## Convertir la variable catégorielle en numérique (0 ou 1)
  condition <- as.numeric(cat_var) - 1  # "A" sera 0, "B" sera 1
  
  ## Simulation des erreurs et reconstitution de Y
  
  epsilon<-rnorm(n)
  Y<-Xbeta+beta[p+2]*condition+epsilon
  
  ajust<-Ajustement_function(y=Y,X=X,condition=condition)
  chisq1 <- qchisq(1-ajust[,1],1) ; 
  lambda1[j,1] = median(chisq1)/qchisq(0.5,1)
  
  }
return(lambda1) 
  
}

## Apply Apply 

lambda_Poisson<-Simulation_Poisson(B=100, n=100, p=300, lambda0=5, clayParam=2)


## Boxplot des lambda   

boxplot(lambda_Poisson)

## Créer un dataframe
df <- data.frame(lambda = lambda_Poisson)

## Créer le boxplot avec ggplot2

ggplot(df, aes(x = "", y = lambda)) +
  geom_boxplot() +
  labs(title = "", y = "Lambda") +
  theme_pubr()+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",lwd=1.5)




## Binomiale negative

## Function ajustment Y~X**beta

Ajustement_function<-function(y,X,condition){
  X<-X[,-1]
  m.D <- vegdist(X, "euclidean") ; 
  result_GENERAL <- pcoa(m.D )
  data_res<-data.frame(Y=y)
  data_res$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data_res$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  data_res$PCoA3<-result_GENERAL$vectors[,"Axis.3"]
  
  data_res$condition<-condition
  nG<-ncol(X)
  result<-matrix(0,ncol=2,nrow = nG)
  colnames(result)<-c("pvalue","logFC")
  
  for(i in 1:nG){
    data_res$Xg<-log10(X[,i]+1)
    
    limma_design <- with(data_res, model.matrix(~ Xg+relevel(as.factor(condition),
                                                       ref = "0")+PCoA1+PCoA2+PCoA3))
    limma_fit <- lmFit(data_res$Y, limma_design)
    out_ebayes <- eBayes(limma_fit)
    
    if(inherits(limma_fit, "try-error")){ cat("NULL") }
    else {
      para <- limma_fit$coefficients ; 
      pval <- out_ebayes$p.value ;
      result[i,1]<-pval[2]
      result[i,2]<-topTable(out_ebayes,coef=2)$logFC
    }
    
  } 
  return(result)
  
}

##


simulation_NB=function(B, # nombre de simulation,
                    n, # Taille d'échantillon
                    p, # nombre de metagene
                    size,
                    mu,
                    clayParam){ # parametre de la loi de Poisson
  beta<-rnorm(p+2)
  lambda1<-matrix(0,ncol=1,nrow = B)
  for( j in 1 : B) {
    
    ## Créer une copule de Clayton pour deux variables
    claytonCop <- claytonCopula(param = clayParam, dim = p)
    
    ## Simulation de données à partir de la copule
    sim_data <- rCopula(n, claytonCop)
    
    ## Transformation des données en données discrètes (par exemple, en Poisson)
    X <- apply(sim_data, 2, function(x) qnbinom(x, size = size, mu = mu))
    X<-data.matrix(cbind(rep(1,n), X))
    
    ## X**BETA
    Xbeta<-X%*% matrix(beta[1:(p+1)],ncol=1,nrow=p+1)
    
    ## Condition
    condition <- factor(sample(c("A", "B"), n, replace = TRUE))
    condition <- as.numeric(condition) - 1  # "A" sera 0, "B" sera 1
    
    ## Simulation des erreurs et reconstitution de Y
    
    epsilon<-rnorm(n)
    Y<-Xbeta+beta[p+2]*condition+epsilon
    #return("ca marche")
    ajust<-Ajustement_function(y=Y,X=X,condition=condition)

    chisq1 <- qchisq(1-ajust[,1],1) ; 
    lambda1[j,1] = median(chisq1)/qchisq(0.5,1)
    
    ## Calcul des FP et TN
    
    
  }
  return(lambda1) 
  
}


### Distribution hurdle













































## Descriptive analysis


# Installer les bibliothèques si nécessaire
install.packages("ggplot2")
library(ggplot2)

# Fonction pour simuler une distribution Gamma-Multinomiale
simulate_gamma_multinomial <- function(n, alpha) {
  k <- length(alpha)  # Nombre de catégories
  gamma_sample <- rgamma(k, shape = alpha, rate = 1)  # Échantillons de la distribution Gamma
  p <- gamma_sample / sum(gamma_sample)  # Convertir en probabilités pour la Multinomiale
  rmultinom(1, n, p)  # Multinomiale avec la somme des comptes égale à n
}

# Paramètres à varier
n <- 100  # Taille d'échantillon
alpha_list <- list(c(2, 2, 2), c(5, 1, 1), c(1, 5, 1), c(1, 1, 5))  # Différentes configurations de alpha

# Simuler des données pour chaque configuration de alpha
sim_data <- data.frame()
for (i in 1:length(alpha_list)) {
  alpha <- alpha_list[[i]]
  sim <- simulate_gamma_multinomial(n, alpha)
  sim_data <- rbind(sim_data, data.frame(sim = sim, alpha_group = paste0("alpha = ", paste(alpha, collapse = ","))))
}

# Convertir en format long pour ggplot2
sim_data_long <- data.frame(category = rep(1:nrow(sim_data), each = 3), 
                            count = as.vector(t(sim_data$sim)),
                            alpha_group = rep(sim_data$alpha_group, times = 3))

# Tracer les densités pour chaque configuration de paramètres alpha

grphik_GM<-ggplot(sim_data_long, aes(x = count, fill = alpha_group)) +
  geom_density(alpha = 0.5) +
  labs(title = "", x = "", y = "Density") +
  theme_pubr()



## Hurdle poisson
# Installer les bibliothèques si nécessaire
install.packages("ggplot2")
library(ggplot2)

# Fonction pour simuler une distribution Hurdle Poisson
simulate_hurdle_poisson <- function(n, p_zero, lambda) {
  # Générer des zéros avec une probabilité p_zero
  is_zero <- rbinom(n, size = 1, prob = p_zero) 
  
  # Générer des valeurs non nulles selon une Poisson tronquée (sans zéro)
  non_zero <- rpois(n, lambda)
  non_zero[non_zero == 0] <- 1  # Remplacer les zéros de la Poisson par 1
  
  # Combiner les zéros et les valeurs non nulles
  return(ifelse(is_zero == 1, 0, non_zero))
}

# Paramètres à varier
n <- 1000  # Taille d'échantillon
params <- list(
  list(p_zero = 0.3, lambda = 2),
  list(p_zero = 0.5, lambda = 3),
  list(p_zero = 0.7, lambda = 1)
)

# Simuler des données pour chaque configuration de paramètres
sim_data <- data.frame()
for (i in 1:length(params)) {
  p_zero <- params[[i]]$p_zero
  lambda <- params[[i]]$lambda
  sim <- simulate_hurdle_poisson(n, p_zero, lambda)
  sim_data <- rbind(sim_data, data.frame(sim = sim, param_group = paste0("p_zero = ", p_zero, ", lambda = ", lambda)))
}

# Convertir les données pour ggplot2
sim_data_long <- data.frame(sim_data)

# Tracer les densités pour chaque configuration de paramètres p_zero et lambda
grphik_Hurdle_Poisson<-ggplot(sim_data_long, aes(x = sim, fill = param_group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Densités de la distribution Hurdle Poisson", x = "Valeurs simulées", y = "Densité") +
  theme_pubr()
