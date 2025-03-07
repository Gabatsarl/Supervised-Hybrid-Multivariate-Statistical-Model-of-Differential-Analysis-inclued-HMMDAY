## load the packages

library(ape) ; library(vegan)
library(limma)
library("statmod")
library(edgeR)
library(ggplot2)
library(Rcpp)
library(devtools)
library(metaSPARSim)
library(cowplot)
library(plyr)
library(ggpubr)
library(PRROC)

## load the data
gen_data <- load("/Users/timakpo/Downloads/Data1_Paper2.RData")

groups<-data$groups
GENOM<-data$GENOM

num_zeros <- apply(GENOM == 0, 2, sum)

z=(num_zeros)*100/ncol(GENOM)    # proportion de chaque cas

rm(data)
ind<-which(z==0)
gen_data<-GENOM[,ind]
gen_data<-t(gen_data)




## Distribution binomiale negative

# Charger les donnees deja disponibel
data <- as.matrix(t(gen_data))  # Assurez-vous que vos données sont sous forme de matrice

# Estimation des paramètres pour la distribution de Poisson
lambda_estimates <- colMeans(data)  # Moyenne pour chaque taxon

# Estimation des paramètres pour la distribution binomiale négative
mu_estimates <- colMeans(data)  # Moyenne pour chaque taxon
var_estimates <- apply(data, 2, var)  # Variance pour chaque taxon

# Calculer le paramètre de dispersion pour la binomiale négative
size_estimates <- mu_estimates^2 / (var_estimates - mu_estimates)
size_estimates[size_estimates < 0] <- NA  # Éviter les valeurs négatives


generate_nb_counts <- function(n_samples, n_taxa, size, mu) {
  # Générer les comptages pour chaque échantillon et chaque taxon
  counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
  for (i in 1:n_taxa) {
    counts[, i] <- rnbinom(n_samples, size = size, mu = mu)
  }
  
  return(counts)
}


# Créer des données simulées
B<-5
LAMBDA<-AUC<-matrix(0,nrow=B,ncol = 7)
n_taxa<-nrow(gen_data)
PVAL<-matrix(0,nrow=n_taxa,ncol = B)
colnames(LAMBDA)<-colnames(AUC)<-c("Poisson","NB","Hurdle","MHG","ZIP","DM","BM")
for (k in 1:B){
 
## With a Poisson distribution
generate_poisson_counts <- function(n_samples, n_taxa, lambda_mean, lambda_sd) {
  
  # Générer des taux d'événements (lambda) pour chaque taxon
  lambda <- rnorm(n_taxa, mean = lambda_mean, sd = lambda_sd)
  lambda[lambda < 0] <- 0  # Assurer que lambda est non négatif
  
  # Générer les comptages pour chaque échantillon et chaque taxon
  counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
  for (i in 1:n_taxa) {
    counts[, i] <- rpois(n_samples, lambda = lambda[i])
  }
  
  return(counts)
}
  
  ## Create a simulate data
  lambda_estimates <- rowMeans(gen_data)  # Moyenne pour chaque taxon

    sim_data_poisson <- generate_poisson_counts(nrow(gen_data), ncol(gen_data),
                                              lambda_estimates,
                                              sqrt(lambda_estimates))
  m.D <- vegdist(t(sim_data_poisson), "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  ## ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  ## Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =sim_data_poisson)
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_pois<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,1]<-median(qchisq(1-pval_pois,1))/qchisq(0.5,1)
  PVAL[,k]<-pval_pois
  
  ## Calcul for AUC for this courbe
  pval_pois<-pval_pois
  #labels <- ifelse(PVAL[,4] > 0.05, 0, 1)
  
  proc<-roc.curve(scores.class0 = pval_pois, scores.class1 = pval_pois,
                  weights.class0 = pval_pois, weights.class1 = 1-pval_pois, curve = TRUE)
  AUC[k,1]<-proc$auc
  
  ## Distribution binomiale negative
  
  # Charger les donnees deja disponibel
  data <- as.matrix(t(gen_data))  # Assurez-vous que vos données sont sous forme de matrice
  

  
  # Estimation des paramètres pour la distribution binomiale négative
  mu_estimates <- colMeans(data)  # Moyenne pour chaque taxon
  var_estimates <- apply(data, 2, var)  # Variance pour chaque taxon
  
  # Calculer le paramètre de dispersion pour la binomiale négative
  size_estimates <- mu_estimates^2 / (var_estimates - mu_estimates)
  size_estimates[size_estimates < 0] <- NA  # Éviter les valeurs négatives
  
  
  generate_nb_counts <- function(n_samples, n_taxa, size, mu) {
    # Générer les comptages pour chaque échantillon et chaque taxon
    counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
    for (i in 1:n_taxa) {
      counts[, i] <- rnbinom(n_samples, size = size, mu = mu)
    }
    
    return(counts)
  }
  
  sim_data_nb <- generate_nb_counts(nrow(gen_data), ncol(gen_data), 
                                    size_estimates, mu_estimates)
  
  
  m.D <- vegdist(t(sim_data_nb), "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  # ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =sim_data_nb)
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_NB<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,2]<-median(qchisq(1-pval_NB,1))/qchisq(0.5,1)
  PVAL[,2]<-mean(c(PVAL[,2],pval_NB))
  
  
  
  
  # DISTRIBUTION 3 :::: ----- Distribution de Hurdle
  
  # Générer les données de présence/absence
  prob_zero <- 0.3  # Probabilité que le comptage soit zéro
  n_samples<-ncol(gen_data) ;   n_taxa<-nrow(gen_data)
  presence_absence <- rbinom(n_samples * n_taxa, 1, 1 - prob_zero)
  presence_absence_matrix <- matrix(presence_absence, nrow = n_samples,
                                    ncol = n_taxa)
  
  # Générer les comptages pour les valeurs non-nulles
  lambda_mean<-lambda_estimates
  sim_data_Hurdle <- matrix(0, nrow = n_samples, ncol = n_taxa)
  for (j in 1:n_taxa) {
    sim_data_Hurdle[, j] <- rpois(n_samples, lambda = lambda_mean) * 
      presence_absence_matrix[, j]
  }
  print(dim(sim_data_Hurdle))
  m.D <- vegdist(sim_data_Hurdle, "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  # ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  # Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =t(sim_data_Hurdle))
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_Hurdle<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,3]<-median(qchisq(1-pval_Hurdle,1))/qchisq(0.5,1)
  
  
  
  ## Loi geometric Multivariate hypergeometric (MHG) : 
  
  depths <- rowSums(gen_data)
  
  ## Normaliser les données en divisant chaque comptage par la profondeur de séquençage de l'échantillon
  
  gen_data_norm <- t(apply(gen_data, 1, function(x) x / sum(x) * median(depths)))
  
  indpop2<-list(ZCarex_OSPW=which(groups=="ZCarex_OSPW"),
               Non_Planted_OSPW=which(groups=="Non_Planted_OSPW"))
  
  
  params2<-estimate_parameter_from_data(gen_data, gen_data_norm, 
                                       indpop2,perc_not_zeros=.2)
  
  
  
  names(params2)<-names(indpop2)
  sim_data_GMH<-metaSPARSim(params2)
  
  m.D <- vegdist(t(sim_data_GMH$counts), "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  ## ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  ## Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
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
  
  ## ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  ## Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =sim_data_ZIP)
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_ZIP<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,5]<-median(qchisq(1-pval_ZIP,1))/qchisq(0.5,1)

  
  
  ## Direchelot multinomial
  # Paramètres alpha pour la distribution Dirichlet
  alpha <- rep(1, n_taxa)  # On peut ajuster ces valeurs selon le cas
  
  # Générer les proportions Dirichlet
  library(MCMCpack)
  dirichlet_proportions <- rdirichlet(n_samples, alpha)
  
  # Simuler les comptages Multinomiaux à partir des proportions Dirichlet
  sim_data_DM <- matrix(0, nrow = n_samples, ncol = n_taxa)
  
  for (i in 1:n_samples) {
    sim_data_DM[i, ] <- rmultinom(1, size = sum(gen_data[, i]),
                                  prob = dirichlet_proportions[i, ])
  }
  m.D <- vegdist(sim_data_DM, "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  # ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  # Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =t(sim_data_DM))
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_DM<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,6]<-median(qchisq(1-pval_DM,1))/qchisq(0.5,1)
  
## BETA multinomial : http://prob140.org/textbook/content/Chapter_21/02_Beta_Binomial_Distribution.html

  # Paramètres alpha et beta pour la distribution bêta-binomiale
  alpha <- runif(n_taxa, 0, 1)  # Paramètres alpha
  beta <- runif(n_taxa, 0, 1)    # Paramètres beta
  
  ## Matrice pour stocker les données simulées
  sim_data_BB <- matrix(0, nrow = n_samples, ncol = n_taxa)
  
  ## Simuler les données bêta-binomiales
  #size<-colSums(gen_data)
  for (i in 1:n_samples) {
    ## Générer des proportions bêta-binomiales pour chaque taxon
     proportions <- rbeta(n_taxa, alpha[i], beta[i])
    
    ## Normaliser les proportions pour qu'elles somment à la taille totale
     size_i<-sum(gen_data[, i]) 
     proportions <- proportions / sum(proportions) * size_i
    
    ## Générer des comptages binomiaux
     sim_data_BB[i, ] <-rmultinom(1, size = sum(gen_data[, i]),
                                  prob = proportions) #rbinom(n_taxa, size = sum(gen_data[, i]), prob = proportions)
  }
  
  m.D <- vegdist(sim_data_BB, "manhattan") ; 
  result_GENERAL <- pcoa(m.D )
  
  # ================================================##
  
  data<-data.frame(groups=groups)
  data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
  data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
  
  # Design 
  
  design <- model.matrix(~groups+PCoA1+PCoA2,data)
  d <- DGEList(counts =t(sim_data_BB))
  y <- voom(d, design, plot = FALSE)
  limma_fit <- lmFit(y,design=design)
  tmp <- contrasts.fit(limma_fit, coef = 2) 
  out_ebayes <- eBayes(tmp)
  
  # result<-topTable(out_ebayes,1)
  result_Limma <- topTable(out_ebayes, number = Inf)
  pval_BB<-result_Limma$P.Value # $adj.P.Val
  
  LAMBDA[k,7]<-median(qchisq(1-pval_BB,1))/qchisq(0.5,1)
  
}



## Transformer les données au format long avec la fonction melt

library(reshape2)
data_long <- melt(LAMBDA)

# Tracer le boxplot avec ggplot2

ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot() +
  labs(title = "", x = "", y = "Lambda") +
  theme_pubr()+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",lwd=1)+
  scale_fill_discrete(guide = "none") 


ggsave("/Users/timakpo/Downloads/simul.png",plot=gplot,
       width = 7, height = 8, units = "in", dpi = 300,bg = "white")




## Recal precision avec un seuil=0.05


library(PRROC)
m.D <- vegdist(t(gen_data), "manhattan") ; 
result_GENERAL <- pcoa(m.D )

## ================================================##

data<-data.frame(groups=groups)
data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]

## Design 

design <- model.matrix(~groups+PCoA1+PCoA2,data)
d <- DGEList(counts =gen_data)
y <- voom(d, design, plot = FALSE)
limma_fit <- lmFit(y,design=design)
tmp <- contrasts.fit(limma_fit, coef = 2) 
out_ebayes <- eBayes(tmp)

# result<-topTable(out_ebayes,1)
result_Limma <- topTable(out_ebayes, number = Inf)
pval<-result_Limma$P.Value # $adj.P.Val


scores<-PVAL[,4]
labels <- ifelse(PVAL[,4] > 0.05, 0, 1)

proc<-roc.curve(scores.class0 = scores, scores.class1 = scores,
                 weights.class0 = scores, weights.class1 = 1-scores, curve = TRUE)
plot(proc)

pr <- pr.curve(scores.class0 = scores[labels == 0],
               scores.class1 = scores[labels == 1],
               curve = TRUE)

plot(pr,rand.plot = TRUE, fill.area = TRUE,
      fill.color = rgb(0.8,1,0.8), maxminrand.col = "blue")
