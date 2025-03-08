## load the packages

library(ape) ; library(vegan) ; library(limma) ; library("statmod") ; library(edgeR) ; library(ggplot2) ; library(Rcpp) ; 
library(devtools) ; library(metaSPARSim) ; library(cowplot) ; library(plyr) ; library(ggpubr) ; library(PRROC)

## Create the function of load the data

simulation_NB=function(B, # nombre de simulation,
                       n_samples, # Taille d'échantillon
                       n_taxa, # nombre de metagene
                       lambda_estimates, # parametres de la poisson
                       size, # longeur de l'abondance (dim=n_sample multinomiale et size)
                       mu_nb, # moyenne de la Poisson qui est celui du Gamma
                       sigma, # probabilité de zero dans les données (dimemsion)
                       groups){
  
  ## parametre de la loi de Poisson
  LAMBDA<-AUC<-matrix(0,nrow=B,ncol = 7)
  #n_taxa<-n_taxa
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
    #lambda_estimates <- rowMeans(gen_data)  # Moyenne pour chaque taxon
    
    sim_data_poisson <- generate_poisson_counts(n_taxa, n_samples,
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

    generate_nb_counts <- function(n_samples, n_taxa, size, mu) {
      # Générer les comptages pour chaque échantillon et chaque taxon
      counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
      for (i in 1:n_taxa) {
        counts[, i] <- rnbinom(n_samples, size = size, mu = mu)
      }
      
      return(counts)
    }
    size_nb<-runif(n_taxa,0,1)
    mu_nb<-sample(n_taxa,n_taxa)
    sim_data_nb <- generate_nb_counts(n_taxa, n_samples, 
                                      size=size_nb, mu=mu_nb)
    
    
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
    n_samples<-n_samples ;   n_taxa<-n_taxa
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
    
    #m.D <- vegdist(t(sim_data_MHG), "manhattan") ; 
    #result_GENERAL <- pcoa(m.D )
    
    ## ================================================##
    
    #data<-data.frame(groups=groups)
    #data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    #data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design 
    
    #design <- model.matrix(~groups+PCoA1+PCoA2,data)
    #d <- DGEList(counts =sim_data_GMH$counts)
    #y <- voom(d, design, plot = FALSE)
    #limma_fit <- lmFit(y,design=design)
    #tmp <- contrasts.fit(limma_fit, coef = 2) 
    #out_ebayes <- eBayes(tmp)
    
    ## result<-topTable(out_ebayes,1)
    #result_Limma <- topTable(out_ebayes, number = Inf)
    #pval_GMH<-result_Limma$P.Value # $adj.P.Val
    
    #LAMBDA[k,4]<-median(qchisq(1-pval_GMH,1))/qchisq(0.5,1)
    
    
    
    ## Zero-Inflated Poisson (ZIP)
    #pi_inflate <- 0.3
    #n <-n_samples*n_taxa
    # Créer la matrice de microbiome simulée
    sim_data_ZIP <- matrix(0, nrow = n_taxa, ncol = n_samples)
    
    ## Remplir la matrice avec les données ZIP
    
    library(gamlss.dist)
    for (i in 1:n_taxa) {
      sim_data_ZIP[i, ] <- rZIP(n_samples, mu =mu_nb[i]  , 
                                sigma = 0.3) # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
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
    alpha <- runif(n_taxa, 0, 1)  # Paramètres alpha
    
    # Générer les proportions Dirichlet
    library(MCMCpack)
    dirichlet_proportions <- rdirichlet(n_samples, alpha)
    
    # Simuler les comptages Multinomiaux à partir des proportions Dirichlet
    sim_data_DM <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    for (i in 1:n_samples) {
      sim_data_DM[i, ] <- rmultinom(1, size = size,
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
    beta <- runif(n_taxa, 0, 1)    # Paramètres beta
    
    ## Matrice pour stocker les données simulées
    sim_data_BB <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    ## Simuler les données bêta-binomiales
    #size<-colSums(gen_data)
    for (i in 1:n_samples) {
      ## Générer des proportions bêta-binomiales pour chaque taxon
      proportions <- rbeta(n_taxa, alpha[i], beta[i])
      
      ## Générer des comptages binomiaux
      sim_data_BB[i, ] <-rmultinom(1, size =size,
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
  return(list(LAMBDA=LAMBDA, AUC=AUC))
}

n_samples<-30
n_taxa<-50
lambda_estimates=runif(n_taxa,1,10)
size<-sample(n_sample:(2*n_sample),n_sample)
mu_nb<-sample(n_taxa,n_taxa)
sigma<-0.4
groups=sample(c("Non_Planted_OSPW","ZCarex_OSPW"),size=n_samples,replace = TRUE)

RESULT<-simulation_NB(B=10, # nombre de simulation,
                               n_samples=n_samples, # Number of sample
                               n_taxa=n_taxa, # Number of taxa
                               lambda_estimates=lambda_estimates, # parametres de la poisson
                               size=size, # longeur de l'abondance (dim=n_sample multinomiale et size)
                               mu_nb=mu_nb, # moyenne de la Poisson qui est celui du Gamma
                               sigma=sigma, # probabilité de zero dans les données (dimemsion)
                               groups=groups)
  

library(reshape2)
LAMBDA<-RESULT[[1]]
data_long <- melt(LAMBDA)

# Tracer le boxplot avec ggplot2

gplot<-ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot() +
  labs(title = "", x = "", y = "Lambda") +
  theme_pubr()+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red",lwd=1)+
  scale_fill_discrete(guide = "none") 


ggsave("/Users/timakpo/Downloads/simul_nodata.png",plot=gplot,
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
