## Ce fichier contient le code de génération des données X qui est une copie des données de base X_obs suivant plusieurs distribution
## Input :
##
## Output : 


## Load the packages

pacman::p_load(ape, vegan, limma, statmod, edgeR, ggplot2, Rcpp, devtools, metaSPARSim, cowplot, plyr, ggpubr, PRROC, gamlss.dist, 
               doParallel, reshape2, MCMCpack, variancePartition, BiocParallel, tidyr, dplyr)

## Créer des données simulées

Simul_with_data<-function(B, #Nombre de replication
                          data, # Matrix de "n ligne et 2 colonnes" contenant la "Condition=group" et "Y"
                          gen_data # contient les abondances 

                          ){
  n_taxa<-nrow(gen_data)
  PVAL<-matrix(0,nrow=n_taxa,ncol = B)
  LAMBDA<-matrix(0,nrow=B,ncol = 6)
  colnames(LAMBDA)<-c("Poisson","NB","Hurdle","ZIP","DM","BM")
  
  ESTIMATION1<-ESTIMATION2<-ESTIMATION3<-ESTIMATION4<-ESTIMATION5<-ESTIMATION6<-
    PVALUE1<-PVALUE2<-PVALUE3<-PVALUE4<-PVALUE5<-PVALUE6<-
    LOGFC1<-LOGFC2<-LOGFC3<-LOGFC4<-LOGFC5<-LOGFC6<-matrix(0,ncol=B,nrow = n_taxa)
  
  
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
    
    X_Poisson <- generate_poisson_counts(nrow(gen_data), ncol(gen_data),
                                                lambda_estimates,
                                                sqrt(lambda_estimates))

    ## PCoA for X
    m.D <- vegdist(t(X_Poisson), "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    ## ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design of the model
    
    for(i in 1:n_taxa){
      
      data$Xg<-log10(X_Poisson[i,]+1)
      design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
      
      ## Fitting the model 
      
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      ## pcoa
      para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
      ESTIMATION1[i,k] <- limma_fit$coefficients[2] ; 
      PVALUE1[i,k] <- pval[2] ; 
      LOGFC1[i,k] <- out_ebayes$p.value[2]
      
      
    }
    
    LAMBDA[k,1]<-median(qchisq(1-PVALUE1[,k],1),na.rm = "TRUE")/qchisq(0.5,1)

    
    ## Charger les donnees deja disponibel
    gen_data2 <- as.matrix(t(gen_data))  # Assurez-vous que vos données sont sous forme de matrice
    
    
    ## Estimation des paramètres pour la distribution binomiale négative
    mu_estimates <- colMeans(gen_data2)  # Moyenne pour chaque taxon
    var_estimates <- apply(gen_data2, 2, var)  # Variance pour chaque taxon
    rm(gen_data2)
    
    ## Calculer le paramètre de dispersion pour la binomiale négative
    size_estimates <- mu_estimates^2 / (var_estimates - mu_estimates)
    size_estimates[size_estimates < 0] <- 0.001  # Éviter les valeurs négatives
    
    
    generate_nb_counts <- function(n_samples, n_taxa, size, mu) {
      ## Générer les comptages pour chaque échantillon et chaque taxon
      counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
      for (i in 1:n_taxa) {
        counts[, i] <- rnbinom(n_samples, size = size, mu = mu)
      }
      
      return(counts)
    }
    
    X_NB <- generate_nb_counts(nrow(gen_data), ncol(gen_data),  size_estimates, mu_estimates)
    
    
    m.D <- vegdist(t(X_NB), "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]

    ## Design of the model
    
    for(i in 1:n_taxa){
      data$Xg <- log10(X_NB[i,] + 1)
      design <- model.matrix(~Xg + groups + PCoA1 + PCoA2, data)
      
      # Assurez-vous que les dimensions correspondent avant de procéder
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      para <- limma_fit$coefficients
      pval <- out_ebayes$p.value
      
      ESTIMATION2[i, k] <- limma_fit$coefficients[2]
      PVALUE2[i, k] <- pval[2]
      LOGFC2[i, k] <- out_ebayes$p.value[2]
      
    }
    
    LAMBDA[k,2]<-median(qchisq(1-PVALUE2[,k],1),na.rm = "TRUE")/qchisq(0.5,1)

    ## DISTRIBUTION 3 ::::------ Distribution de Hurdle
    
    ## Générer les données de présence/absence
    prob_zero <- 0.3  # Probabilité que le comptage soit zéro
    n_samples<-ncol(gen_data) ;   n_taxa<-nrow(gen_data)
    presence_absence <- rbinom(n_samples * n_taxa, 1, 1 - prob_zero)
    presence_absence_matrix <- matrix(presence_absence, nrow = n_samples,
                                      ncol = n_taxa)
    
    ## Générer les comptages pour les valeurs non-nulles
    lambda_mean<-lambda_estimates
    X_Hurdle <- matrix(0, nrow = n_samples, ncol = n_taxa)
    for (j in 1:n_taxa) {
      X_Hurdle[, j] <- rpois(n_samples, lambda = lambda_mean) * 
        presence_absence_matrix[, j]
    }
    m.D <- vegdist(X_Hurdle, "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    ## ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design of the model
    for(i in 1:n_taxa){
      
      data$Xg<-log10(X_Hurdle[,i]+1)
      design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
      
      ## Fitting the model 
      
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      ## pcoa
      para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
      
      ESTIMATION3[i,k] <- limma_fit$coefficients[2] ; 
      PVALUE3[i,k] <- pval[2] ; 
      LOGFC3[i,k] <- out_ebayes$p.value[2]
      
      
    }
    
    LAMBDA[k,3]<-median(qchisq(1-PVALUE3[,k],1),na.rm = TRUE)/qchisq(0.5,1)

    
    ## Model : Zero-Inflated Poisson (ZIP)
    
    ## Génération de données ZIP
    lambda_est <- mean(gen_data[gen_data > 0])  # Moyenne des comptages non-nuls

    # Créer la matrice de microbiome simulée
    X_ZIP <- matrix(0, nrow = n_taxa, ncol = n_samples)

    ## Remplir la matrice avec les données ZIP
    
    for (i in 1:n_taxa) {
      sigma<-sd(gen_data[i,],na.rm = TRUE)
      if(sigma==0){
      X_ZIP[i, ] <- rZIP(n_samples, mu =mean(gen_data[i,])  , 
                                sigma =)} # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
      else {
        X_ZIP[i, ] <- rZIP(n_samples, mu =mean(gen_data[i,])  , 
                                  sigma =0.001)} # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
    
      }
    
    
    ## PCoA from matrix X
    
    m.D <- vegdist(t(X_ZIP), "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    ## ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design of the model
    for(i in 1:n_taxa){
      
      data$Xg<-log10(X_ZIP[i,]+1)
      design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
      
      ## Fitting the model 
      
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      ## pcoa
      para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
      
      ESTIMATION4[i,k] <- limma_fit$coefficients[2] ; 
      PVALUE4[i,k] <- pval[2] ; 
      LOGFC4[i,k] <- out_ebayes$p.value[2]
      
      
    }
    LAMBDA[k,4]<-median(qchisq(1-PVALUE4[,k],1),na.rm = TRUE)/qchisq(0.5,1)

    
    ## Direchelot multinomial
    alpha <- rep(1, n_taxa)  # On peut ajuster ces valeurs selon le cas
    
    # Générer les proportions Dirichlet
    dirichlet_proportions <- rdirichlet(n_samples, alpha)
    
    # Simuler les comptages Multinomiaux à partir des proportions Dirichlet
    X_DM <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    for (i in 1:n_samples) {
      X_DM[i, ] <- rmultinom(1, size = sum(gen_data[, i]),
                                    prob = dirichlet_proportions[i, ])
    }

    ## PCoA from matrix X
    m.D <- vegdist(X_DM, "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design of the model
    for(i in 1:n_taxa){
      
      data$Xg<-log10(X_DM[,i]+1)
      design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
      
      ## Fitting the model 
      
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      ## Result
      
      para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
      
      ESTIMATION5[i,k] <- limma_fit$coefficients[2] ; 
      PVALUE5[i,k] <- pval[2] ; 
      LOGFC5[i,k] <- out_ebayes$p.value[2]
      
      
    }
    
    pval_DM<-PVALUE5[,k]
    LAMBDA[k,5]<-median(qchisq(1-pval_DM,1),na.rm = "TRUE")/qchisq(0.5,1)

    
    ## BETA multinomial : http://prob140.org/textbook/content/Chapter_21/02_Beta_Binomial_Distribution.html
    
    # Paramètres alpha et beta pour la distribution bêta-binomiale
    alpha <- runif(n_taxa, 0, 1)  # Paramètres alpha
    beta <- runif(n_taxa, 0, 1)    # Paramètres beta
    
    ## Matrice pour stocker les données simulées
    X_BB <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    ## Simuler les données bêta-binomiales
    #size<-colSums(gen_data)
    for (i in 1:n_samples) {
      ## Générer des proportions bêta-binomiales pour chaque taxon
      proportions <- rbeta(n_taxa, alpha[i], beta[i])
      
      ## Normaliser les proportions pour qu'elles somment à la taille totale
      size_i<-sum(gen_data[, i]) 
      proportions <- proportions / sum(proportions) * size_i
      
      ## Générer des comptages binomiaux
      X_BB[i, ] <-rmultinom(1, size = sum(gen_data[, i]),
                                   prob = proportions) #rbinom(n_taxa, size = sum(gen_data[, i]), prob = proportions)
    }
    ## PCoA from matrix X

    m.D <- vegdist(X_BB, "manhattan",na.rm = TRUE) ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ## Design of the model
    for(i in 1:n_taxa){
      
      data$Xg<-log10(X_BB[,i]+1)
      design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
      
      ## Fitting the model 
      
      limma_fit  <- lmFit(data$Y, design)
      out_ebayes <- eBayes(limma_fit)
      
      ## pcoa
      para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
      ESTIMATION6[i,k] <- limma_fit$coefficients[2] ; 
      PVALUE6[i,k] <- pval[2] ; 
      LOGFC6[i,k] <- out_ebayes$p.value[2]
      
      
    }

    LAMBDA[k,6]<-median(qchisq(1-PVALUE6[,k],1),na.rm = "TRUE")/qchisq(0.5,1)    
    print(k)
    
  }
  
  return(list(LAMBDA=LAMBDA,
              ESTIMATION1=ESTIMATION1,ESTIMATION2=ESTIMATION2,ESTIMATION3=ESTIMATION3,
              ESTIMATION4=ESTIMATION4, ESTIMATION5=ESTIMATION5,ESTIMATION6=ESTIMATION6,
              PVALUE1=PVALUE1,PVALUE2=PVALUE2, PVALUE3=PVALUE3, PVALUE4=PVALUE4, PVALUE5=PVALUE5, PVALUE6=PVALUE6,
              LOGFC1=LOGFC1,LOGFC2=LOGFC2,LOGFC3=LOGFC3,LOGFC4=LOGFC4,LOGFC5=LOGFC5,LOGFC6=LOGFC6))
  
}

## Qqplot

myqqplot2 = function(pval,seuil){
  obs = -log10(sort(pval, decreasing = FALSE))
  expect = -log10(1:length(obs) / length(obs))
  df = data.frame(Observed = obs, Expected = expect, pval = pval)
  
  # Création du graphique avec points rouges pour y > 3
  plot1 = ggplot(df, aes(x = Expected, y = Observed, color = Observed > seuil)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("black", "red")) +  # Couleurs pour les points
    scale_x_continuous(name = expression(-log[10](italic(p)) ~~ Expected)) +
    scale_y_continuous(name = expression(-log[10](italic(p)) ~~ Observed), breaks = 1:6) +
    theme_pubr()+
    theme( axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          legend.position = "none"  # Ne pas afficher la légende de couleur
         )+
  geom_abline(intercept = 0, slope = 1, color="red",lwd=1,size=1.5)
  
  return(plot1)

  
}

## Volcanoplot

volcanoplot <- function(x, y, tit,col) {
  df <- data.frame(x = x, y = y)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = ifelse(y > 5, "Above 5", "Below 5")), size = 2) +
    scale_color_manual(name = "Values", values = c("Above 5" = col, "Below 5" = "#808080")) +
    scale_x_continuous(name = expression(Log[2] ~ Fold~Change)) +
    scale_y_continuous(name = expression(-Log[10](italic(p))~Observed),
                       breaks = c(1:9)) +
    labs(title = tit) +
    theme_classic() +
    theme(text = element_text(family = "Times New Roman"),
          axis.text.x = element_text(size = 12, angle = 360, hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          legend.title =element_blank(),
          legend.position ="none")
}



