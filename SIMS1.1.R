-- Load the packages

library(ape) ; library(vegan) ; library(limma) ; library("statmod")
library(edgeR) ; library(ggplot2) ; library(Rcpp) ; library(devtools)
library(metaSPARSim) ; library(cowplot) ; library(plyr) ; library(ggpubr)
library(PRROC) ; library(gamlss.dist) ; library(doParallel) ; library(reshape2) ; library(MCMCpack)


##-- Create the function of load the data -----------##

simulation_NB<-function(B, # Number of replication,
                        beta, alpha, # parameters of modelling 
                        n_samples, # Sample size
                        n_taxa, # Number of metagens
                        lambda_estimates, # Parameters of distribution
                        size, # Size of abundance (dim=n_sample multinomiale et size)
                        mu_nb, # Mean of Poisson et Gamma distributions
                        sigma, # Probability of zero
                        groups){
  
  ## Parameters of Poisson distribution
  LAMBDA<-AUC<-matrix(0,nrow=B,ncol = 6)
  n0<-length(which(beta!=0)) # Number of significatif genes
  
  #n_taxa<-n_taxa
  #PVAL<-matrix(0,nrow=n_taxa,ncol = B)
  
  ESTIMATION1<-ESTIMATION2<-ESTIMATION3<-ESTIMATION4<-ESTIMATION5<-ESTIMATION6<-
  PVALUE1<-PVALUE2<-PVALUE3<-PVALUE4<-PVALUE5<-PVALUE6<-
  LOGFC1<-LOGFC2<-LOGFC3<-LOGFC4<-LOGFC5<-LOGFC6<-matrix(0,ncol=B,nrow = n_taxa)
  #colnames(ESTIMATION)<-c("beta.1") #
  #colnames(PVALUE)<-c("pval.1") #pval.i :=pvalue de beta.i
  #colnames(LOGFC)<-c("logFCbeta.1")
  
  colnames(LAMBDA)<-colnames(AUC)<-c("Poisson","NB","Hurdle","ZIP","DM","BM")
  for (k in 1:B){
    
    ## With a Poisson distribution
    generate_poisson_counts <- function(n_samples, n_taxa, lambda_mean, lambda_sd) {
      
      ## Generate of number of situation for taxon
      lambda <- rnorm(n_taxa, mean = lambda_mean, sd = lambda_sd)
      lambda[lambda < 0] <- 0 # Negatif distribution
      
      ## Generate les comptages pour chaque échantillon et chaque taxon
      counts <- matrix(0, nrow = n_samples, ncol = n_taxa)
      for (i in 1:n_taxa) {
        counts[, i] <- rpois(n_samples, lambda = lambda[i])
      }
      
      return(counts)
    }
    
    ## Create a simulate data
    #lambda_estimates <- rowMeans(gen_data) # Moyenne pour chaque taxon
    
    sim_data_poisson <- generate_poisson_counts(n_taxa, n_samples,
                                                lambda_estimates,
                                                sqrt(lambda_estimates)) # X matrix
    
    m.D <- vegdist(t(sim_data_poisson), "manhattan") ; 
    result_GENERAL <- pcoa(m.D)
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-t(sim_data_poisson)%*%beta
    eps<-rt(n_samples,df=5)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_poisson[i,]+1)
        print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        print(colnames(design))
        ## Fitting the model 
        
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        ## pcoa
        #d <- DGEList(counts =t(sim_data_Hurdle))
        #y <- voom(d, design, plot = FALSE)
        #limma_fit <- lmFit(y,design=design)
        #tmp <- contrasts.fit(limma_fit, coef = 2) 
        #out_ebayes <- eBayes(tmp)
        
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION1[i, k] <- NA ; PVALUE1[i, k] <- NA ; LOGFC1[i, k] <- NA
      }
      else {
        ESTIMATION1[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE1[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC1[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_pois<-PVALUE1[,k]
    LAMBDA[k,1]<-median(qchisq(1-pval_pois,1),na.rm = "TRUE")/qchisq(0.5,1)
    
    ## Calcul for AUC for this courbe
    labels <- c(rep(1,n0),rep(0,n_taxa-n0))
    proc<-roc.curve(scores.class0 = pval_pois, weights.class0 = labels)
                      
    #proc<-roc.curve(scores.class0 = pval_pois[1:(n0)], scores.class1 = 
     #                 pval_pois[(n0+1):n_taxa],curve = TRUE)
                  #  weights.class0 = pval_pois[], weights.class1 = 1-pval_pois, curve = TRUE)
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
    
    
    ##===========Calculate Y===================## 
    
    Xbeta<-t(sim_data_nb)%*%beta
    #eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_nb[i,]+1)
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        
        ## Fitting the model 
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION2[i, k] <- NA ; PVALUE2[i, k] <- NA ; LOGFC2[i, k] <- NA
      }
      else {
        ESTIMATION2[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE2[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC2[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_NB<-PVALUE2[,k]
    LAMBDA[k,2]<-median(qchisq(1-pval_NB,1),na.rm = "TRUE")/qchisq(0.5,1)
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_NB[1:(n0)], scores.class1 = 
                      pval_NB[(n0+1):n_taxa],curve = TRUE)
    #  weights.class0 = pval_pois[], weights.class1 = 1-pval_pois, curve = TRUE)
    AUC[k,2]<-proc$auc
    
    ## DISTRIBUTION 3 :::: ----- Distribution de Hurdle
    
    ## Generate the data presence/absence
    prob_zero <- 0.3 # Probability que le comptage soit zéro
    n_samples<-n_samples ; n_taxa<-n_taxa
    presence_absence <- rbinom(n_samples * n_taxa, 1, 1 - prob_zero)
    presence_absence_matrix <- matrix(presence_absence, nrow = n_samples,
                                      ncol = n_taxa)
    
    ## Generate the count of non-nuls
    lambda_mean<-lambda_estimates
    sim_data_Hurdle <- matrix(0, nrow = n_samples, ncol = n_taxa)
    for (j in 1:n_taxa) {
      sim_data_Hurdle[, j] <- rpois(n_samples, lambda = lambda_mean) * 
        presence_absence_matrix[, j]
    }
    m.D <- vegdist(sim_data_Hurdle, "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    
    ##===========Calculate Y===================## 
    
    Xbeta<-sim_data_Hurdle%*%beta
    #eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_Hurdle[,i]+1)
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        
        ## Fitting the model 
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION3[i, k] <- NA ; PVALUE3[i, k] <- NA ; LOGFC3[i, k] <- NA
      }
      else {
        ESTIMATION3[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE3[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC3[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_Hurdle<-PVALUE3[,k]
    LAMBDA[k,3]<-median(qchisq(1-pval_Hurdle,1))/qchisq(0.5,1)
    
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_Hurdle[1:(n0)], scores.class1 = 
                      pval_Hurdle[(n0+1):n_taxa],curve = TRUE)
    #  weights.class0 = pval_pois[], weights.class1 = 1-pval_pois, curve = TRUE)
    AUC[k,3]<-proc$auc
    
    
    ## Zero-Inflated Poisson (ZIP)
    #pi_inflate <- 0.3
    #n <-n_samples*n_taxa
    # Créer la matrice de microbiome simulée
    sim_data_ZIP <- matrix(0, nrow = n_taxa, ncol = n_samples)
    
    ## Remplir la matrice avec les données ZIP
    
    for (i in 1:n_taxa) {
      sim_data_ZIP[i, ] <- rZIP(n_samples, mu =mu_nb[i] , 
                                sigma = 0.3) # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
    }
    
    m.D <- vegdist(t(sim_data_ZIP), "manhattan") ; 
    result_GENERAL <- pcoa(m.D)
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-t(sim_data_ZIP)%*%beta
    #eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_ZIP[i,]+1)
        print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        print(colnames(design))
        ## Fitting the model 
        
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION4[i, k] <- NA ; PVALUE4[i, k] <- NA ; LOGFC4[i, k] <- NA
      }
      else {
        ESTIMATION4[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE4[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC4[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_ZIP<-PVALUE4[,k]
    LAMBDA[k,4]<-median(qchisq(1-pval_ZIP,1))/qchisq(0.5,1)
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_ZIP[1:(n0)], scores.class1 = 
                      pval_ZIP[(n0+1):n_taxa],curve = TRUE)
    #  weights.class0 = pval_pois[], weights.class1 = 1-pval_pois, curve = TRUE)
    AUC[k,4]<-proc$auc
    
    ## Direchelot multinomial
    # Paramètres alpha_direchlet pour la distribution Dirichlet
    alpha_direchlet <- runif(n_taxa, 0, 1) # Paramètres alpha_ of direchlet
    
    # Générer les proportions Dirichlet
    dirichlet_proportions <- rdirichlet(n_samples, alpha_direchlet)
    
    # Simuler les comptages Multinomiaux à partir des proportions Dirichlet
    sim_data_DM <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    for (i in 1:n_samples) {
      sim_data_DM[i, ] <- rmultinom(1, size = size,
                                    prob = dirichlet_proportions[i, ])
    }
    m.D <- vegdist(sim_data_DM, "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-sim_data_DM%*%beta
    #eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_DM[,i]+1)
        print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        print(colnames(design))
        ## Fitting the model 
        
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION5[i, k] <- NA ; PVALUE5[i, k] <- NA ; LOGFC5[i, k] <- NA
      }
      else {
        ESTIMATION5[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE5[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC5[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_DM<-PVALUE5[,k]
    LAMBDA[k,5]<-median(qchisq(1-pval_DM,1),na.rm = "TRUE")/qchisq(0.5,1)
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_DM[1:(n0)], scores.class1 = 
                      pval_DM[(n0+1):n_taxa],curve = TRUE)
    AUC[k,5]<-proc$auc
    
    ## BETA multinomial : http://prob140.org/textbook/content/Chapter_21/02_Beta_Binomial_Distribution.html
    
    # Paramètres alpha et beta pour la distribution bêta-binomiale
    beta_beta <- runif(n_taxa, 0, 1) # Parameters beta_beta
    
    ## Matrix 
    sim_data_BB <- matrix(0, nrow = n_samples, ncol = n_taxa)
    
    ## Simulation of data 
    #size<-colSums(gen_data)
    for (i in 1:n_samples) {
      ## Generate of proportion beta-binomial
      proportions <- rbeta(n_taxa, alpha_direchlet[i], beta_beta[i])
      
      ## Generate of binomial distribution
      sim_data_BB[i, ] <-rmultinom(1, size =size,
                                   prob = proportions) #rbinom(n_taxa, size = sum(gen_data[, i]), prob = proportions)
    }
    
    m.D <- vegdist(sim_data_BB, "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-sim_data_BB%*%beta
    #eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_BB[,i]+1)
        print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        print(colnames(design))
        ## Fitting the model 
        
        limma_fit  <- lmFit(data$Y, design)
        out_ebayes <- eBayes(limma_fit)
        
        if(inherits(limma_fit, "try-error")){ cat("NULL") }
        else {
          para <- limma_fit$coefficients ; pval <- out_ebayes$p.value ;
          
          print(c(para[2],pval[2],topTable(out_ebayes,coef=2)$logFC))
        }
      }
      
    )
    
    ## Arreter le cluster de calcul
    stopCluster(cl)
    
    
    ## Integrer les solutions
    
    for(i in 1:n_taxa){
      if(is.null(Three_beta_pval[[i]])==TRUE){
        
        ESTIMATION6[i, k] <- NA ; PVALUE6[i, k] <- NA ; LOGFC6[i, k] <- NA
      }
      else {
        ESTIMATION6[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE6[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC6[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    pval_BB<-PVALUE6[,k]
    LAMBDA[k,6]<-median(qchisq(1-pval_BB,1),na.rm = "TRUE")/qchisq(0.5,1)
    
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_BB[1:(n0)], scores.class1 = 
                      pval_BB[(n0+1):n_taxa],curve = TRUE)
    AUC[k,6]<-proc$auc
    
    print(k)
  }
  
  
  
  

  
  
  return(list(LAMBDA=LAMBDA,AUC=AUC,
         ESTIMATION1=ESTIMATION1,ESTIMATION2=ESTIMATION2,ESTIMATION3=ESTIMATION3,
         ESTIMATION4=ESTIMATION4, ESTIMATION5=ESTIMATION5,ESTIMATION6=ESTIMATION6,
         PVALUE1=PVALUE1,PVALUE2=PVALUE2, PVALUE3=PVALUE3, PVALUE4=PVALUE4, PVALUE5=PVALUE5, PVALUE6=PVALUE6,
         LOGFC1=LOGFC1,LOGFC2=LOGFC2,LOGFC3=LOGFC3,LOGFC4=LOGFC4,LOGFC5=LOGFC5,LOGFC6=LOGFC6))
  
}



## Application of the multiple    ##
##
##
## ============================== ##

n_samples<-30
n_taxa<-3000
n_signif<-30
beta<-c(rnorm(n_signif,0,1), rep(0,n_taxa-n_signif)) ## n_signif est le nombre de genes significatif
alpha<-runif(1)
lambda_estimates=runif(n_taxa,1,10)
size<-sample(n_samples:(2*n_samples),n_samples)
mu_nb<-sample(n_taxa,n_taxa)
sigma<-0.4 
groups<-sample(rep(c("Non_Planted_OSPW", "ZCarex_OSPW"), each = n_samples/2)) 

RESULT<-simulation_NB(B=300, # nombre de simulation,
                      beta=beta,
                      alpha=alpha,
                      n_samples=n_samples, # Number of sample
                      n_taxa=n_taxa, # Number of taxa
                      lambda_estimates=lambda_estimates, # parametres de la poisson
                      size=size, # longeur de l'abondance (dim=n_sample multinomiale et size)
                      mu_nb=mu_nb, # moyenne de la Poisson qui est celui du Gamma
                      sigma=sigma, # probabilité de zero dans les données (dimemsion)
                      groups=groups)
save(list = ls(), file ="/Users/fatimaakpo/Downloads//simul_no_data_Small_Scenario_ST1.1.RData")

#boxplot(RESULT$LAMBDA[,1])    

load("/Users/fatimaakpo/Downloads/Student/simul_no_data_Small_Scenario_ST1.1.RData")


LAMBDA<-RESULT$LAMBDA
data_long <- melt(LAMBDA)

# Tracer le boxplot avec ggplot2

gplot1.1 <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y

gplot1.1


ggsave("/Users/fatimaakpo/Downloads/simul1_nodata_ST1.1.pdf",plot=gplot1.1,
       width = 4.2, height = 3.2, units = "in", dpi = 300,bg = "white")



## QQPlot de quelques valeurs.

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






## Load the data


LAMBDA<-RESULT$LAMBDA
pval<-matrix(0,nrow=n_taxa,ncol = 6)
colnames(pval)<-c("Poisson","NB","Hurdle","ZIP","DM","BM")


## pval 1 : Poisson distribution
ind<-which(round(LAMBDA[,1],2)==1.00)[1]#which((RESULT$LAMBDA[,1]<1.05)&(RESULT$LAMBDA[,1]>0.99))[1]
pval[,1]<-RESULT$PVALUE1[,ind]

## pval 2 : NB distribution
ind<-which(round(LAMBDA[,2],2)==1.00)[1] #which((RESULT$LAMBDA[,2]<1.05)&(RESULT$LAMBDA[,2]>0.99))[1]
pval[,2]<-RESULT$PVALUE1[,ind]

## pval 3 : Hurdle distribution
ind<-which(round(LAMBDA[,3],2)==min(round(LAMBDA[,3],2)))[1] # which((RESULT$LAMBDA[,3]<1.05)&(RESULT$LAMBDA[,3]>0.99))[1]
pval[,3]<-RESULT$PVALUE1[,ind]

## pval 4 : ZIP distribution
ind<-which(round(LAMBDA[,4],2)==1.00)[1] #which((RESULT$LAMBDA[,4]<1.05)&(RESULT$LAMBDA[,4]>0.99))[1]
pval[,4]<-RESULT$PVALUE1[,ind]

## pval 5 : DM distribution
ind<-which(round(LAMBDA[,5],2)==1.00)[1] #which((RESULT$LAMBDA[,5]<1.05)&(RESULT$LAMBDA[,5]>0.99))[1]
pval[,5]<-RESULT$PVALUE1[,ind]

## pval 6 : BM distribution
ind<-which(round(LAMBDA[,6],2)==1.00)[1] #which((RESULT$LAMBDA[,6]<1.05)&(RESULT$LAMBDA[,6]>0.99))[1]
pval[,6]<-RESULT$PVALUE1[,ind]


##=================================================================##
##
##
##
##
##====================== Tracer les fig ==========================##

library(tidyr)
library(dplyr)
data_qqplot<-as.matrix(pval)
rownames(data_qqplot)<-NULL


data_long <- as.data.frame(data_qqplot) %>%
  pivot_longer(everything(), names_to = "method", values_to = "pvalue")

# Ajouter une colonne pour les indices des méthodes
data_long$method_index <-data_long$method #as.factor(rep(1:ncol(data_qqplot), each = nrow(data_qqplot)))

data_long <- data_long %>%
  group_by(method) %>%
  arrange(pvalue) %>%
  mutate(index = row_number(), # Ajouter un index pour chaque p-value
         length_obs = n(), # Nombre total d'observations par méthode
         theoretical_quantile = -log10(index / length_obs)) %>%
  ungroup()
data_long <- data_long %>% data.frame()


# Créer le graphique
result_plot <- ggplot(data_long, aes(x = theoretical_quantile,
                                     y = -log10(pvalue), color = method)) +
  geom_point() + # Ajouter les points
  geom_line(aes(group = method)) + # Ajouter des lignes pour chaque méthode
  labs(x = "-log10(p) Expected", 
       y = "-log10(p) Observed", 
       title = "") +
  scale_color_discrete(name = "") +
  theme_pubr() +
  theme(legend.position = "top")+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed",lwd=1.2)+
  scale_linetype_manual(name = "y=x", values = c("dashed"))+
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))# A

# Afficher le graphique
print(result_plot)
ggsave("/Users/fatimaakpo/Downloads/qqplot_scenario_ST1.2.pdf",plot=result_plot,
       width = 4.2, height = 3, units = "in", dpi = 300,bg = "white")


###======= Calcul des biais ============ ###


function_RMSE <- function(x) { 
  sqrt(mean((x - beta)^2,na.rm=TRUE)) 
}

## Appliquer la fonction RMSE data


RMSE<-matrix(0,nrow=n_taxa,ncol = 6)
colnames(RMSE)<-c("Poisson","NB","Hurdle","ZIP","DM","BM")
RMSE[,1] <- apply(RESULT$ESTIMATION1, 2, function_RMSE)
RMSE[,2] <- apply(RESULT$ESTIMATION2, 2, function_RMSE)
RMSE[,3] <- apply(RESULT$ESTIMATION3, 2, function_RMSE)
RMSE[,4] <- apply(RESULT$ESTIMATION4, 2, function_RMSE)
RMSE[,5] <- apply(RESULT$ESTIMATION5, 2, function_RMSE)
RMSE[,6] <- apply(RESULT$ESTIMATION6, 2, function_RMSE)

## Transformation des variables 

data_long <- melt(RMSE)

## Tracer le boxplot avec ggplot2

gplot1.1_MSE <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "RMSE") +
  theme_pubr() +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5))+  # Centrer le titre
  coord_cartesian(ylim = c(0, 7000))  # Fixer les limites de l'axe y


gplot1.1_MSE


ggsave("/Users/fatimaakpo/Downloads/RMSE_scenario_ST1.1.pdf",plot=gplot1.1_MSE,
       width = 4, height = 2.8, units = "in", dpi = 300,bg = "white")


## R2 Lasso

#function_R2 <- function(pred) {
#  ss_res <- sum((beta - pred)^2)  # Somme des carrés des résidus
#  ss_tot <- sum((beta - mean(beta))^2)  # Somme totale des carrés
#  R2 <- 1 - (ss_res / ss_tot)
#  return(R2)
#}
## Appliquer la fonction RMSE data
function_R2<- function(col) cor(beta, col, use = "complete.obs")^2

R2_HMMDAY<-matrix(0,nrow=n_taxa,ncol = 6)
colnames(R2_HMMDAY)<-c("Poisson","NB","Hurdle","ZIP","DM","BM")
R2_HMMDAY[,1] <- apply(RESULT$ESTIMATION1, 2, function_R2)
R2_HMMDAY[,2] <- apply(RESULT$ESTIMATION2, 2, function_R2)
R2_HMMDAY[,3] <- apply(RESULT$ESTIMATION3, 2, function_R2)
R2_HMMDAY[,4] <- apply(RESULT$ESTIMATION4, 2, function_R2)
R2_HMMDAY[,5] <- apply(RESULT$ESTIMATION5, 2, function_R2)
R2_HMMDAY[,6] <- apply(RESULT$ESTIMATION6, 2, function_R2)

## Transformation des variables 

data_long <- melt(R2_HMMDAY)

## Tracer le boxplot avec ggplot2

gplot1.1_R2 <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "R2") +
  theme_pubr() +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5))+  # Centrer le titre
  coord_cartesian(ylim = c(0, 0.02))  # Fixer les limites de l'axe y


gplot1.1_R2

ggsave("/Users/fatimaakpo/Downloads/R2_scenario_ST1.1.pdf",plot=gplot1.1_R2,
       width = 4, height = 2.8, units = "in", dpi = 300,bg = "white")

