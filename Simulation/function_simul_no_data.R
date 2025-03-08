## Load the packages

pacman::p_load(ape, vegan, limma, statmod, edgeR, ggplot2, Rcpp, devtools, metaSPARSim, cowplot, plyr, ggpubr, PRROC, gamlss.dist, 
               doParallel, reshape2, MCMCpack, variancePartition, BiocParallel, tidyr, dplyr)

## Create the function of load the data

## Create the function
## Cette simule des donnees d'abondance X suivant plusieurs distribution discrète (ZIP, NB, etc.), calcul y et fixe des paramètres beta
## Input : 

## Output :

simulation_NB<-function(B, # Number of replication,
                        beta, alpha, # parameters of modelling 
                        n_samples, # Sample size
                        n_taxa, # Number of metagens
                        lambda_estimates, # Parameters of distribution
                        size, # Size of abundance (dim=n_sample multinomiale et size)
                        mu_nb, # Mean of Poisson et Gamma distributions
                        sigma, # Probability of zero
                        groups){

  ## Nombre de gènes significatif(tau=t)
  n0<-length(which(beta!=0)) # Number of significatif genes pour les gènes

  
  ## Matrice des facteurs d'impact (IF)
  LAMBDA<-  LAMBDA_X<-LAMBDA_Y<-LAMBDA_Z<-matrix(0,nrow=B,ncol = 4)

  colnames(LAMBDA)<-colnames(AUC)<-c("NB","Hurdle","ZIP","NB_Hurdle_ZIP")
  colnames(LAMBDA_X)<-c("Edge_NB","Edge_Hurdle","Edge_ZIP","Edge_NB_Hurdle_ZIP")
  colnames(LAMBDA_Y)<-c("Limma_NB","Limma_Hurdle","Limma_ZIP","Limma_NB_Hurdle_ZIP")
  colnames(LAMBDA_Z)<-c("Dream_NB","Dream_Hurdle","Dream_ZIP","Dream_NB_Hurdle_ZIP")

   ## Parameters estimates
  ESTIMATION1<-ESTIMATION2<-ESTIMATION3<-ESTIMATION4<-
    PVALUE1<-PVALUE2<-PVALUE3<-PVALUE4<-
    LOGFC1<-LOGFC2<-LOGFC3<-LOGFC4<-
  LOGFC_NB_EDGE<-LOGFC_Hurdle_EDGE<-LOGFC_ZIP_EDGE<-LOGFC_Mel_EDGE<-
  LOGFC_NB_Limma<-LOGFC_Hurdle_Limma<-LOGFC_ZIP_Limma<-LOGFC_Mel_Limma<-
  LOGFC_NB_Dream<-LOGFC_Hurdle_Dream<-LOGFC_ZIP_Dream<-LOGFC_Mel_Dream<-matrix(0,ncol=B,nrow = n_taxa)
  

  ## La boucle de la replication
  for (k in 1:B){
    
   
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
    X_NB <- generate_nb_counts(n_taxa, n_samples, 
                                      size=size_nb, mu=mu_nb) #ligne =taxa, colonne=echantillon
    
    
    m.D <- vegdist(t(X_NB), "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    
    ##===========Calculate Y===================## 
    
    Xbeta<-t(X_NB)%*%beta
    eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(X_NB[i,]+1)
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
        
        ESTIMATION1[i, k] <- NA ; PVALUE1[i, k] <- NA ; LOGFC1[i, k] <- NA
      }
      else {
        ESTIMATION1[i,k] <- Three_beta_pval[[i]][[1]] ; 
        PVALUE1[i,k] <- Three_beta_pval[[i]][[2]] ; 
        LOGFC1[i,k] <- Three_beta_pval[[i]][[3]]
      }
      
    }
    
    LAMBDA[k,1]<-median(qchisq(1-PVALUE1[,k],1),na.rm = "TRUE")/qchisq(0.5,1)
    
    ##=======================================
    #
    #
    # Application de  model "edge"
    
    metaData<-data.frame(groups=groups)
    
    # Préparation des données
    d  <- DGEList(counts = X_NB, group = metaData$groups)
    y_d <- calcNormFactors(d)
    
    # Définition du design
    design <- model.matrix(~ groups)
    y_estimate <- estimateDisp(y_d, design)
    
    # Test de l'expression différentielle
    fit <- glmQLFit(y_estimate, design)
    result <- glmQLFTest(fit, coef = 2)
    
    
    ## take the FDR
    
    pval<-result$table$PValue
    LAMBDA_X[k,1]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_NB_EDGE[,k]<-result$table$logFC
    
    
    ## Limma application
    ##
    ##===============================##
    

    #d <- DGEList(X_NB, group = metaData$groups)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y, design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    result<-topTable(out_ebayes,number = Inf)
    
    
    # Take the result
    
    LAMBDA_Y[k,1]<-median(qchisq(1-result$P.Val,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_NB_Limma[,k]<-result$logFC

    
    ## Model Dream
    
    # estimate weights using linear mixed model of dream
    
    vobjDream <- voomWithDreamWeights(y_d, formula =~groups, data=metaData)
    
    # otherwise it uses the Satterthwaite approximation
    fitmm <- dream(vobjDream, form=~groups, data=metaData)
    fitmm <- eBayes(fitmm)
    
    # Extracte result
    result<-topTable(fitmm, coef = 2, number = Inf)
    
    # Take the result
    pval_dream<-result$P.Value
    LAMBDA_Z[k,1]<-median(qchisq(1-pval_dream,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_NB_Dream[,k]<-result$logFC
        
    ## DISTRIBUTION 2 :::: ----- Distribution de ""Hurdle""
    
    ## Generate the data presence/absence
    prob_zero <- 0.3 # Probability que le comptage soit zéro
    n_samples<-n_samples ; n_taxa<-n_taxa
    presence_absence <- rbinom(n_samples * n_taxa, 1, 1 - prob_zero)
    presence_absence_matrix <- matrix(presence_absence, nrow = n_samples,
                                      ncol = n_taxa)
    
    ## Generate the count of non-nuls
    lambda_mean<-lambda_estimates
    X_Hurdle <- matrix(0, nrow = n_samples, ncol = n_taxa)
    for (j in 1:n_taxa) {
      X_Hurdle[, j] <- rpois(n_samples, lambda = lambda_mean) * 
        presence_absence_matrix[, j]
    }
    
    m.D <- vegdist(X_Hurdle, "manhattan") ; 
    result_GENERAL <- pcoa(m.D )
    
    # ================================================##
    
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    
    ##===========Calculate Y===================## 
    
    Xbeta<-X_Hurdle%*%beta
    eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(X_Hurdle[,i]+1)
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
    
    pval_Hurdle<-PVALUE2[,k]
    LAMBDA[k,2]<-median(qchisq(1-pval_Hurdle,1))/qchisq(0.5,1)
    
    ##=======================================================================##
    #
    #
    # Application de edge
    # Préparation des données
    #========================================================================##
    d <- DGEList(counts = t(X_Hurdle), group = metaData$groups)
    y_d <- calcNormFactors(d)
    
    # Définition du design
    y_estimate <- estimateDisp(y_d, design)
    
    # Test de l'expression différentielle
    fit <- glmQLFit(y_estimate, design)
    result <- glmQLFTest(fit, coef = 2)
    
    
    ## take the FDR
    
    pval<-result$table$PValue
    LAMBDA_X[k,2]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Hurdle_EDGE[,k]<-result$table$logFC
    
    ## Model : "Limma" application
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y, design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    result<-topTable(out_ebayes,number = Inf)
    
    
    # Take the result
    pval<-result$P.Val
    LAMBDA_Y[k,2]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Hurdle_Limma[,k]<-result$logFC
    
    # estimate weights using linear mixed model of dream
    
    
    vobjDream <- voomWithDreamWeights(y_d, 
                                      formula =~groups, 
                                      data=metaData)
    
    # otherwise it uses the Satterthwaite approximation
    fitmm <- dream(vobjDream, 
                   form=~groups,
                   data=metaData)
    fitmm <- eBayes(fitmm)
    
    #
    result<-topTable(fitmm, coef = 2, number = Inf)
    
    # Take the result
    pval_dream<-result$P.Value
    LAMBDA_Z[k,2]<-median(qchisq(1-pval_dream,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Hurdle_Dream[,k]<-result$logFC    
    
    ## Zero-Inflated Poisson (ZIP)
    #pi_inflate <- 0.3
    #n <-n_samples*n_taxa
    # Créer la matrice de microbiome simulée
    X_ZIP <- matrix(0, nrow = n_taxa, ncol = n_samples)
    
    ## Remplir la matrice avec les données ZIP
    
    for (i in 1:n_taxa) {
      X_ZIP[i, ] <- rZIP(n_samples, mu =mu_nb[i] , 
                                sigma = 0.3) # sum(gen_data[i,] < 20) / n_samplesrpois(1, lambda_est)
    }
    
    m.D <- vegdist(t(X_ZIP), "manhattan") ; 
    result_GENERAL <- pcoa(m.D)
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-t(X_ZIP)%*%beta
    eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(X_ZIP[i,]+1)
        #print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        #print(colnames(design))
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
    
    pval_ZIP<-PVALUE3[,k]
    LAMBDA[k,3]<-median(qchisq(1-pval_ZIP,1))/qchisq(0.5,1)
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_ZIP[1:(n0)], scores.class1 = 
                      pval_ZIP[(n0+1):n_taxa],curve = TRUE)
    AUC[k,3]<-proc$auc
    
    ##=======================================
    #
    #
    # Application de edge
    # Préparation des données
    d <- DGEList(counts = X_ZIP, group = metaData$groups)
    y_d <- calcNormFactors(d)
    
    # Définition du design
    y_estimate <- estimateDisp(y_d, design)
    
    # Test de l'expression différentielle
    fit <- glmQLFit(y_estimate, design)
    result <- glmQLFTest(fit, coef = 2)
    
    
    ## take the FDR
    
    pval<-result$table$PValue
    LAMBDA_X[k,3]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_ZIP_EDGE[,k]<-result$table$logFC

    
    ## Limma application ==================================
    
    #d <- DGEList(X_NB, group = metaData$groups)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y, design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    result<-topTable(out_ebayes,number = Inf)
    
    
    # Take the result
    pval<-result$P.Val
    LAMBDA_Y[k,3]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_ZIP_Limma[,k]<-result$logFC
    
    
    ## Model Dream

    #d  <- DGEList(counts = X_NB, group = metaData$groups)
    #y_d <- calcNormFactors(d)
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    #param <- SnowParam(4, "SOCK", progressbar = TRUE)
    
    # estimate weights using linear mixed model of dream
    
    
    vobjDream <- voomWithDreamWeights(y_d, 
                                      formula =~groups, 
                                      data=metaData)
    
    # otherwise it uses the Satterthwaite approximation
    fitmm <- dream(vobjDream, 
                   form=~groups,
                   data=metaData)
    fitmm <- eBayes(fitmm)
    
    #
    result<-topTable(fitmm, coef = 2, number = Inf)
    
    # Take the result
    pval_dream<-result$P.Value
    LAMBDA_Z[k,3]<-median(qchisq(1-pval_dream,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_ZIP_Dream[,k]<-result$logFC    
    
    #print("Deam3")
    
    
    
    
    

    ###=============================================================
    ##
    ## DIstribution de melange "NB", "Hurdle" et "ZIP"
    X_NB<-t(X_NB)
    X_Hurdle<-X_Hurdle
    X_ZIP<-t(X_ZIP)
    
    indice_gene<-1:n_taxa
    gen_nb<-sample(indice_gene, 1000) # les genes tirés dans la premieres matrices avec Distribution NB

    vect_res<-indice_gene[-gen_nb]

    gen_Hurdle_sample<-sample(vect_res, 1000)

    gen_ZIP<-indice_gene[-c(gen_nb,gen_Hurdle_sample)]

    data_NB<-X_NB[,gen_nb]
    data_Hurdle<-X_Hurdle[,gen_Hurdle_sample]
    data_ZIP<-X_NB[,gen_ZIP]
    
    sim_data_melang<- cbind(data_NB, data_Hurdle, data_ZIP)

    m.D <- vegdist(sim_data_melang, "manhattan") ; 
    result_GENERAL <- pcoa(m.D)
    
    ##=================================================##
    data<-data.frame(groups=groups)
    data$PCoA1<-result_GENERAL$vectors[,"Axis.1"]
    data$PCoA2<-result_GENERAL$vectors[,"Axis.2"]
    
    ##===========Calculate Y===================## 
    
    Xbeta<-sim_data_melang%*%beta
    eps<-rnorm(n_samples,0,1)
    Y<-Xbeta+alpha*as.numeric(factor(groups))+eps
    data$Y<-as.numeric(Y)
    sim_data_melang<-t(sim_data_melang)
    ## Design of the model
    
    cl <- detectCores() %>% -1 %>% makeCluster ### N'utiliser pas tous les coeurs (-1)
    registerDoParallel(cl)
    
    system.time(
      Three_beta_pval<-foreach( i= 1:n_taxa,.packages=c("limma")) %dopar% {
        #for(i in 1:n){
        data$Xg<-log10(sim_data_melang[i,]+1)
        #print(dim(data))
        design <- model.matrix(~Xg+groups+PCoA1+PCoA2,data)
        #print(colnames(design))
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
    
    pval_melang<-PVALUE4[,k]
    LAMBDA[k,4]<-median(qchisq(1-pval_melang,1),na.rm = "TRUE")/qchisq(0.5,1)
    ## Calcul for AUC for this courbe
    #labels <- c(rep(1,n_taxa-n0),rep())
    
    proc<-roc.curve(scores.class0 = pval_melang[1:(n0)], scores.class1 = 
                      pval_melang[(n0+1):n_taxa],curve = TRUE)
    AUC[k,4]<-proc$auc
    
    ##=======================================
    #
    #
    # Application de edge
    # Préparation des données
    d <- DGEList(counts = sim_data_melang, group = metaData$groups)
    y_d <- calcNormFactors(d)
    
    # Définition du design
    y_estimate <- estimateDisp(y_d, design)
    
    # Test de l'expression différentielle
    fit <- glmQLFit(y_estimate, design)
    result <- glmQLFTest(fit, coef = 2)
    
    
    ## take the FDR
    
    pval<-result$table$PValue
    LAMBDA_X[k,4]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Mel_EDGE[,k]<-result$table$logFC
    
    ## Limma application

    #d <- DGEList(X_NB, group = metaData$groups)
    y <- voom(d, design, plot = FALSE)
    limma_fit <- lmFit(y, design)
    tmp <- contrasts.fit(limma_fit, coef = 2) 
    out_ebayes <- eBayes(tmp)
    
    result<-topTable(out_ebayes,number = Inf)
    
    
    # Take the result
    pval<-result$P.Val
    LAMBDA_Y[k,3]<-median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Mel_Limma[,k]<-result$logFC
    
    
    ## Model Dream

    #d  <- DGEList(counts = X_NB, group = metaData$groups)
    #y_d <- calcNormFactors(d)
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    #param <- SnowParam(4, "SOCK", progressbar = TRUE)
    
    # estimate weights using linear mixed model of dream
    
    
    vobjDream <- voomWithDreamWeights(y_d, 
                                      formula =~groups, 
                                      data=metaData)
    
    # otherwise it uses the Satterthwaite approximation
    fitmm <- dream(vobjDream, 
                   form=~groups,
                   data=metaData)
    fitmm <- eBayes(fitmm)
    
    #
    result<-topTable(fitmm, coef = 2, number = Inf)
    
    # Take the result
    pval_dream<-result$P.Value
    LAMBDA_Z[k,4]<-median(qchisq(1-pval_dream,1),na.rm = "TRUE")/qchisq(0.5,1)
    LOGFC_Mel_Dream[,k]<-result$logFC    
    
    #print("Deam4")
    
    
    
    print(k)
  }
  

  
  
  
  return(list(LAMBDA=LAMBDA, AUC=AUC,
              ESTIMATION1=ESTIMATION1,ESTIMATION2=ESTIMATION2,ESTIMATION3=ESTIMATION3,
              ESTIMATION4=ESTIMATION4,
              PVALUE1=PVALUE1,PVALUE2=PVALUE2, PVALUE3=PVALUE3, PVALUE4=PVALUE4,
              LOGFC1=LOGFC1,LOGFC2=LOGFC2,LOGFC3=LOGFC3,LOGFC4=LOGFC4,
              ##Ajouter les resultats du model Edge
              LAMBDA_X=LAMBDA_X, LAMBDA_Y=LAMBDA_Y, LAMBDA_Z=LAMBDA_Z,
              LOGFC_NB_EDGE=LOGFC_NB_EDGE,LOGFC_Hurdle_EDGE=LOGFC_Hurdle_EDGE,
              LOGFC_ZIP_EDGE=LOGFC_ZIP_EDGE, LOGFC_Mel_EDGE=LOGFC_Mel_EDGE,
              
              LOGFC_NB_Limma=LOGFC_NB_Limma, LOGFC_Hurdle_Limma=LOGFC_Hurdle_Limma,
              LOGFC_ZIP_Limma=LOGFC_ZIP_Limma, LOGFC_Mel_Limma=LOGFC_Mel_Limma,
              LOGFC_NB_Dream=LOGFC_NB_Dream,LOGFC_Hurdle_Dream=LOGFC_Hurdle_Dream,
              LOGFC_ZIP_Dream=LOGFC_ZIP_Dream, LOGFC_Mel_Dream=LOGFC_Mel_Dream
              ))
  
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


save(list = ls(), file ="/Users/fatimaakpo/Desktop/Simulation-melange/simul_no_data_Small_Scenario-melange1.1.RData")

#boxplot(RESULT$LAMBDA[,1])    



## Model HMMDAY

LAMBDA_X<-RESULT$LAMBDA_X
data_long <- melt(LAMBDA_X)

# Tracer le boxplot avec ggplot2

gplot1.1_X <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1.1_X


ggsave("/Users/fatimaakpo/Desktop/Simulation-melange/simul1_nodata_1.1_edge.pdf",plot=gplot1.1_X,
       width = 4.2, height = 3.2, units = "in", dpi = 300,bg = "white")

# Tracer le boxplot avec ggplot2 "Limma"
LAMBDA_Y<-RESULT$LAMBDA_Y
data_long <- melt(LAMBDA_Y)

gplot1.1_Y <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1.1_Y


ggsave("/Users/fatimaakpo/Desktop/Simulation-melange/simul1_nodata_1.1_edge.pdf",plot=gplot1.1_Y,
       width = 4.2, height = 3.2, units = "in", dpi = 300,bg = "white")


# Tracer le boxplot avec ggplot2
LAMBDA_Z<-RESULT$LAMBDA_Z
data_long <- melt(LAMBDA_Z)

gplot1.1_Z <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1.1_Z


## Model HMMDAY

LAMBDA<-RESULT$LAMBDA
colnames(LAMBDA)<-c("HMMDAY_NB","HMMDAY_Hurdle","HMMDAY_ZIP","HMMDAY_NB_Hurdle_ZIP")

data_long <- melt(LAMBDA)
# Tracer le boxplot avec ggplot2

gplot1 <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1


save(list = ls(), file ="/Users/fatimaakpo/Downloads//simul_no_data_Scenario_melagane1.1.RData")

AUC<-RESULT$AUC
data_long <- melt(AUC)


# Tracer le boxplot avec ggplot2

gplot2<-ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = "", x = "", y = "AUC") +
  theme_pubr()+
  scale_fill_discrete(guide = "none") 
gplot2



## QQPlot de quelques valeurs.

####qqplot comparison---------------

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



## ============================== ##
##
##
##
##. DESEQUILIBRE



n_samples<-30  # taille d'échantillon
n_taxa<-3000 # nombre de taxa
n_signif<-30
beta<-c(rnorm(n_signif,0,1), rep(0,n_taxa-n_signif)) ## n_signif est le nombre de genes significatif
alpha<-runif(1)
lambda_estimates=runif(n_taxa,1,10)
size<-sample(n_samples:(2*n_samples),n_samples)
mu_nb<-sample(n_taxa,n_taxa)
sigma<-0.4 
groups<-sample(c("Non_Planted_OSPW", "ZCarex_OSPW"), size = n_samples,replace = TRUE) 

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




save(list = ls(), file ="/Users/fatimaakpo/Desktop/Simulation-melange/simul_no_data_Desequilibre_Scenario-melange1.1.RData")


##--------------------------------------#
#
#
#--------------------------------------#
load("/Users/fatimaakpo/Downloads/Simulation_Part2/simul_no_data_Desequilibre_Scenario1.1.RData")

LAMBDA<-RESULT$LAMBDA
data_long <- melt(LAMBDA)

# Tracer le boxplot avec ggplot2

gplot1.1_des <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1.1_des

ggsave("/Users/fatimaakpo/Downloads/simul1_nodata_deseq_1.1.pdf",plot=gplot1.1_des,
       width = 4.2, height = 3.2, units = "in", dpi = 300,bg = "white")



load("/Users/fatimaakpo/Downloads/Simulation_Part2/simul_no_data_Desequilibre_Scenario1.1.RData")

LAMBDA_X<-RESULT$LAMBDA_X
data_long <- melt(LAMBDA_X)

# Tracer le boxplot avec ggplot2

gplot1.1_des_edge <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "Impact Factor") +
  theme_pubr() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", lwd = 1) +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +  # Centrer le titre
  coord_cartesian(ylim = c(0, 1.8))  # Fixer les limites de l'axe y


gplot1.1_des_edge

##=========================================##

X_NB<-X_NB;
metaData<-data.frame(groups=groups)
library(edgeR)
# Préparation des données
y <- DGEList(counts = X_NB, group = metaData$groups)
y <- calcNormFactors(y)

# Définition du design
groups = metaData$groups
design <- model.matrix(~ groups)

# Estimation de la dispersion
y <- estimateDisp(y, design)

# Test de l'expression différentielle
fit <- glmQLFit(y, design)
result <- glmQLFTest(fit, coef = 2)


## take the FDR

pval<-result$table$PValue
logFC<-result$table$logFC
FDR <- p.adjust(pval, method="BH") # Benjamini-Hojberg
boxplot(FDR,col="skyblue", outline = FALSE)
median(qchisq(1-pval,1),na.rm = "TRUE")/qchisq(0.5,1)


v_edge<-volcanoplot(x=logFC,y=-log10(pval),tit="",col="#FF69B4")
myqqplot2(pval)
v_edge

## RMSE

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

gplot_des_1.1_MSE <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "RMSE") +
  theme_pubr() +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5))+  # Centrer le titre
  coord_cartesian(ylim = c(0, 6000))  # Fixer les limites de l'axe y


gplot_des_1.1_MSE


ggsave("/Users/fatimaakpo/Downloads/RMSE_scenario_des_1.1.pdf",plot=gplot_des_1.1_MSE,
       width = 4.2, height = 3, units = "in", dpi = 300,bg = "white")


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
gplot_des_1.1_R2 <- ggplot(data_long, aes(x = Var2, y = value, fill = Var2)) + 
  geom_boxplot(outlier.shape = NA) +
  labs(title = expression("" * tau * " = 0.01, " * n/p * " = 1%"), 
       x = "", 
       y = "R2") +
  theme_pubr() +
  scale_fill_discrete(guide = "none") +
  theme(plot.title = element_text(hjust = 0.5))+  # Centrer le titre
  coord_cartesian(ylim = c(0, 0.02))  # Fixer les limites de l'axe y


gplot_des_1.1_R2

ggsave("/Users/fatimaakpo/Downloads/R2_scenario_des_1.1.pdf",plot=gplot_des_1.1_R2,
       width = 4.2, height = 3, units = "in", dpi = 300,bg = "white")






## load the packages

library(ape) ; library(vegan) ; library(limma) ; library(statmod) ; library(edgeR)
library(ggplot2) ; library(Rcpp) ; library(devtools) ; library(metaSPARSim)
library(cowplot) ; library(plyr) ; library(ggpubr) ; library(PRROC) ; library(MCMCpack)
library(doParallel) ;   library(gamlss.dist)



