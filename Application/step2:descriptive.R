##----------------------- Load the packages -------------------------------------------##

pacman::p_load(ape, BiocParallel, caret, cowplot, devtools, doParallel, dplyr, edgeR, gamlss.dist, ggplot2, ggpubr, ggvenn, haven, lattice, limma, 
               MASS, MCMCpack, metaSPARSim, plyr, PRROC, Rcpp, reshape2, statmod, tidyr, variancePartition, vegan, VennDiagram)

###-------------- Recharger le fichier des resultats precedents ------------------------##

load("/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/data_sediment_complet.RData")

n<-ncol(GENOM_S)-4  ## 4 variables sont supprimes : Time, Meso, Groups et NA_pmm. Dans GENOM_S "Compartment" est deja supprimer a cause de 
print(n) ;
print(GENOM_S[,n+1])

###----- Data final using------------------------------

library(doParallel)
cl <- detectCores() %>% -1 %>% makeCluster  ###n'utiliser pas tous les coeurs (-1)
registerDoParallel(cl)

system.time(
  DES_DATA<-foreach( i= 1:n) %dopar% {
    maxim<-max(GENOM_S[,I])
    media<-median(GENOM_S[,I])
    sde<-sd(GENOM_S[,I])
    print(c(maxim,media,sde))
}
)

## Arreter le cluster de calcul

stopCluster(cl)



##--------Reunir les resultats------------##

maxim<-vector(mode = "numeric", length = n) ; media<-vector(mode = "numeric", length = n) ; sde<-vector(mode = "numeric", length = n)

for(i in 1:n){
              maxim[i]<-DES_DATA[[i]][[1]] ;
              media[i]<-DES_DATA[[i]][[2]] ;
              sde[i]<-DES_DATA[[i]][[3]]
   }



par(mfrow=c(1,2))
hist(media,100,col="skyblue",xlab="Maximum", ylab="Frequence")
hist(log10(media+1),30,col="skyblue",xlab="Log (maximum)", ylab="Frequence")


########--Construction des graphiques---------------------###


##maximum
pmax<-ggplot() +
      geom_boxplot(aes(y = maxim), fill = "lightblue", color = "black", alpha = 0.7, outlier.color = "red") +
      labs(x = "", y = "Valeur ppm") +
      ggtitle("Boxplot des valeurs maximales")

ggsave("/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/BOXPLOY_MAX.jpg", plot=pmax)

###mediane

pmedia<-ggplot() +
geom_boxplot(aes(y = media), fill = "lightblue", color = "black", alpha = 0.7, outlier.color = "red") +
  labs(x = "", y = "Valeur ppm") +
  ggtitle("Boxplot des valeurs medianes des genes")

ggsave("/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/BOX_MEDIAN.jpg", plot=pmedia)

####--Standard----
psd<-ggplot() +
     geom_boxplot(aes(y = sde), fill = "lightblue", color = "black", alpha = 0.7, outlier.color = "red") +
     labs(x = "", y = "Valeur ppm") +
     ggtitle("Boxplot des ecart-types")


ggsave("/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/BOX_SD.jpg", plot=psd)



##----Grid of the plot-------

grid<-plot_grid(pmax, pmedia,psd, ncol = 3)
ggsave("/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/GRILL.jpg", plot=gris)


## Save the rdata

save.image(file = "/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/Descriptive_sediment_Complet.RData")


