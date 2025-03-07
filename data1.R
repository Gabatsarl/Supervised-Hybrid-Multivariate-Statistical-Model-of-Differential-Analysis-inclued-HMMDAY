##---------------------------------------------------------------------------------------------##
#
#
#----------------------------------------------------------------------------------------------##

library(elasticnet) ; library(kernlab) ; library(xgboost) ; library(glmnet) ; library(RSNNS); library(MASS) ; library(dplyr) ; 
require(haven) ; require(ggplot2) ; library (reshape2) ; require(lattice) ; require(caret) ; library(gbm) ; library(e1071) 

###------ Charger les donnees metagenomic --------------###

XX= read.csv('/lustre03/project/6074237/gabat/GENOM/merged_gene_abundance.tsv', sep='\t', header = TRUE, row.names = 1, comment.char="*")

XXT<-XX
XXT<-as.data.frame(t(XXT))
row.names(XXT)<-NULL
XXT$ID=colnames(XX)

rm(XX)



##------load NA Degradation data------------------------###

mapping_TRAIN <- read.csv("/lustre03/project/6074237/gabat/GENOM/mapping_TRAIN.tsv", sep=";", comment.char="*", header=TRUE)


###---Correction d'une variables sous------------------------------------------

groups<-NA

groups[which((mapping_TRAIN$Sample_type=="No_plant")&(mapping_TRAIN$Water_type=="OSPW"))]<-"OSPW"

groups[which(mapping_TRAIN$Sample_type=="Carex")]<-"OSPW+Carex"
groups[which((mapping_TRAIN$Sample_type=="No_plant")&(mapping_TRAIN$Water_type=="Artificial_OSPW"))]<-"Control"


mapping_TRAIN$groups<-groups



####----------supprimer les colonnes non moins importantes -------##

mapping<-mapping_TRAIN[,-c(2,5,6,8)]

mapping<-mapping[,c(1:3,6,4,5)]
names(mapping)[1] <- "ID"
#mapping<-mapping[which(mapping$groups!="Control"),]



##-----------Ajouter des variables-------------##
##prendre une partie des donnees

g1<- gsub("\\.\\.\\.", "---", XXT[,ncol(XXT)])  ##remplacer le ... de ID par ---

g1 <- unlist(lapply(strsplit(g1, "---"), function(x) x[1]))
g2<-unlist(lapply(strsplit(mapping[,1], "---"), function(x) x[1]))




###-----Remplacer les colonnes par l'importance--------##
XXT$ID<-g1 ; mapping$ID<-g2


##---Merger les deux datas-----------------------------##
GENOM_S<-merge(XXT,mapping,by="ID")
GENOM_S$ID<-NULL


##-----Supprimer les variables avec des zeros partout--------##
remove <- nearZeroVar(GENOM_S)
GENOM_S <- GENOM_S[, -remove]


##prendre ce qui est important---------

#GENOM_S<-GENOM[which(GENOM$Compartment=="Sediments"),]
#GENOM_S<- subset(GENOM_S, select = -c(Compartment))
#GENOM_S<-GENOM_S[which(GENOM_S$Time!="D0"),]

##---Afficher le nombre de variable vide, nulle ou constante-----## 
print(length(remove))
print(dim(GENOM_S))

###------Supprimer des data inutiles--------------##
rm(remove)
rm(XXT)
rm(mapping)
rm(mapping_TRAIN)



###-------ENREGISTRER LES DONNEE SEDIMENTS-----------##

save.image(file = "/home/gabat/scratch/GENOM/LOGICIEL-R/STEP0/SOLUTION/data1_complet.RData")


