library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)

expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1,sep = "\t", header = T,check.names = FALSE)
SereData <- read.table("SERE.txt", row.names = 1,sep = "\t", header = T,check.names = FALSE)
Gene<- read.table("gene.txt",sep = "\t", header = T,check.names = FALSE,fill=TRUE)


#### Correspondance Gene ####
#Nom de gene différents entre article et sortie galaxy, correspondance effectuée sur Yeastmine mais deux noms de gènes
#modifiés DUR1,2 -> DUR12 02-oct -> OCT1

#renommer les colonnes par celle de l'article
colnames(expData)<-colnames(SereData)
#Correspance entre les noms de genes de l'alignement et de nom de l'article
expData[,ncol(expData)+1] <- Gene$input[match(rownames(expData), Gene$secondaryIdentifier)]
#Enlever les lignes sans noms pour renommer les lignes
expData2<-expData[complete.cases(expData), ]
rownames(expData2)<-expData2$V51
expData2<-expData2[,-51]


#### Z-Score Matrix ####


expData3<-expData2
expData2<-as.matrix(expData2)
for (i in 1:nrow(expData2)) {
  avg_var <- mean(expData2[i,])
  sd_var  <- sd(expData2[i,])
  for (j in 1:ncol(expData3)) {
    
    z_sc<-((expData2[i, j] - avg_var) / sd_var)
    expData3[i, j] <- z_sc
  }
}



#### Heatmap ####
####Brut####

#ajout d'un nom de colonne pour lignes (pour tidy/ggplot)
expData3<-rownames_to_column(as.data.frame(expData3))
#Fixer les lignes selon l'ordre de la publication pour le graphique et supprimer lignes vides
expData3<-expData3[match(rownames(SereData), expData3$rowname), ]
expData3<-expData3[complete.cases(expData3), ]
expData3$rowname <- factor(expData3$rowname,level=expData3$rowname)
#Transformer les colonnes en lignes
expData4<-pivot_longer(expData3,c(2:ncol(expData3)))
#Fixer les anciennes colonnes pour conserver leurs ordres lors du graphiques
expData4$name <- factor(expData4$name,level=unique(expData4$name))
#Graphique
ggplot(expData4,aes(x=name,y=rowname,fill=value))+geom_tile()+
  #geom_tile = heatmap
  scale_fill_gradient2(low = "#03f9f9", high = "#f4f402", mid = "#050302",midpoint = 0, limit = c(-1.5,1.5),oob=squish)+
  #choix couleur gradient et "oob" permet de gérer les outliers des l'échelle choisies
  scale_x_discrete(breaks=c("0","50","100","150","200"))+
  #couper l'échelle en valeurs voulues
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(face="italic"))+
  #Permet de modifier les paramêtres graphiques comme le titre et enlever les valeurs en Y trop nombreuses
  xlab("Times (minutes)")+ylab(paste("Top Periodic Genes (",nrow(expData3),")"))+ggtitle("Saccharomyces cerevisiae")

#ggsave("Sere.png", height=10, width=7)

#### Let's go Tidy ####

# la fonction rownames_to_column pourrait être aussi inclus mais pas de solution trouvé pour conserver l'ordre des y
expData3 %>%
  pivot_longer(c(2:ncol(.)))%>%
  ggplot(aes(x=fct_inorder(name),y=fct_inorder(rowname),fill=value))+geom_tile()+
  #fc_inorder permet de garder l'ordre dans la matrice (seulement en tidy)
  scale_fill_gradient2(low = "#03f9f9", high = "#f4f402", mid = "#050302",midpoint = 0, limit = c(-1.5,1.5),oob=squish)+
  scale_x_discrete(breaks=c("0","50","100","150","200"))+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(face="italic"))+
  xlab("Times (minutes)")+ylab(paste("Top Periodic Genes (1196)"))+ggtitle("Saccharomyces cerevisiae")
