library("dynamicTreeCut", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("ade4", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("seriation", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("gplots", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("DendSer", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("stringr", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("ape", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("FactoInvestigate", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("factoextra", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("FactoMineR", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("ape", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("ade4", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("vegan", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")
library("heatmaply", lib.loc="~/R/x86_64-pc-linux-gnu-library/4.0")

#Creation of ordered matrix with reference tree

#Creation of the matrix with species order according to tree
nomFi_matrix="matrice_fold_article_complet"
nomFi="./arbre_ref_binary_ultrametric.nwk"
arbre_binaire_ultra=read.tree(file=nomFi)


dendro=as.hclust(arbre_binaire_ultra) #YOUPI !!!

d=read.table(nomFi_matrix, header = T)


#Order the matrix according to the tree
matrice=matrix(data = 0, nrow=nrow(d), ncol=ncol(d), dimnames=list(rownames(d), dendro$labels))
i=1
for(org in dendro$labels){
  #print(org)
  matrice[,i] =d[,colnames(d)==org]
  i=i+1
}

#Create distance matrix to make the clustering method 7 is for ochiai distance
dMat2 <- dist.binary(matrice, method = 7, diag = FALSE, upper = FALSE)
TdMat2 <- dist.binary(t(as.matrix(matrice)), method = 7, diag = FALSE, upper = FALSE)

#compute the hierarchical clustering with mcquitty method
HC2 <- hclust(dMat2, method = "mcquitty", members=NULL)
#THC2_clust <- hclust(TdMat2, method = "mcquitty", members=NULL)
#write.nexus(as.phylo(THC2_clust), file="clustering.nex")


#Seriation on the hierarchical clustering for fold and on the reference tree for species
OrdSer1 <- DendSer(HC2,dMat2, cost = costARc)
OrdSer2 <- DendSer(dendro,TdMat2, cost = costARc)



#reorder the matrix according to the seriation.
datamatOrdered <- as.matrix(matrice[get_order(OrdSer1), get_order(OrdSer2)])
heatmap.2(datamatOrdered,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')




#Create a vector to choose color of clusters for matalt
vectcolor=function(ClusterVec,OrdSer1){
  a=ClusterVec[get_order(OrdSer1)]
  a=as.vector(a)
  newvect=vector()
  number=1
  actu=a[1]
  for (i in a){
    if (i==0){
      newvect=append(newvect,0)
    }
    else{
      if( i==actu){
        newvect=append(newvect,number)}
      else{
        if (number==1){
          number=2
          newvect=append(newvect,number)
          actu=i}
        else{number=1
        newvect=append(newvect,number)
        actu=i
        }}}}
  return (newvect)
}
newvect=vectcolor(ClusterVec,OrdSer1)







#creation d'une matrice alternée
matalt=function(matrix_fold_spe_full,newvect){
  vecteurside=as.vector(newvect[get_order(OrdSer1)])
  couleurintercale=c("#FFFFFF","#333333","#FF3300","#3333FF")
  i=1
  matrix2=as.matrix(matrix_fold_spe_full[get_order(OrdSer1), get_order(OrdSer2)])
  while (i<(length(newvect)+1)){
    for (k in 1:length(matrix2[i,])){
      if (matrix2[i,k]==1){
        matrix2[i,k]=1+newvect[i]
      }
    }
    i=i+1
  }
  return (matrix2)
}
#utilisation fonction creation matrice alterné
matrixalt=matalt(matrice,newvect )





#PCA
#Utilisation de la pca
data2=t(matrice)
res.pca2=PCA(data2)
variablepca=get_pca_var(res.pca2)

matricepca=datamatOrdered
coordpca=variablepca$coord
coordpcaord=coordpca[get_order(OrdSer1), ]
contribpca=variablepca$contrib
contribpcaord=contribpca[get_order(OrdSer1), ]

for (i in 1:length(coordpca[,1])){
  a=sqrt(coordpcaord[i,"Dim.1"]^2+coordpcaord[i,"Dim.2"]^2)
  if (a>0.80){
    if (coordpcaord[i,1]>0.80){
      valeurmat=2
    }else{
      if (coordpcaord[i,1]>0 & coordpcaord[i,2]>0){
        valeurmat=3
      }
      if (coordpcaord[i,1]>0 & coordpcaord[i,2]<0){
        valeurmat=4
      }
      if (coordpcaord[i,1]<0 & coordpcaord[i,2]>0){
        valeurmat=5
      }
    }
    for (j in 1:210){
      if (matricepca[i,j]==1){
        matricepca[i,j]=valeurmat
      }
    }
  }
}

heatmap.2(matricepca,dendrogram='none',col=c("white","#999999","#0066CC","#CC3366","#660066","#FF6633"), Rowv=FALSE, Colv=FALSE,trace='none')



plot.new()
plot(1, type="n", xlab="", ylab="", xlim=c(-1, 1), ylim=c(-1, 1),asp=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
draw.circle(0,0,1)


for (i in 1:length(coordpca[,1])){
  a=sqrt(coordpcaord[i,"Dim.1"]^2+coordpcaord[i,"Dim.2"]^2)
  coordx=coordpcaord[i,"Dim.1"]
  coordy=coordpcaord[i,"Dim.2"]
  if (a>0.80){
    if (coordpcaord[i,1]>0.75){
      valeurmat=2
      arrows(0,0,coordx,coordy,col="#0066CC")
    }else{
      if (coordpcaord[i,1]>0 & coordpcaord[i,2]>0){
        valeurmat=3
        arrows(0,0,coordx,coordy,col="#CC3366")
      }
      if (coordpcaord[i,1]>0 & coordpcaord[i,2]<0){
        valeurmat=4
        arrows(0,0,coordx,coordy,col="#660066")
      }
      if (coordpcaord[i,1]<0 & coordpcaord[i,2]>0){
        valeurmat=5
        arrows(0,0,coordx,coordy,col="#FF6633")
      }
    }
    for (j in 1:210){
      if (matricepca[i,j]==1){
        matricepca[i,j]=valeurmat
      }
    }
  }
}


liste=c(1,2,3,1,2)
hist(liste)


