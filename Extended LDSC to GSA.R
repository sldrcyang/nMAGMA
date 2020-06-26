setwd("/your path/")
g=read.table("your .genes.out file", header = T,sep = "\t",stringsAsFactors = F)
gs=readLines("your gene set annot file") #with each row corresponding to a gene set:name of the gene set followed by the gene IDs
matrix=data.frame(g$GENE,stringsAsFactors = F)

library(dplyr)
geneset <- c()
genelist <- list()
generows <- list()

#Create a gene-geneset matrix
for(i in 1:length(gs)){
  temp <- unlist(strsplit(gs[i], " "))
  geneset[i] <- temp[1]
  matrix[,i]=0
  genelist[[i]] <- temp[-1]
  genelist[[i]]=as.vector(genelist[[i]])
  #find the location of each gene in the gene set
  generows[[i]]=as.vector(na.omit(as.vector(match(genelist[[i]],g$GENE))))
  if(length(generows[[i]]>0)){
  matrix[generows[[i]],i]=1
  }
}

#remove gene sets which have no genes in .genes.out file 
m.vars=apply(matrix,2,sum)
matrix=matrix[,which(m.vars>0)] 

#Scale each column of the matrix
matrix.scale=matrix
for (j in 1:ncol(matrix)){
  matrix.scale[,j]=scale(matrix[,j],center=T,scale=T) 
}

#Calculate the LD score of each gene set
ld=c()
for(i in 1:ncol(matrix.scale)){
  ld[i]=0
  for(j in 1:ncol(matrix.scale)){
    if(length(intersect(generows[[i]],generows[[j]]))!=0){
      rij=mean(matrix.scale[,i]*matrix.scale[,j])
      ld[i]=ld[i]+(rij^2)
    }
  }
}

#Calculate the chisq^2-Statistics of each gene set
x <- scale(g$ZSTAT,center=T,scale=T)
chi <- c()
for(i in 1:ncol(matrix.scale)){
  chi[i] <- (t(matrix.scale[,i])%*%x)^2/nrow(g)
}
f1=data.frame(ld=ld,chi=chi,stringsAsFactors = FALSE)

#Calculate the regression weights
w=c()
for(i in 1:length(ld)){
  w[i]=1/max(ld[i],1)
}
data=data.frame(chi=chi,ld=ld,w=w,stringsAsFactors = F)
data=data[!(is.na(data$chi)),]
data=data[!(is.na(data$ld)),]
data=data[-which(data$chi>150),]  #Remove extreme values

#Regression
lm.reg <- lm((data$chi)~(data$ld),weights = data$w)  #加权最小二乘法
lm.reg$coefficients