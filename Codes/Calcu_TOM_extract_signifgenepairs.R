#====================
#Data preparation
#====================
setwd('/your path/')
expro <- read.table("your expression data file name", header = T,  sep = "\t",check.names=F)
options(stringsAsFactors = FALSE)
dim(expro)

#Select genes in the top fifth of the variance across samples
m.vars <- apply(expro,1,var) 
expro.upper <- expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.2))[5]),]
dim(expro.upper)
write.table(expro.upper,file="geneInput_variancetop0.2.txt",sep='\t ',quote=F,row.names=T)

#Transpose the expression matrix with rows are genes and columns are samples
datExpr0 <- as.data.frame(t(expro.upper))
# Load the WGCNA package 
library(WGCNA)
#Qualify the expression matrix
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#======================================
# Sample clustering to detect outliers
#======================================
sampleTree <- hclust(dist(datExpr0), method = "average"); 
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches 
# The user should change the dimensions if the window is too large or too small. 
sizeGrWindow(12,9) 
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9); 
par(cex = 0.45); par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) 

# Plot a line to show the cut 
abline(h = 110000, col = "red"); 
# Determine cluster under the line 
clust <- cutreeStatic(sampleTree, cutHeight = 110000, minSize = 10) 
table(clust) 
# clust 1 contains the samples we want to keep. 
keepSamples <- (clust==1) 
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr) 
nSamples <- nrow(datExpr) 

#=========================
# Network constrution
#=========================
# The following setting is important, do not omit. 
options(stringsAsFactors = FALSE); 
# Allow multi-threading within WGCNA. At present this call is necessary. 
# Any error here may be ignored but you may want to update WGCNA if you see one. 
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above. 
enableWGCNAThreads() 
# Load the data saved in the first part 
#lnames = load(file = "G-01-dataInput.RData"); 
#The variable lnames contains the names of loaded variables. 
#lnames 
#==============================================
#  1) choose the most suitable soft threshold
#==============================================
# Choose a set of soft-thresholding powers 
powers <- c(c(1:10), seq(from = 12, to=24, by=2)) 
# Call the network topology analysis function 
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results: 
sizeGrWindow(9, 5) 
par(mfrow = c(1,2)); 
cex1 = 0.9; 
# Scale-free topology fit index as a function of the softthresholding power 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")); text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2],
       labels=powers,cex=cex1,col="red"); # this line corresponds to using an R^2 cut-off of h 
abline(h=0.85,col="red")
sft$powerEstimate #your suitable soft-thresholding

#Mean connectivity
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],     
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",     
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# here we define the adjacency matrix using soft thresholding with beta
beta <- sft$powerEstimate
ADJ1 <- abs(cor(datExpr,use="p"))^beta
# When you have relatively few genes (<5000) use the following code 
k <- as.vector(apply(ADJ1,2,sum, na.rm=T))#pick one of the two methods
# When you have a lot of genes use the following code 
k <- softConnectivity(datE=datExpr,power=beta) 
# Plot a histogram of k and a scale free topology plot 
sizeGrWindow(10,5) 
par(mfrow=c(1,2)) 
hist(k) 
scaleFreePlot(k, main="Check scale free topology\n") 
#=======================================
#  2) Calculate the adjacency matrix
#=======================================
softPower <- beta; 
adjacency <- adjacency(datExpr, power = softPower) 
#======================================================
#  3) Calculate the final topological overlap matrix
#======================================================
TOM <- TOMsimilarity(adjacency);
dimnames(TOM) <- dimnames(adjacency)
TOM <- as.data.frame(TOM)

#==========================================================
#Extract strongly interacted gene pairs with TOM >= 0.15
#==========================================================
result <- data.frame(stringsAsFactors = FALSE) 
for (i in 1:(nrow(TOM)-1)){
  for (j in (i+1):nrow(TOM)){
    if (TOM[i,j] >= 0.15){
      tmp.result <- data.frame(rownames(TOM[i,]),rownames(TOM[j,]),TOM[i,j],i,j,stringsAsFactors = FALSE)
      names(tmp.result) <- c('gene1','gene2','TOM','i','j')
      result <- rbind(result, tmp.result)
      names(result) <- c('gene1','gene2','TOM','i','j')
      as.numeric(as.character(result[nrow(result),3])) 
      as.numeric(as.character(result[nrow(result),4]))
      as.numeric(as.character(result[nrow(result),5]))
    }
  }
}
write.table(result, "Tissue.0.15.TOM.signifgenepairs.txt",row.names = F,col.names=T,sep="\t",quote=F)
