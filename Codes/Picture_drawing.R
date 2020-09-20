#===================
#Boxplot
#===================
Orig <- read.table("gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)

re.Cortex <- read.table("Cortex.orig+re.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.Hipp <- read.table("Hippocampus.orig+re.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.Liver <- read.table("Liver.orig+re.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.SmaBow <- read.table("Small_Bowel.orig+re.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)

eqtl.Cortex <- read.table("Cortex.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
eqtl.Hipp <- read.table("Hippocampus.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
eqtl.Liver <- read.table("Liver.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
eqtl.SmaBow <- read.table("Small_Bowel.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)

re.eqtl.Cortex <- read.table("Cortex.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.eqtl.Hipp <- read.table("Hippocampus.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.eqtl.Liver <- read.table("Liver.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
re.eqtl.SmaBow <- read.table("Small_Bowel.orig+re+eqtl.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)

nMAGMA.Cortex <- read.table("Cortex.orig+re+eqtl+wgcna.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
nMAGMA.Hipp <- read.table("Hippocampus.orig+re+eqtl+wgcna.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
nMAGMA.Liver <- read.table("Liver.orig+re+eqtl+wgcna.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)
nMAGMA.SmaBow <- read.table("Small_Bowel.orig+re+eqtl+wgcna.gene_analysis.mean.genes.out",header = T,sep = "\t",check.names=F,stringsAsFactors = F)


x <- par(mar=c(5.1,5.1,4.1,2.1))
orig1 <- data.frame(orig[,5]);orig1$tissue="Cortex";
orig2 <- data.frame(orig[,5]);orig2$tissue="Hippocampus";
orig3 <- data.frame(orig[,5]);orig3$tissue="Liver";
orig4 <- data.frame(orig[,5]);orig4$tissue="Small_Bowel";
names(orig2) <- names(orig1)
names(orig3) <- names(orig1)
names(orig4) <- names(orig1)
data <- rbind(orig1,orig2)
data <- rbind(data,orig3)
data <- rbind(data,orig4)
names(data)=c("Nsnps","Tissue")
boxplot(Nsnps ~ Tissue , data = data,boxwex = 0.1,at = 1:4 - 0.7,col="#8DD3C7",ylim=c(0,1000),outline=F,xaxt = "n",xlim = c(0, 4.3))

re1 <- data.frame(re.Cortex[,5]);re1$tissue="Cortex"
re2 <- data.frame(re.Hipp[,5]);re2$tissue="Hippocampus"
re3 <- data.frame(re.Liver[,5]);re3$tissue="Liver"
re4 <- data.frame(re.SmaBow[,5]);re4$tissue="Small_Bowel"
names(re2) <- names(re1)
names(re3) <- names(re1)
names(re4) <- names(re1)
data <- rbind(re1,re2)
data <- rbind(data,re3)
data <- rbind(data,re4)
names(data)=c("Nsnps","Tissue")
boxplot(Nsnps ~ Tissue, data = data,
        boxwex = 0.1, at = 1:4-0.5,add = TRUE,col="#FFFFB3",outline=F)

eqtl1 <- data.frame(eqtl.Cortex[,5]);eqtl1$tissue="Cortex"
eqtl2 <- data.frame(eqtl.Hipp[,5]);eqtl2$tissue="Hippocampus"
eqtl3 <- data.frame(eqtl.Liver[,5]);eqtl3$tissue="Liver"
eqtl4 <- data.frame(eqtl.SmaBow[,5]);eqtl4$tissue="Small_Bowel"
names(eqtl2) <- names(eqtl1)
names(eqtl3) <- names(eqtl1)
names(eqtl4) <- names(eqtl1)
data <- rbind(eqtl1,eqtl2)
data <- rbind(data,eqtl3)
data <- rbind(data,eqtl4)
names(data)=c("Nsnps","Tissue")
boxplot(Nsnps ~ Tissue, data = data,
        boxwex = 0.1, at = 1:4 -0.3,add = TRUE,xaxt = "n",col="#BEBADA",outline=F)

re.eqtl1 <- data.frame(re.eqtl.Cortex[,5]);re.eqtl1$tissue="Cortex"
re.eqtl2 <- data.frame(re.eqtl.Hipp[,5]);re.eqtl2$tissue="Hippocampus"
re.eqtl3 <- data.frame(re.eqtl.Liver[,5]);re.eqtl3$tissue="Liver"
re.eqtl4 <- data.frame(re.eqtl.SmaBow[,5]);re.eqtl4$tissue="Small_Bowel"
names(re.eqtl2) <- names(re.eqtl1)
names(re.eqtl3) <- names(re.eqtl1)
names(re.eqtl4) <- names(re.eqtl1)
data <- rbind(re.eqtl1,re.eqtl2)
data <- rbind(data,re.eqtl3)
data <- rbind(data,re.eqtl4)
names(data) <- c("Nsnps","Tissue")
boxplot(Nsnps ~ Tissue, data = data,
        boxwex = 0.1, at = 1:4 -0.1,add = TRUE,xaxt = "n",col="#BEBADA",outline=F)

nMAGMA1 <- data.frame(nMAGMA.Cortex[,5]);nMAGMA1$tissue="Cortex"
nMAGMA2 <- data.frame(nMAGMA.Hipp[,5]);nMAGMA2$tissue="Hippocampus"
nMAGMA3 <- data.frame(nMAGMA.Liver[,5]);nMAGMA3$tissue="Liver"
nMAGMA4 <- data.frame(nMAGMA.SmaBow[,5]);nMAGMA4$tissue="Small_Bowel"
names(nMAGMA2) <- names(nMAGMA1)
names(nMAGMA3) <- names(nMAGMA1)
names(nMAGMA4) <- names(nMAGMA1)
data <- rbind(nMAGMA1,nMAGMA2)
data <- rbind(data,nMAGMA3)
data <- rbind(data,nMAGMA4)
names(data)=c("Nsnps","Tissue")
boxplot(Nsnps ~ Tissue, data = data,
        boxwex = 0.1, at = 1:4 +0.1,add = TRUE,xaxt = "n",col="#FB8072",outline=F)

#add legend
legend(3.4, 1040,col="black", c("MAGMA", "MAGMA+eqtl", "MAGMA+Hi-C","MAGMA+Hi-C+eQTL","nMAGMA"),
       fill = c("#6ED4C2", "#FFFF90","#B1ABDA","#FC4733", #5BB0D6),cex=0.7,text.font=2,bty="n") 


#===============
#Venn digram
#===============
library(VennDiagram)
x1 <- read.table("signifcant genes in MAGMA.txt", header = T,  sep = "\t",stringsAsFactors = F)
x2 <- read.table("signifcant genes in nMAGMA.txt", header = T,  sep = "\t",stringsAsFactors = F)
venn_data=list("MAGMA"=x1,"nMAGMA"=x2)
venn.diagram(venn_data,
             filename = 'Vennplot.tiff', 
             resolution = 800, 
             imagetype = 'tiff', 
             fontfamily ="serif", 
             main = paste('Liver'), 
             #main.col = 'black', 
             #main.fontface = 'bold', 
             #main.fontfamily = 'Helvetica', 
             #main.cex = 2, 
             main.just=c(3,3),
             col = c('#4292C6','red'), 
             fill = c('#FFFFB3','#BEBADA'), 
             scaled = TRUE,
             label.col = c('black','black','black'), 
             cex=0.7,
             cat.col = c('#4292C6','red'), 
             cat.cex = 0.9,
             cat.dist=0.05,
             cat.pos = c(-145,145)
)


#===========
#QQ plot
#===========
library("qqman")
f1 <- read.table("Cortex.merge.geneset_analysis.mean.gsa.out", header = T,  sep = "\t",stringsAsFactors = F)
qq(f1$P,ylim=c(0,10),col="blue",main="QQ-plot of dataset1 using nMAGMA",cex.main=1)