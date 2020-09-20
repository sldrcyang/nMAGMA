setwd("/your path/")
#The first file
snp_to_gene1 <- readLines("file1.genes.annot")
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene1)){
  temp <- unlist(strsplit(snp_to_gene1[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot1 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp1 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#The second file
snp_to_gene2 <- readLines("file2.genes.annot")
library(dplyr)
gene <- c()
bp <- c()
snp <- c()
for(i in 1:length(snp_to_gene2)){
  temp <- unlist(strsplit(snp_to_gene2[i], "\t"))
  gene[i] <- temp[1]
  bp[i] <- temp[2]
  snp[i] <- paste(temp[-2:-1],collapse="\t")
}
annot2 <- data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
gene_bp2 <- data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#Merge
annot_comb <- rbind(annot1,annot2)
#If the annotate files contains some non-coding genes, and you want retain only coding genes
gene_bp=read.table("protein_coding_gene.loc", header = F,  sep = "\t",stringsAsFactors = F)
library(tidyr)
gene_bp=tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
gene_bp=tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
gene_bp=gene_bp[,c(1,2)]
#If you have retained coding genes when you generating annotation files
#gene_bp <- rbind(gene_bp1,gene_bp2)
#gene_bp <- gene_bp[!duplicated(gene_bp$gene),]

annot <- annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
annot <- merge(gene_bp,annot,by="gene",all=F)
write.table(annot,"Intermediate.genes.annot",row.names = F,col.names=F,sep="\t",quote=F)

#Remove duplicate SNPs on the same gene
data <- readLines("Intermediate.genes.annot")
data <- as.list(data)
temp <- list()
for(i in 1:length(data)){
  temp[[i]] <- unlist(strsplit(data[[i]], "\t"))
  temp[[i]] <- as.vector(temp[[i]])
  temp[[i]] <- unique(temp[[i]])
}
file.remove('Intermediate.genes.annot')
for(i in 1:length(temp)){
  cat(temp[[i]], file = "Merge.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  cat("\n", file = "Merge.genes.annot", sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
}
