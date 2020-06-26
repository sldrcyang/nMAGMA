setwd("/your path/")
#The first file
snp_to_gene1=readLines("file1.genes.annot")
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
gene_bp1=data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#The second file
snp_to_gene2=readLines("file2.genes.annot")
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
gene_bp2=data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)

#Merge
annot_comb=rbind(annot1,annot2)
gene_bp=rbind(gene_bp1,gene_bp2)
gene_bp=gene_bp[!duplicated(gene_bp$gene),]
annot <- annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
annot=merge(gene_bp,annot,by="gene")
write.table(annot,"Merge.genes.annot",row.names = F,col.names=F,sep="\t",quote=F)