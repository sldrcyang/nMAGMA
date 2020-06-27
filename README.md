# T-MAGMA
T-MAGMA is an approach to intergrate Hi-C, eQTL and TOM matrix of gene expression calculated by WGCNA into original physical location based MGAMA annotation.
## List of files
* __RE_to_genes.zip__: *3DIV_DLPFC_RE_to_genes.txt,3DIV_Hippocampus_RE_to_genes.txt,3DIV_Liver_RE_to_genes.txt,3DIV_Small_Bowell_RE_to_genes.txt.* This file contains RE-gene regulatory relationships obtained from Hi-C.
* __Signifgenepairs.zip__: *Cortex_signifgenepairs.txt,Hippocampus_signifgenepairs.txt,Liver_signifgenepairs.txt,Small_Bowel_signifgenepairs.txt.* This file contains significant gene pairs with TOM >= 0.15.
* __Hi-C_annot__: *DLPFC_assig_throu_re_dataset1.genes.annot, DLPFC_assig_throu_re_dataset2.genes.annot, Hippocampus_assig_throu_re_dataset1.genes.annot, Hippocampus_assig_throu_re_dataset2.genes.annot, Liver_assig_throu_re_dataset1.genes.annot, Liver_assig_throu_re_dataset2.genes.annot, Small_Bowel_assig_throu_re_dataset1.genes.annot, Small_Bowel_assig_throu_re_dataset2.genes.annot.* This folder contains .genes.annot files of two EUR GWAS summary datasets we used through Hi-C. Snps were first mapped to REs and then mapped to genes via the RE-gene relationships.
* __eqtl_annot__: *Cortex_eqtl.genes.annot, Hippocampus_eqtl.genes.annot, Liver_eqtl.genes.annot, Small_Bowel_eqtl.genes.annot.* This folder contains .genes.annot files through eQTL. Snps were directly mapped to their eGenes.
* __TOM_annot__: *Cortex_TOM.zip, Hippocampus_TOM.zip, Liver_TOM.zip, Small_Bowel_TOM.zip.* This folder contains .genes.annot files through assigning eQTLs to genes highly connected with its target gene (TOM >= 0.15).
## Scripts
* __Calcu_TOM_extract_signifgenepairs.R__: This script can calaulate a TOM matrix using R package 'WGCNA' and extract siginificant gene regulatory pairs with TOM >=0.15. Your can also define a threshold you prefer.
* __Merge_annot_files.R__: Your can use this script to merge to .genes.annot file in order to expand the annotation. This result file can directly used as input for gene-analysis in MGAMA. By merging the original MAGMA annotation file and three kind of .annot files we provides, you can get the final T-MAGMA annotation file used in our paper.
* __Extended LDSC to GSA.R__: This is the LD Score regression extended to gene-set analysis. First a gene-geneset matrix is created and then the LD score and chisq statistics of each gene set are calculated. Your can detect from the regression intercept that whether the inflation of statistics in gene-set analysis are caused by confounding factors or the complex pathogenesis of the disease. 
## Reference