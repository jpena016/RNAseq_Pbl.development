library(DESeq2)
library(jsonlite)
library(ggplot2)
library(gplots)
library(dendextend)
library(dynamicTreeCut)
library(tidyverse)
library(cluster)
library(viridis)
library(reshape)
library(gridExtra)
library(bnlearn)
library(cowplot)
library(directlabels)
library(WGCNA)
library(data.table)
library(dplyr)
library(magrittr)
library(factoextra)
library(NbClust)
library(tximport)
library(DEGreport)
library(rhdf5)
library(pheatmap)

#read in abundance files 
samplefil<-read.table(file.path("samples.file.csv"), header=TRUE, sep = ",")
samplefil$FILE<-paste0(samplefil$Sample.ID,".", samplefil$Rep) #add column with file name construct
fileskall<-file.path("kallisto",samplefil$FILE, "abundance.h5")
names(fileskall) <- paste0(samplefil$FILE)
txi.kallisto <- tximport(fileskall, type = "kallisto", txOut = TRUE)

samples <- read.csv("samples.file.csv",header=TRUE)
exprnames <- do.call(paste,c(samples[c("Sample.ID","Rep")],sep="."))
#exprnames=samples$Sample.ID #added line
exprnames <- sub(".([123])$","r\\1",exprnames,perl=TRUE) #removed "." before r\\1
files <- file.path("kallisto",exprnames,"abundance.h5")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE, countsFromAbundance = )
head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance,file.path("reports","kallisto","kallisto_single.TPM.csv"))
write.csv(txi.kallisto$counts,file.path("reports","kallisto","kallisto_single.counts.csv"))

# DEseq2 analyses
geno = factor( samples$Sample.ID)
treatment = factor (samples$Stage)

sampleTable <- data.frame(condition = treatment,
                          genotype = geno)

rownames(sampleTable) = exprnames
dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable,~ condition)
dds<-estimateSizeFactors(dds)
normalized_counts<-counts(dds, normalized=TRUE)
dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable,~ condition)
dds <-DESeq(dds, test="LRT", reduced=~1)

#count transformation for downstream visualization
rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)

# Subset the LRT results to return genes with padj < 0.01
res<-results(dds)
res0<-as.data.frame(res)
res0<-res0[rowSums(res0[])>0,]
res0<-res0[complete.cases(res0), ]

sig_res_LRT <- res %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>% 
               as_tibble() %>% 
               filter(padj < .01) #padj <.01 or .05
 
# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
                
length(sigLRT_genes)

#subset rlog DEGs for clustering 
cluster_rlog<-rld_mat[sig_res_LRT$gene, ]

#hierarchical clustering of RLOG Norm counts
ZA_unmated.RL=rowMeans(cluster_rlog[,1:3])
ZA_early.RL=rowMeans(cluster_rlog[,4:6])
ZA_late.RL=rowMeans(cluster_rlog[,7:9])
z_ave.RL<-cbind(ZA_unmated.RL,ZA_early.RL,ZA_late.RL)

#scaledata.RL<-cluster_rlog[complete.cases(cluster_rlog),]
scaledata.RL<-t(scale(t(cluster_rlog))) #centers and scales daa
#scaledata.RL<-cluster_rlog
hr.RL=hclust(as.dist(1-cor(t(scaledata.RL),method = "pearson")),method="complete") #cluster rows by pearson correlation
hc.RL<-hclust(as.dist(1-cor(scaledata.RL, method="spearman")),method="complete") #cluster columns by spearman correlation

#with AVERAGES
scaledata_AVE.RL<-t(scale(t(z_ave.RL))) #centers and scales daa
#scaledata_AVE.RL<-z_ave.RL
#scaledata_AVE.RL<-z_ave.RL[complete.cases(z_ave.RL),]
hr_AVE.RL=hclust(as.dist(1-cor(t(scaledata_AVE.RL),method = "pearson")),method="complete") #cluster rows by pearson correlation
hc_AVE.RL<-hclust(as.dist(1-cor(scaledata_AVE.RL, method="spearman")),method="complete") #cluster columns by spearman correlation



png("plots/rlog.clustering.replicates.png", width = 230, height = 345, units='mm', res = 300)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
heatmap.2(z2_labeled.RL, Rowv=as.dendrogram(hr.RL), Colv=reorder(as.dendrogram(hc.RL), 1:9), col=viridis(100), scale="row", labRow=FALSE, trace="none", srtCol = 0, cexCol=1,keysize=0.6, key.par = list(cex=0.5))
dev.off()

png("plots/rlog.clustering.AVE.png", width = 230, height = 345, units='mm', res = 300)
par(
  mar      = c(5, 5, 5, 5),
              #b,L,T,r
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2
)
heatmap.2(z_ave_labeled.RL, Rowv=as.dendrogram(hr_AVE.RL), Colv=as.dendrogram(hc_AVE.RL), col=viridis(100), scale="row", labRow=FALSE, trace="none", srtCol = 45, cexCol=1,
          #lhei=c(.1,.2), lwid=c(.2,.35), 
          keysize=0.6, key.par = list(cex=0.5))
dev.off()

