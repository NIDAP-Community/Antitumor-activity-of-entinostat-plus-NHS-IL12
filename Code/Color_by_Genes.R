# Load libraries
library(Seurat) # 3.1.5
library(ggplot2) # 3.3.3
library(tidyverse) # 1.3.0
library(gridExtra) # 2.3

# Color by Gene Code
SO_input <- SO_merge

ColorbyGene <- function(SO_input,sample_to_include,genes_to_view,nrow_plot){

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

SO = SO_input
print(SO)
doCiteSeq <- FALSE

samples = sample_to_include

if (length(samples) == 0) {
  samples = unique(SO@meta.data$sample_name)
}

if("active.ident" %in% slotNames(SO)){
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
} else {
  sample_name = as.factor(SO@meta.data$orig.ident)
  names(sample_name)=names(SO@active.ident)
  SO@active.ident <- as.factor(vector())
  SO@active.ident <- sample_name
  SO.sub = subset(SO, ident = samples)
}

gene = genes_to_view

nogene = gene[!gene %in% rownames(SO.sub$SCT@scale.data)]

if(!is.null(nogene)){
  print("Gene(s) missing from dataset:")
  print(nogene)
}

gene = gene[gene %in% rownames(SO.sub$SCT@scale.data)]

if(length(gene)>0){
  plotgene <- function(gene){
    gene.mat=SO.sub$SCT@scale.data[gene,]
    gene.quant=quantile(gene.mat[gene.mat>1],probs=c(.1,.5,.90))
    gene.mat[gene.mat>gene.quant[3]]=gene.quant[3]
    gene.mat[gene.mat<gene.quant[1]]=0
    
    Reduction = "umap"
    
    if (!(doCiteSeq)) {
      if(Reduction == "tsne"){
        p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, gene=gene.mat)
      }
      else if(Reduction == "umap"){
        p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, gene=gene.mat)
      }
      else{
        p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, gene=gene.mat)
      } #if CITEseq is chosen then:
    } else {
      if(Reduction == "tsne"){
        p1 <- DimPlot(SO.sub, reduction = "protein_tsne", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, gene=gene.mat)
      }
      else if(Reduction == "umap"){
        p1 <- DimPlot(SO.sub, reduction = "protein_umap", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, gene=gene.mat)
      }
      else{
        p1 <- DimPlot(SO.sub, reduction = "protein_pca", group.by = "ident")
        clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, gene=gene.mat)
      }
    }
    
    clusmat %>% dplyr::arrange(gene) -> clusmat
    g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      ggtitle(gene) +
      geom_point(aes(colour=gene),alpha=0.5,shape=16, size=1) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),legend.text=element_text(size=rel(0.5)) )+
      scale_color_gradient(limits = c(0, gene.quant[3]),low = "lightgrey", high = "red") +
      xlab("umap-1") + ylab("umap-2")  
    return(g)
  }
  
  grob <- lapply(seq_along(gene),function(x) plotgene(gene[x]))
  Final_figure <- gridExtra::grid.arrange(grobs=grob,nrow=nrow_plot,newpage=F)
  return(Final_figure)
}
else{
  print("No genes found in dataset")
  return(NULL)
}}

# Gene Arguments for ColorbyGene function:

# Main Figures
Figure_3b_genes <- c("Cd3e","Cd8a","Cd4","Foxp3","Ncr1","Ly6g","S100a9","G0s2","Csf3r","Itgam","Ly6c1","Itgax","H2-Ab1","Cst3","Siglec1","Adgre1","Arg1","Apoe","Mrc1","Cd38","Nos2")
Figure_5a_genes <- c("Ccl3","Ccl4","Ccl5","Cxcl9","Cxcl10")
Figure_5e_genes <- c("Ifng","Csf2","Ifnb1")
Figure_5g_genes <- c("Il1a","Il1b")
Figure_7e_genes <- c("Cd38","Nos2","Mrc1","Arg1")
Figure_7m_genes <- c("Ccl3","Ccl4","Ccl5","Cxcl9","Cxcl10")

# Supplementary Figures
Figure_S5b_genes <- c("Itgam","Adgre1","Cd68","Siglec1","Csf1r","Apoe",
                      "Mafb","F13a1","Lgals3","Cxcl9","Cxcl10","Arg1","Mmp13","Mmp12",
                      "Cd38","Nos2","Mrc1","Cd163","Cd40","Cd86")

Figure_S5c_genes <- c("Cd38","Nos2","Mrc1","Arg1")

# Figure Generation:

# Main Figures:
Figure_3b <- ColorbyGene(SO_input,sample_to_include = c("CD8dep","Combo","ENT","NHSIL12","PBS"),Figure_3b_genes,3)

Figure_5a_PBS <- ColorbyGene(SO_input,sample_to_include = "PBS",Figure_5a_genes,1)
Figure_5a_ENT <- ColorbyGene(SO_input,sample_to_include = "ENT",Figure_5a_genes,1)
Figure_5a_NHS_IL12 <- ColorbyGene(SO_input,sample_to_include = "NHSIL12",Figure_5a_genes,1)
Figure_5a_ENT_and_NHS_IL12 <- ColorbyGene(SO_input,sample_to_include = "Combo",Figure_5a_genes,1)

Figure_5e <- ColorbyGene(SO_input,sample_to_include = c("CD8dep","Combo","ENT","NHSIL12","PBS"),Figure_5e_genes,1)

Figure_5g <- ColorbyGene(SO_input,sample_to_include = c("CD8dep","Combo","ENT","NHSIL12","PBS"),Figure_5g_genes,1)

Figure_7e_CD8_Combo <- ColorbyGene(SO_input,sample_to_include = "Combo",Figure_7e_genes,1)
Figure_7e_CD8_Depl <- ColorbyGene(SO_input,sample_to_include = "CD8dep",Figure_7e_genes,1)

Figure_7m_CD8_Combo <- ColorbyGene(SO_input,sample_to_include = "Combo",Figure_7m_genes,1)
Figure_7m_CD8_Depl <- ColorbyGene(SO_input,sample_to_include = "CD8dep",Figure_7m_genes,1)

# Supplementary Figures:
Figure_S5b <- ColorbyGene(SO_input,sample_to_include = c("CD8dep","Combo","ENT","NHSIL12","PBS"),Figure_S5b_genes,4)

Figure_S5c_PBS <- ColorbyGene(SO_input,sample_to_include = "PBS",Figure_S5c_genes,2)
Figure_S5c_Combo <- ColorbyGene(SO_input,sample_to_include = "Combo",Figure_S5c_genes,2)

